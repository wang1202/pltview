/*
 * pltview.c - Fast AMReX plotfile viewer in C
 * Similar to ncview, using X11/Athena Widgets for GUI
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Box.h>
#include <X11/Xaw/Scrollbar.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/Simple.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/AsciiText.h>
#include <X11/Xaw/Viewport.h>

#define MAX_VARS 128
#define MAX_BOXES 1024
#define MAX_PATH 512
#define MAX_LINE 1024
#define MAX_TIMESTEPS 1024
#define MAX_LEVELS 10
#define VIDEO_FRAME_WIDTH 1000
#define VIDEO_FRAME_HEIGHT 700
#define VIDEO_CACHE_LIMIT ((size_t)2 * 1024 * 1024 * 1024)

/* Data structures */
typedef struct {
    int lo[3];
    int hi[3];
    char filename[64];
    long long offset;
} Box;

/* Per-level data storage for multi-level overlay rendering */
typedef struct {
    int grid_dims[3];       /* Grid dimensions for this level */
    int level_lo[3];        /* Lower index bounds in level's coordinates */
    int level_hi[3];        /* Upper index bounds in level's coordinates */
    Box boxes[MAX_BOXES];   /* Box definitions for this level */
    int n_boxes;            /* Number of boxes at this level */
    double *data;           /* Variable data for this level */
    double *z_phys_data;    /* Cached physical Z coordinate for terrain-following grids */
    int loaded;             /* Flag: 1 if data is loaded, 0 otherwise */
    int z_phys_loaded;      /* Flag: 1 if z_phys_data is loaded */
} LevelData;

typedef struct {
    char plotfile_dir[MAX_PATH];
    char variables[MAX_VARS][64];
    int n_vars;
    int ndim;
    double time;
    int grid_dims[3];
    int level_lo[3];    /* Current level's lower index bounds */
    int level_hi[3];    /* Current level's upper index bounds */
    Box boxes[MAX_BOXES];
    int n_boxes;
    double *data;  /* Current variable data */
    int current_var;
    int slice_axis;
    int slice_idx;
    int colormap;  /* 0=viridis, 1=jet, 2=turbo, 3=plasma */
    int current_level;  /* Current AMR level */
    int n_levels;       /* Number of AMR levels */
    double prob_lo[3];  /* Domain lower bounds */
    double prob_hi[3];  /* Domain upper bounds */
    /* Multi-level overlay data */
    LevelData levels[MAX_LEVELS];  /* Per-level data for overlay rendering */
    int ref_ratio[MAX_LEVELS];     /* Refinement ratio between levels */
    int overlay_mode;              /* 0=single level, 1=overlay all levels */
    int map_mode;                  /* 0=normal view, 1=longitude-latitude map view */
    int use_z_phys;                /* Use 3D z_phys for X/Y-slice vertical coordinates */
    double *z_phys_data;           /* Cached z_phys data for the current level */
    int z_phys_level;              /* Level represented by z_phys_data */
} PlotfileData;

typedef struct {
    unsigned char *rgb;
    double time;
    int timestep_number;
} VideoFrame;

typedef struct VideoState {
    Widget shell, canvas_widget, status_label, frame_label, scrubber;
    Widget fps_text, play_button, stop_button, save_button;
    Window canvas;
    Pixmap frame_pixmap;
    GC gc;
    XImage *display_image;
    unsigned char *display_rgb;
    VideoFrame *frames;
    int n_frames, current_frame, frame_width, frame_height;
    double global_vmin, global_vmax, fps;
    int loading, playing, closing, load_phase, load_index;
    XtWorkProcId work_id;
    XtIntervalId timer_id;
    char variable_name[64];
    int slice_axis, requested_level, overlay_mode, map_mode, colormap, scale_mode;
    double slice_position;
    Widget save_shell, save_text;
} VideoState;

/* Colormap */
typedef struct {
    unsigned char r, g, b;
} RGB;

/* Plot data for popup window */
typedef struct {
    double *data;
    double *x_values;  /* X-axis coordinate values */
    int n_points;
    double vmin, vmax;
    double xmin, xmax;  /* X-axis range */
    char title[128];
    char xlabel[64];    /* X-axis label (or Y-axis for horizontal plots) */
    char vlabel[64];    /* Value axis label (X-axis for horizontal plots) */
} PlotData;

/* Popup window data */
typedef struct {
    Widget shell;
    PlotData *plot_data_array[3];
} PopupData;

#define MAX_SDM_VARS 32
#define SDM_SUBDIR "super_droplets_moisture"

/* Y-axis metric types for SDM histogram */
#define SDM_METRIC_PARTICLE_COUNT   0  /* Sum of multiplicity per bin */
#define SDM_METRIC_SD_COUNT         1  /* Number of super droplets per bin */
#define SDM_METRIC_CONCENTRATION    2  /* Particle count / domain volume */
#define SDM_METRIC_MASS             3  /* Sum of (mass * multiplicity) per bin */
#define SDM_METRIC_MEAN_MULT        4  /* Mean multiplicity per bin */
#define SDM_N_METRICS               5

typedef struct {
    int n_particles;
    int n_real_comps;   /* Number of real components (from Header, excluding x,y,z) */
    int n_int_comps;    /* Number of int components (from Header, excluding id,cpu) */
    char real_comp_names[MAX_SDM_VARS][64];
    char int_comp_names[MAX_SDM_VARS][64];
    int ndim;
    double *radius;        /* Extracted radius array */
    double *multiplicity;  /* Extracted multiplicity array */
    double *mass;          /* Extracted particle_mass array */
    int radius_idx;        /* Index of "radius" in real comp names */
    int mult_idx;          /* Index of "multiplicity" */
    int mass_idx;          /* Index of "particle_mass" */
    double domain_volume;  /* For number concentration */
    int current_metric;    /* Current y-axis metric (SDM_METRIC_*) */
    int log_x;              /* 0=linear, 1=log10 x-axis */
    int log_y;              /* 0=linear, 1=log10 y-axis */
    int log_bin;            /* 0=linear bin width, 1=log-spaced bins */
    double cutoff_radius;   /* Cutoff in um (0 = no cutoff) */
    double custom_bin_width; /* Custom bin width in um (0 = auto/Sturges) */
    double xlim_min;        /* X-axis min in um (0 = auto) */
    double xlim_max;        /* X-axis max in um (0 = auto) */
    int pdf_mode;           /* PDF mode: normalize by sum and bin width */
    /* Per-grid info from particle Header */
    int n_grids;
    int grid_file_num[MAX_BOXES];
    int grid_count[MAX_BOXES];
    long grid_offset[MAX_BOXES];
} ParticleData;

/* SBM (Spectral Bin Microphysics) mode definitions */
#define MAX_SBM_BINS 64
#define SBM_BIN_INFO_FILE "bin_info.txt"

/* SBM metric types */
#define SBM_METRIC_QC_MASS     0  /* Cloud mass per bin */
#define SBM_METRIC_QI_MASS     1  /* Ice mass per bin */
#define SBM_METRIC_TOTAL_MASS  2  /* Cloud + Ice mass per bin */
#define SBM_METRIC_QC_NUM      3  /* Cloud number per bin */
#define SBM_METRIC_QI_NUM      4  /* Ice number per bin */
#define SBM_METRIC_TOTAL_NUM   5  /* Cloud + Ice number per bin */
#define SBM_N_METRICS          6

typedef struct {
    int n_bins;                         /* Number of bins (from bin_info.txt) */
    double bin_radius_um[MAX_SBM_BINS]; /* Radius in micrometers per bin */
    double bin_diameter_um[MAX_SBM_BINS]; /* Diameter in micrometers per bin */
    double bin_values[MAX_SBM_BINS];    /* Summed values per bin (current metric) */
    double domain_volume;               /* For concentration */
    int current_metric;                 /* Current y-axis metric */
    int log_x;                          /* Log x-axis toggle */
    int log_y;                          /* Log y-axis toggle */
    int pdf_mode;                       /* PDF mode: normalize by sum and bin width */
    double xlim_min;                    /* X-axis min in um (0 = auto) */
    double xlim_max;                    /* X-axis max in um (0 = auto) */
    char plotfile_dir[MAX_PATH];
} SBMData;

/* 1D Profile mode definitions */
#define PROFILE_FILE_SURF    0
#define PROFILE_FILE_MEAN    1
#define PROFILE_FILE_FLUX    2
#define PROFILE_FILE_SUBGRID 3
#define N_PROFILE_FILES      4
#define MAX_PROFILE_ROWS     500000
#define MAX_PROFILE_COLS     32

typedef struct {
    char filename[256];   /* basename of the file (surf/mean/flux/subgrid) */
    int has_z;            /* 0 for surf (time-series), 1 for mean/flux/subgrid */
    int ncols;            /* total number of data columns */
    char col_names[MAX_PROFILE_COLS][32];
    double *data;         /* [nrows * ncols] row-major */
    int nrows;
    /* For has_z=1: indexed by timestep */
    double *times;        /* unique sorted time values [ntimes] */
    int ntimes;
    int nz;               /* z levels per timestep */
    double z_min, z_max;
} ProfileFile;

typedef struct {
    char dir[MAX_PATH];
    ProfileFile files[N_PROFILE_FILES];
    int loaded[N_PROFILE_FILES];   /* 1 if file exists and was read */
} ProfileData;

/* X11 globals */
Display *display;
Widget toplevel, form, canvas_widget, var_box, info_label;
Widget var_viewport, var_scrollbar;
Widget axis_box, nav_box, colorbar_widget, layer_label;
Window canvas, colorbar;
GC gc, text_gc, colorbar_gc;
XImage *ximage;
int screen;
unsigned long *pixel_data;
size_t pixel_data_size = 0;
int canvas_width = 800;
int canvas_height = 600;
Pixmap pixmap, colorbar_pixmap;
XFontStruct *font;
double current_vmin = 0, current_vmax = 1;

/* Current slice rendering info for mouse interaction */
double *current_slice_data = NULL;
double *current_z_phys_slice = NULL;
int slice_width = 0, slice_height = 0;
int render_offset_x = 0, render_offset_y = 0;
int render_width = 0, render_height = 0;
double render_phys_ymin = 0.0, render_phys_ymax = 1.0;
int render_uses_z_phys = 0;
int *map_hover_lookup = NULL;       /* Canvas pixel -> base-slice cell index */
size_t map_hover_lookup_size = 0;
int map_hover_lookup_active = 0;
char hover_value_text[256] = "";

/* Zoom and scroll state */
double zoom_level = 1.0;        /* 1.0 = no zoom */
int zoom_scroll_x = 0;          /* Pixel offset into zoomed view (horizontal) */
int zoom_scroll_y = 0;          /* Pixel offset into zoomed view (vertical) */
int zoom_dragging = 0;          /* 1 = mouse drag in progress */
int zoom_drag_start_x = 0;      /* Mouse X at drag start */
int zoom_drag_start_y = 0;      /* Mouse Y at drag start */
int zoom_drag_scroll_x0 = 0;    /* scroll_x at drag start */
int zoom_drag_scroll_y0 = 0;    /* scroll_y at drag start */
int vis_area_x = 0, vis_area_y = 0, vis_area_w = 0, vis_area_h = 0;  /* Visible clip region on screen */
int initial_focus_set = 0;  /* Flag for setting keyboard focus on first expose */
int dialog_active = 0;  /* Flag to track when a dialog is open */
Widget active_text_widget = NULL;  /* Text widget in active dialog */

/* Custom colorbar range */
int use_custom_range = 0;  /* Flag for using custom min/max */
double custom_vmin = 0.0;
double custom_vmax = 1.0;

/* Scale mode for colorbar: 0=Linear, 1=Log(+), 2=Log(-) */
int scale_mode = 0;

/* Multi-timestep support */
char *timestep_paths[MAX_TIMESTEPS];  /* Array of plotfile paths */
int timestep_numbers[MAX_TIMESTEPS];   /* Numerical values for sorting */
int timestep_levels[MAX_TIMESTEPS];    /* Number of levels at each timestep */
int n_timesteps = 0;                   /* Number of timesteps found */
int current_timestep = 0;              /* Current timestep index */
int max_levels_all_timesteps = 1;      /* Max levels across all timesteps */
Widget time_label = NULL;              /* Time step display label */
static VideoState *active_video = NULL;

/* Data structure for histogram expose handler (forward declaration for SDM) */
typedef struct {
    double *bin_counts;
    double *bin_centers;
    double *bin_edges;    /* Optional: n_bins+1 edges for non-uniform bins (NULL if uniform) */
    int n_bins;
    double count_max;
    double bin_min, bin_max;
    double total;         /* Total sum of all bin values (before PDF normalization) */
    char title[128];
    char xlabel[64];
    double mean, std, skewness;
    double kurtosis;
} HistogramData;

/* SDM mode globals */
int sdm_mode = 0;                /* Flag for SDM mode */
ParticleData *global_pd = NULL;  /* Global particle data pointer */
Widget sdm_canvas_widget = NULL;
Widget sdm_metric_buttons[SDM_N_METRICS];
Widget sdm_info_label = NULL;
Window sdm_canvas = 0;
HistogramData *sdm_hist_data = NULL;  /* Persistent histogram data for SDM canvas */
Widget sdm_settings_text_cutoff = NULL;
Widget sdm_settings_text_binwidth = NULL;
Widget sdm_settings_text_xlim_min = NULL;
Widget sdm_settings_text_xlim_max = NULL;
int sdm_dialog_active = 0;
Widget sdm_active_text_widget = NULL;
int sdm_active_field = 0;  /* 0=cutoff, 1=binwidth, 2=xlim_min, 3=xlim_max */
Widget sdm_dialog_shell = NULL;
Widget overlay_button = NULL;  /* Overlay toggle button */

/* SBM mode globals */
int sbm_mode = 0;                /* Flag for SBM mode */
SBMData *global_sbm = NULL;      /* Global SBM data pointer */
Widget sbm_canvas_widget = NULL;
Widget sbm_metric_buttons[SBM_N_METRICS];
Widget sbm_info_label = NULL;
Window sbm_canvas = 0;
HistogramData *sbm_hist_data = NULL;  /* Persistent histogram data for SBM canvas */
Widget sbm_settings_text_xlim_min = NULL;
Widget sbm_settings_text_xlim_max = NULL;
int sbm_dialog_active = 0;
Widget sbm_active_text_widget = NULL;
int sbm_active_field = 0;  /* 0=xlim_min, 1=xlim_max */
Widget sbm_dialog_shell = NULL;
/* 1D Profile mode globals */
int profile_mode = 0;
ProfileData *global_profile = NULL;
Widget profile_canvas_widget = NULL;
Window profile_canvas = 0;
Widget profile_info_label = NULL;
int profile_current_file = PROFILE_FILE_MEAN;   /* which file to show */
int profile_current_col = 2;                     /* which column (0=time,1=z,2=first data) */
int profile_current_time_idx = 0;
int profile_log_x = 0;
Widget profile_file_buttons[N_PROFILE_FILES];
Widget *profile_var_buttons = NULL;
int profile_n_var_buttons = 0;
Widget profile_var_box_widget = NULL;   /* viewport's inner box for var buttons */
Widget profile_var_viewport = NULL;     /* scrollable viewport for var buttons */
int profile_contour_mode = 0;           /* 0=profile line, 1=time-height contour */

Widget map_dialog_shell = NULL;
Widget map_unavailable_shell = NULL;
Widget map_unavailable_label = NULL;
int map_color_option = 0;  /* 0=black, 1=red, 2=gray, 3=white */
unsigned long map_color_pixel = 0;
int map_coastlines_enabled = 1;
double map_last_lon_min = 0.0, map_last_lon_max = 0.0, map_last_lat_min = 0.0, map_last_lat_max = 0.0;
int map_has_bounds = 0;
int map_auto_detected = 0;
static char map_layers_dir[MAX_PATH] = "";
static int map_layers_available = 0;

#define MAX_COASTLINES 64
typedef struct {
    char filename[MAX_PATH];
    char label[128];
    int enabled;
    int bbox_loaded;
    double lon_min, lon_max, lat_min, lat_max;
    Widget button;
} CoastlineEntry;

static CoastlineEntry coastlines[MAX_COASTLINES];
static int n_coastlines = 0;

/* Function prototypes */
int detect_levels(PlotfileData *pf);
int detect_levels_for_path(const char *plotfile_dir);
void show_level_warning(int level);
int read_header(PlotfileData *pf);
int read_cell_h(PlotfileData *pf);
int read_variable_data(PlotfileData *pf, int var_idx);
void extract_slice(PlotfileData *pf, double *slice, int axis, int idx);
void extract_slice_level(LevelData *ld, double *slice, int axis, int idx);
/* Multi-level overlay functions */
int read_cell_h_level(PlotfileData *pf, int level);
int read_variable_data_level(PlotfileData *pf, int var_idx, int level);
int load_all_levels(PlotfileData *pf, int var_idx);
void free_all_levels(PlotfileData *pf);
void apply_colormap(double *data, int width, int height, 
                   unsigned long *pixels, double vmin, double vmax, int cmap_type);
RGB viridis_colormap(double t);
RGB jet_colormap(double t);
RGB turbo_colormap(double t);
RGB plasma_colormap(double t);
RGB hot_colormap(double t);
RGB cool_colormap(double t);
RGB gray_colormap(double t);
RGB magma_colormap(double t);
RGB get_colormap_rgb(double t, int cmap_type);
const char *get_variable_unit(const char *varname);
void draw_colorbar(double vmin, double vmax, int cmap_type, const char *varname);
void cmap_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void colormap_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void colorbar_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void colorbar_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void init_gui(PlotfileData *pf, int argc, char **argv);
void render_slice(PlotfileData *pf);
void update_info_label(PlotfileData *pf);
void var_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void axis_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void z_phys_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void update_z_phys_button(PlotfileData *pf);
void map_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_map_settings_dialog(PlotfileData *pf);
void show_map_unavailable_dialog(const char *message);
void map_color_callback(Widget w, XtPointer client_data, XtPointer call_data);
void map_remove_callback(Widget w, XtPointer client_data, XtPointer call_data);
void map_unavailable_ok_callback(Widget w, XtPointer client_data, XtPointer call_data);
void map_coastline_toggle_callback(Widget w, XtPointer client_data, XtPointer call_data);
void render_map_overlay(PlotfileData *pf, double lon_min, double lon_max, double lat_min, double lat_max);
void level_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void overlay_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void nav_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void jump_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void range_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void profile_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_slice_statistics(PlotfileData *pf);
void distribution_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_distribution(PlotfileData *pf);
void distrib_mode_callback(Widget w, XtPointer client_data, XtPointer call_data);
void update_distribution_histogram(int mode);
void draw_histogram(Display *dpy, Window win, GC plot_gc, double *bin_counts, double *bin_centers,
                    int n_bins, int width, int height, double count_max,
                    double bin_min, double bin_max, const char *title, const char *xlabel,
                    double mean, double std, double skewness);
void fft2d_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_2dfft(PlotfileData *pf);
void quiver_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_quiver_dialog(PlotfileData *pf);
int find_variable_index(PlotfileData *pf, const char *name);
int find_velocity_component(PlotfileData *pf, const char *primary, char fallback_char);
void get_default_quiver_components(PlotfileData *pf, char *x_comp, char *y_comp);
void quiver_apply_callback(Widget w, XtPointer client_data, XtPointer call_data);
void quiver_close_callback(Widget w, XtPointer client_data, XtPointer call_data);
void quiver_remove_callback(Widget w, XtPointer client_data, XtPointer call_data);
void quiver_density_callback(Widget w, XtPointer client_data, XtPointer call_data);
void quiver_scale_callback(Widget w, XtPointer client_data, XtPointer call_data);
void quiver_color_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_variable_selector(Widget w, XtPointer client_data, XtPointer call_data);
void variable_select_callback(Widget w, XtPointer client_data, XtPointer call_data);
void variable_selector_close_callback(Widget w, XtPointer client_data, XtPointer call_data);
void render_quiver_overlay(PlotfileData *pf);
void draw_arrow(Display *dpy, Drawable win, GC graphics_gc, int x1, int y1, int x2, int y2);
void extract_slice_from_data(double *data, PlotfileData *pf, double *slice, int axis, int idx);
int ensure_z_phys_data(PlotfileData *pf);
int ensure_z_phys_level_data(PlotfileData *pf, int level);
void free_z_phys_cache(PlotfileData *pf);
double z_phys_corner(const double *z_values, int width, int height,
                     int x_boundary, int z_boundary);
void draw_z_phys_cells(const double *z_values, const unsigned long *pixels,
                       const unsigned char *mask, int width, int height,
                       double data_xmin, double data_xmax,
                       double view_xmin, double view_xmax,
                       double view_ymin, double view_ymax,
                       int offset_x, int offset_y, int draw_width, int draw_height);
int canvas_to_data_indices(int mouse_x, int mouse_y, int *data_x, int *data_y);
int prepare_map_hover_lookup(void);
void record_map_hover_rect(int x, int y, int width, int height, int cell_idx);
void update_layer_label(PlotfileData *pf);
void canvas_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void canvas_motion_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void canvas_button_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void canvas_button_release_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void zoom_in_callback(Widget w, XtPointer client_data, XtPointer call_data);
void zoom_out_callback(Widget w, XtPointer client_data, XtPointer call_data);
void zoom_reset(void);
void clamp_zoom_scroll(void);
void show_line_profiles(PlotfileData *pf, int data_x, int data_y);
void cleanup(PlotfileData *pf);
int scan_timesteps(const char *base_dir, const char *prefix);
int scan_sdm_timesteps(const char *base_dir, const char *prefix);
void switch_timestep(PlotfileData *pf, int new_timestep);
void time_nav_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void update_time_label(void);
void time_jump_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void time_series_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_time_series(PlotfileData *pf);
void video_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
static Boolean video_load_workproc(XtPointer client_data);
static void video_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_play_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_stop_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_scrub_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_timer_callback(XtPointer client_data, XtIntervalId *id);
static void video_close_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_save_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_save_gif_callback(Widget w, XtPointer client_data, XtPointer call_data);
static void video_save_mp4_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_time_height_contour(PlotfileData *pf);
void profile_fmt_val(char *buf, int bufsz, double v);
static void time_height_contour_callback(Widget w, XtPointer client_data, XtPointer call_data);

/* SDM functions */
int read_sdm_header(ParticleData *pd, const char *plotfile_dir);
int read_sdm_data(ParticleData *pd, const char *plotfile_dir);
double compute_domain_volume(const char *plotfile_dir);
void compute_sdm_histogram(ParticleData *pd, HistogramData *hist);
void init_sdm_gui(ParticleData *pd, const char *plotfile_dir, int argc, char **argv);
void render_sdm_histogram(ParticleData *pd);
void draw_sdm_histogram(Display *dpy, Window win, GC plot_gc, HistogramData *hist,
                         int width, int height, int log_x, int log_y, const char *ylabel);
void sdm_switch_timestep(ParticleData *pd, int new_timestep);
void update_sdm_info_label(ParticleData *pd, const char *plotfile_dir);
void sdm_logx_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_logy_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_logbin_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_pdf_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_settings_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_settings_apply_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_settings_close_callback(Widget w, XtPointer client_data, XtPointer call_data);
void var_viewport_wheel_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void var_scrollbar_scroll_proc(Widget w, XtPointer client_data, XtPointer call_data);
void var_scrollbar_jump_proc(Widget w, XtPointer client_data, XtPointer call_data);
void var_viewport_configure_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);

/* SBM functions */
int read_sbm_bin_info(SBMData *sbm, const char *plotfile_dir);
int compute_sbm_values(SBMData *sbm, const char *plotfile_dir);
void compute_sbm_histogram(SBMData *sbm, HistogramData *hist);
void init_sbm_gui(SBMData *sbm, const char *plotfile_dir, int argc, char **argv);
void render_sbm_histogram(SBMData *sbm);
void update_sbm_info_label(SBMData *sbm, const char *plotfile_dir);
void sbm_metric_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_logx_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_logy_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_pdf_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_switch_timestep(SBMData *sbm, int new_timestep);
int scan_sbm_timesteps(const char *base_dir, const char *prefix);
void sbm_time_nav_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_canvas_expose_callback(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void sbm_settings_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_settings_apply_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sbm_settings_close_callback(Widget w, XtPointer client_data, XtPointer call_data);

/* 1D Profile functions */
int read_profile_file(ProfileFile *pf, const char *path, int file_type);
void free_profile_file(ProfileFile *pf);
void init_profile_gui(ProfileData *pd, int argc, char **argv);
void render_profile_canvas(ProfileData *pd);
void draw_profile_plot(Display *dpy, Window win, GC plot_gc, ProfileData *pd,
                       int file_idx, int col_idx, int time_idx,
                       int width, int height, int log_x);
void draw_profile_contour(Display *dpy, Window win, GC plot_gc, ProfileData *pd,
                          int file_idx, int col_idx, int width, int height);
void update_profile_info_label(ProfileData *pd);
void profile_file_callback(Widget w, XtPointer client_data, XtPointer call_data);
void profile_var_callback(Widget w, XtPointer client_data, XtPointer call_data);
void profile_time_nav_callback(Widget w, XtPointer client_data, XtPointer call_data);
void profile_logx_callback(Widget w, XtPointer client_data, XtPointer call_data);
void profile_contour_callback(Widget w, XtPointer client_data, XtPointer call_data);
void profile_canvas_expose_callback(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void profile_rebuild_var_buttons(ProfileData *pd);

/* Global pointer for callbacks */
PlotfileData *global_pf = NULL;
Widget map_button_widget = NULL;
Widget z_phys_button_widget = NULL;

/* Quiver state and dialog widgets */
typedef struct {
    Widget shell;
    Widget x_comp_text;
    Widget y_comp_text;
    Widget density_label;
    Widget scale_label;
    int x_comp_index;
    int y_comp_index;
    int enabled;
    int density;     /* 1-5, controls arrow spacing */
    double scale;    /* 0.5-3.0, controls arrow length */
    int color;       /* 0=black, 1=white, 2=red, 3=blue */
} QuiverData;

QuiverData quiver_data = {NULL, NULL, NULL, NULL, NULL, -1, -1, 0, 3, 1.0, 0};

/* Variable selection popup data */
typedef struct {
    Widget shell;
    Widget viewport;
    Widget list_widget;
    int selecting_for_x;  /* 1 if selecting for X component, 0 for Y component */
    Widget *var_buttons;  /* Array of buttons for each variable */
    int n_vars;
} VarSelectData;

VarSelectData var_select_data = {NULL, NULL, NULL, 0, NULL, 0};

/* Normalize dimensions for 2D plotfiles */
static void normalize_plotfile_dims(PlotfileData *pf) {
    int i;
    for (i = pf->ndim; i < 3; i++) {
        pf->grid_dims[i] = 1;
        pf->level_lo[i] = 0;
        pf->level_hi[i] = 0;
    }
}

static void normalize_level_dims(PlotfileData *pf, LevelData *ld) {
    int i;
    for (i = pf->ndim; i < 3; i++) {
        ld->grid_dims[i] = 1;
        ld->level_lo[i] = 0;
        ld->level_hi[i] = 0;
    }
}

static void update_var_scrollbar(void) {
    if (!var_viewport || !var_scrollbar) return;

    Widget child = NULL;
    WidgetList children = NULL;
    Cardinal num_children = 0;
    Dimension viewport_height = 0, child_height = 0;
    Position child_y = 0;

    XtVaGetValues(var_viewport,
                  XtNchildren, &children,
                  XtNnumChildren, &num_children,
                  XtNheight, &viewport_height,
                  NULL);
    if (!children || num_children == 0) return;
    child = children[0];
    XtVaGetValues(child, XtNy, &child_y, XtNheight, &child_height, NULL);

    if (child_height <= viewport_height || viewport_height == 0) {
        XawScrollbarSetThumb(var_scrollbar, 0.0f, 1.0f);
        return;
    }

    float shown = (float)viewport_height / (float)child_height;
    float top = (float)(-child_y) / (float)(child_height - viewport_height);
    if (top < 0.0f) top = 0.0f;
    if (top > 1.0f) top = 1.0f;
    XawScrollbarSetThumb(var_scrollbar, top, shown);
}

/* Mouse wheel handler for variable viewport scrolling */
void var_viewport_wheel_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != ButtonPress) return;
    XButtonEvent *bev = (XButtonEvent *)event;
    if (bev->button != Button4 && bev->button != Button5) return;

    Widget child = NULL;
    WidgetList children = NULL;
    Cardinal num_children = 0;
    Dimension viewport_height = 0, child_height = 0;
    Position child_y = 0;

    XtVaGetValues(w,
                  XtNchildren, &children,
                  XtNnumChildren, &num_children,
                  XtNheight, &viewport_height,
                  NULL);
    if (!children || num_children == 0) return;
    child = children[0];
    XtVaGetValues(child, XtNy, &child_y, XtNheight, &child_height, NULL);

    int scroll = (bev->button == Button4) ? -40 : 40;
    int current_y = (int)(-child_y);
    int new_y = current_y + scroll;
    int max_y = (child_height > viewport_height) ? (int)(child_height - viewport_height) : 0;

    if (new_y < 0) new_y = 0;
    if (new_y > max_y) new_y = max_y;

    XawViewportSetCoordinates(w, 0, new_y);
    update_var_scrollbar();
    if (continue_dispatch) *continue_dispatch = False;
}

void var_scrollbar_scroll_proc(Widget w, XtPointer client_data, XtPointer call_data) {
    int delta = (int)(long)call_data;
    Widget child = NULL;
    WidgetList children = NULL;
    Cardinal num_children = 0;
    Dimension viewport_height = 0, child_height = 0;
    Position child_y = 0;

    XtVaGetValues(var_viewport,
                  XtNchildren, &children,
                  XtNnumChildren, &num_children,
                  XtNheight, &viewport_height,
                  NULL);
    if (!children || num_children == 0) return;
    child = children[0];
    XtVaGetValues(child, XtNy, &child_y, XtNheight, &child_height, NULL);

    int current_y = (int)(-child_y);
    int new_y = current_y + delta;
    int max_y = (child_height > viewport_height) ? (int)(child_height - viewport_height) : 0;

    if (new_y < 0) new_y = 0;
    if (new_y > max_y) new_y = max_y;

    XawViewportSetCoordinates(var_viewport, 0, new_y);
    update_var_scrollbar();
}

void var_scrollbar_jump_proc(Widget w, XtPointer client_data, XtPointer call_data) {
    float top = *(float *)call_data;
    Widget child = NULL;
    WidgetList children = NULL;
    Cardinal num_children = 0;
    Dimension viewport_height = 0, child_height = 0;

    XtVaGetValues(var_viewport,
                  XtNchildren, &children,
                  XtNnumChildren, &num_children,
                  XtNheight, &viewport_height,
                  NULL);
    if (!children || num_children == 0) return;
    child = children[0];
    XtVaGetValues(child, XtNheight, &child_height, NULL);

    int max_y = (child_height > viewport_height) ? (int)(child_height - viewport_height) : 0;
    int new_y = (int)(top * max_y + 0.5f);
    if (new_y < 0) new_y = 0;
    if (new_y > max_y) new_y = max_y;

    XawViewportSetCoordinates(var_viewport, 0, new_y);
    update_var_scrollbar();
}

void var_viewport_configure_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type == ConfigureNotify) {
        update_var_scrollbar();
    }
}


/* Detect number of levels by scanning for Level_X directories */
int detect_levels(PlotfileData *pf) {
    char path[MAX_PATH];
    int level = 0;

    /* Count how many Level_X directories exist */
    while (level < 100) {
        snprintf(path, MAX_PATH, "%s/Level_%d", pf->plotfile_dir, level);
        DIR *dir = opendir(path);
        if (!dir) break;
        closedir(dir);
        level++;
    }

    return level > 0 ? level : 1;
}

/* Detect number of levels for a given plotfile path */
int detect_levels_for_path(const char *plotfile_dir) {
    char path[MAX_PATH];
    int level = 0;

    /* Count how many Level_X directories exist */
    while (level < 100) {
        snprintf(path, MAX_PATH, "%s/Level_%d", plotfile_dir, level);
        DIR *dir = opendir(path);
        if (!dir) break;
        closedir(dir);
        level++;
    }

    return level > 0 ? level : 1;
}

/* Show warning popup when level is not available at current timestep */
/* Callback to close warning popup */
void warning_ok_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Widget shell = (Widget)client_data;
    if (shell) {
        XtDestroyWidget(shell);
    }
}

void show_level_warning(int level) {
    Arg args[10];
    int n;
    char msg[128];

    snprintf(msg, sizeof(msg), "Level %d not available at this timestep", level);

    /* Create popup shell */
    n = 0;
    XtSetArg(args[n], XtNtitle, "Warning"); n++;
    Widget warning_shell = XtCreatePopupShell("levelWarning", transientShellWidgetClass, toplevel, args, n);

    /* Create form container */
    Widget form = XtVaCreateManagedWidget("form", formWidgetClass, warning_shell, NULL);

    /* Create warning label */
    n = 0;
    XtSetArg(args[n], XtNlabel, msg); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget warning_label = XtCreateManagedWidget("warningLabel", labelWidgetClass, form, args, n);

    /* Create OK button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "OK"); n++;
    XtSetArg(args[n], XtNfromVert, warning_label); n++;
    Widget ok_button = XtCreateManagedWidget("okButton", commandWidgetClass, form, args, n);

    /* OK button callback - destroy the popup */
    XtAddCallback(ok_button, XtNcallback, warning_ok_callback, (XtPointer)warning_shell);

    XtPopup(warning_shell, XtGrabExclusive);
}

/* Comparison function for sorting timesteps */
int compare_timesteps(const void *a, const void *b) {
    int idx_a = *(const int *)a;
    int idx_b = *(const int *)b;
    return timestep_numbers[idx_a] - timestep_numbers[idx_b];
}

/* Scan directory for plotfiles and sort them by number */
int scan_timesteps(const char *base_dir, const char *prefix) {
    DIR *dir;
    struct dirent *entry;
    char check_path[MAX_PATH];
    int indices[MAX_TIMESTEPS];
    int temp_levels[MAX_TIMESTEPS];
    int prefix_len = strlen(prefix);

    dir = opendir(base_dir);
    if (!dir) {
        fprintf(stderr, "Error: Cannot open directory %s\n", base_dir);
        return -1;
    }

    n_timesteps = 0;
    max_levels_all_timesteps = 1;

    while ((entry = readdir(dir)) != NULL && n_timesteps < MAX_TIMESTEPS) {
        /* Check if entry starts with the specified prefix */
        if (strncmp(entry->d_name, prefix, prefix_len) == 0) {
            /* Ensure ALL characters after prefix are digits (to avoid plt matching plt2d) */
            const char *suffix = entry->d_name + prefix_len;
            int all_digits = 1;
            if (*suffix == '\0') all_digits = 0;  /* Must have at least one digit */
            for (const char *p = suffix; *p != '\0'; p++) {
                if (!isdigit(*p)) {
                    all_digits = 0;
                    break;
                }
            }
            if (!all_digits) continue;

            /* Check if it's a valid plotfile directory (has Header file) */
            snprintf(check_path, MAX_PATH, "%s/%s/Header", base_dir, entry->d_name);
            FILE *fp = fopen(check_path, "r");
            if (fp) {
                fclose(fp);

                /* Extract number from plotfile name (after prefix) */
                int num = atoi(entry->d_name + prefix_len);

                /* Allocate and store path */
                timestep_paths[n_timesteps] = (char *)malloc(MAX_PATH);
                snprintf(timestep_paths[n_timesteps], MAX_PATH, "%s/%s", base_dir, entry->d_name);
                timestep_numbers[n_timesteps] = num;

                /* Detect levels for this timestep */
                int levels = detect_levels_for_path(timestep_paths[n_timesteps]);
                temp_levels[n_timesteps] = levels;
                if (levels > max_levels_all_timesteps) {
                    max_levels_all_timesteps = levels;
                }

                indices[n_timesteps] = n_timesteps;
                n_timesteps++;
            }
        }
    }

    closedir(dir);

    if (n_timesteps == 0) {
        return -1;
    }

    /* Sort indices by timestep number */
    qsort(indices, n_timesteps, sizeof(int), compare_timesteps);

    /* Reorder arrays based on sorted indices */
    char *temp_paths[MAX_TIMESTEPS];
    int temp_numbers[MAX_TIMESTEPS];
    int sorted_levels[MAX_TIMESTEPS];

    for (int i = 0; i < n_timesteps; i++) {
        temp_paths[i] = timestep_paths[indices[i]];
        temp_numbers[i] = timestep_numbers[indices[i]];
        sorted_levels[i] = temp_levels[indices[i]];
    }

    for (int i = 0; i < n_timesteps; i++) {
        timestep_paths[i] = temp_paths[i];
        timestep_numbers[i] = temp_numbers[i];
        timestep_levels[i] = sorted_levels[i];
    }

    printf("Found %d timesteps, max levels across all: %d\n", n_timesteps, max_levels_all_timesteps);
    return n_timesteps;
}

/* Scan directory for plotfiles with SDM data and sort them by number */
int scan_sdm_timesteps(const char *base_dir, const char *prefix) {
    DIR *dir;
    struct dirent *entry;
    char check_path[MAX_PATH];
    int indices[MAX_TIMESTEPS];
    int prefix_len = strlen(prefix);

    dir = opendir(base_dir);
    if (!dir) {
        fprintf(stderr, "Error: Cannot open directory %s\n", base_dir);
        return -1;
    }

    n_timesteps = 0;

    while ((entry = readdir(dir)) != NULL && n_timesteps < MAX_TIMESTEPS) {
        if (strncmp(entry->d_name, prefix, prefix_len) == 0) {
            const char *suffix = entry->d_name + prefix_len;
            int all_digits = 1;
            if (*suffix == '\0') all_digits = 0;
            for (const char *p = suffix; *p != '\0'; p++) {
                if (!isdigit(*p)) {
                    all_digits = 0;
                    break;
                }
            }
            if (!all_digits) continue;

            /* Check for SDM subdirectory Header */
            snprintf(check_path, MAX_PATH, "%s/%s/%s/Header",
                     base_dir, entry->d_name, SDM_SUBDIR);
            FILE *fp = fopen(check_path, "r");
            if (fp) {
                fclose(fp);

                int num = atoi(entry->d_name + prefix_len);

                timestep_paths[n_timesteps] = (char *)malloc(MAX_PATH);
                snprintf(timestep_paths[n_timesteps], MAX_PATH, "%s/%s", base_dir, entry->d_name);
                timestep_numbers[n_timesteps] = num;
                indices[n_timesteps] = n_timesteps;
                n_timesteps++;
            }
        }
    }

    closedir(dir);

    if (n_timesteps == 0) {
        return -1;
    }

    /* Sort indices by timestep number */
    qsort(indices, n_timesteps, sizeof(int), compare_timesteps);

    char *temp_paths[MAX_TIMESTEPS];
    int temp_numbers[MAX_TIMESTEPS];

    for (int i = 0; i < n_timesteps; i++) {
        temp_paths[i] = timestep_paths[indices[i]];
        temp_numbers[i] = timestep_numbers[indices[i]];
    }

    for (int i = 0; i < n_timesteps; i++) {
        timestep_paths[i] = temp_paths[i];
        timestep_numbers[i] = temp_numbers[i];
    }

    printf("Found %d SDM timesteps\n", n_timesteps);
    return n_timesteps;
}

/* Scan directory for SBM plotfiles (those with bin_info.txt) */
int scan_sbm_timesteps(const char *base_dir, const char *prefix) {
    DIR *dir;
    struct dirent *entry;
    char check_path[MAX_PATH];
    int indices[MAX_TIMESTEPS];
    int prefix_len = strlen(prefix);

    dir = opendir(base_dir);
    if (!dir) {
        fprintf(stderr, "Error: Cannot open directory %s\n", base_dir);
        return -1;
    }

    n_timesteps = 0;

    while ((entry = readdir(dir)) != NULL && n_timesteps < MAX_TIMESTEPS) {
        if (strncmp(entry->d_name, prefix, prefix_len) == 0) {
            const char *suffix = entry->d_name + prefix_len;
            int all_digits = 1;
            if (*suffix == '\0') all_digits = 0;
            for (const char *p = suffix; *p != '\0'; p++) {
                if (!isdigit(*p)) {
                    all_digits = 0;
                    break;
                }
            }
            if (!all_digits) continue;

            /* Check for bin_info.txt in this plotfile */
            snprintf(check_path, MAX_PATH, "%s/%s/%s",
                     base_dir, entry->d_name, SBM_BIN_INFO_FILE);
            FILE *fp = fopen(check_path, "r");
            if (fp) {
                fclose(fp);

                int num = atoi(entry->d_name + prefix_len);

                timestep_paths[n_timesteps] = (char *)malloc(MAX_PATH);
                snprintf(timestep_paths[n_timesteps], MAX_PATH, "%s/%s", base_dir, entry->d_name);
                timestep_numbers[n_timesteps] = num;
                indices[n_timesteps] = n_timesteps;
                n_timesteps++;
            }
        }
    }

    closedir(dir);

    if (n_timesteps == 0) {
        return -1;
    }

    /* Sort indices by timestep number */
    qsort(indices, n_timesteps, sizeof(int), compare_timesteps);

    char *temp_paths[MAX_TIMESTEPS];
    int temp_numbers[MAX_TIMESTEPS];

    for (int i = 0; i < n_timesteps; i++) {
        temp_paths[i] = timestep_paths[indices[i]];
        temp_numbers[i] = timestep_numbers[indices[i]];
    }

    for (int i = 0; i < n_timesteps; i++) {
        timestep_paths[i] = temp_paths[i];
        timestep_numbers[i] = temp_numbers[i];
    }

    printf("Found %d SBM timesteps\n", n_timesteps);
    return n_timesteps;
}

/* Switch to a different timestep */
void switch_timestep(PlotfileData *pf, int new_timestep) {
    if (new_timestep < 0 || new_timestep >= n_timesteps) return;

    zoom_reset();
    current_timestep = new_timestep;

    /* Update plotfile directory */
    strncpy(pf->plotfile_dir, timestep_paths[current_timestep], MAX_PATH - 1);

    /* Always free old overlay data before reading new timestep */
    free_all_levels(pf);
    free_z_phys_cache(pf);

    /* Save overlay_mode and map_mode before read_header (which resets them) */
    int saved_overlay_mode = pf->overlay_mode;
    int saved_map_mode = pf->map_mode;

    /* Re-read header for new timestep */
    read_header(pf);

    /* Restore overlay_mode and map_mode */
    pf->overlay_mode = saved_overlay_mode;
    pf->map_mode = saved_map_mode;
    if (map_button_widget)
        XtVaSetValues(map_button_widget, XtNlabel, pf->map_mode ? "Map: ON" : "Map: OFF", NULL);

    /* Clamp current_level if new timestep has fewer levels */
    if (pf->current_level >= pf->n_levels) {
        pf->current_level = pf->n_levels - 1;
        if (pf->current_level < 0) pf->current_level = 0;
    }

    /* Reset boxes and re-read cell data */
    pf->n_boxes = 0;
    read_cell_h(pf);

    /* Clamp slice_idx if new data has fewer layers */
    int max_idx = pf->grid_dims[pf->slice_axis] - 1;
    if (pf->slice_idx > max_idx) {
        pf->slice_idx = max_idx;
    }

    /* Re-read variable data */
    read_variable_data(pf, pf->current_var);

    /* If overlay mode is on, reload all levels for new timestep */
    /* Don't change overlay_mode or button label - just reload data if needed */
    printf("switch_timestep: overlay_mode=%d, n_levels=%d\n", pf->overlay_mode, pf->n_levels);
    if (pf->overlay_mode && pf->n_levels > 1) {
        printf("switch_timestep: Reloading overlay levels...\n");
        load_all_levels(pf, pf->current_var);
    }

    /* Update UI */
    update_time_label();
    update_layer_label(pf);
    update_z_phys_button(pf);
    update_info_label(pf);
    render_slice(pf);
}

/* Update time step label */
void update_time_label(void) {
    if (time_label && n_timesteps > 0) {
        char text[32];
        snprintf(text, sizeof(text), "%d/%d", current_timestep + 1, n_timesteps);
        Arg args[1];
        XtSetArg(args[0], XtNlabel, text);
        XtSetValues(time_label, args, 1);
    }
}

/* Time navigation button callback */
void time_nav_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int dir = (int)(long)client_data;  /* 0 = prev (<), 1 = next (>) */

    if (global_pf && n_timesteps > 1) {
        int new_timestep = current_timestep;

        if (dir == 1) {
            /* Next: go to next timestep, wrap to 0 if at end */
            new_timestep++;
            if (new_timestep >= n_timesteps) {
                new_timestep = 0;
            }
        } else {
            /* Prev: go to previous timestep, wrap to end if at 0 */
            new_timestep--;
            if (new_timestep < 0) {
                new_timestep = n_timesteps - 1;
            }
        }

        switch_timestep(global_pf, new_timestep);
    }
}

/* Jump to specific timestep - button-based callback */
void time_jump_to_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    long time_type = (long)client_data;

    if (global_pf && n_timesteps > 1) {
        int new_timestep = current_timestep;

        switch (time_type) {
            case 0: new_timestep = 0; break;                    /* First */
            case 1: new_timestep = n_timesteps - 1; break;      /* Last */
            case 2: new_timestep = n_timesteps / 2; break;      /* Middle */
            case 3: new_timestep = n_timesteps / 4; break;      /* 1/4 */
            case 4: new_timestep = 3 * n_timesteps / 4; break;  /* 3/4 */
        }

        if (new_timestep >= 0 && new_timestep < n_timesteps) {
            switch_timestep(global_pf, new_timestep);
        }
    }

    /* Close the dialog */
    Widget shell = XtParent(XtParent(w));
    XtPopdown(shell);
    XtDestroyWidget(shell);
}

/* Structure to pass both text widget and shell to callback */
typedef struct {
    Widget text_widget;
    Widget dialog_shell;
} TimeJumpDialogData;

/* Jump to typed timestep number */
void time_jump_to_typed_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    TimeJumpDialogData *data = (TimeJumpDialogData *)client_data;

    if (global_pf && data && n_timesteps > 1) {
        String value;
        Arg args[1];
        XtSetArg(args[0], XtNstring, &value);
        XtGetValues(data->text_widget, args, 1);

        if (value && strlen(value) > 0) {
            int timestep = atoi(value);

            /* Convert from 1-indexed to 0-indexed and clamp */
            timestep = timestep - 1;
            if (timestep < 0) timestep = 0;
            if (timestep >= n_timesteps) timestep = n_timesteps - 1;

            switch_timestep(global_pf, timestep);
        }

        /* Close the dialog */
        XtPopdown(data->dialog_shell);
        XtDestroyWidget(data->dialog_shell);
        free(data);
    }
}

/* Close time jump dialog */
void time_jump_dialog_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Widget shell = (Widget)client_data;
    XtPopdown(shell);
    XtDestroyWidget(shell);
}

/* Time Jump button callback */
void time_jump_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && n_timesteps > 1) {
        Arg args[10];
        int n;
        Widget dialog_shell, form, label, button, text_widget, text_label;
        char msg[128];

        snprintf(msg, sizeof(msg), "Jump to timestep (1-%d)", n_timesteps);

        n = 0;
        XtSetArg(args[n], XtNtitle, "Jump to Timestep"); n++;
        dialog_shell = XtCreatePopupShell("timeJumpDialog", transientShellWidgetClass, toplevel, args, n);

        n = 0;
        form = XtCreateManagedWidget("form", formWidgetClass, dialog_shell, args, n);

        /* Title label */
        n = 0;
        XtSetArg(args[n], XtNlabel, msg); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("label", labelWidgetClass, form, args, n);

        /* Text input section */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNlabel, "Type timestep:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        text_label = XtCreateManagedWidget("textLabel", labelWidgetClass, form, args, n);

        n = 0;
        XtSetArg(args[n], XtNfromVert, text_label); n++;
        XtSetArg(args[n], XtNwidth, 100); n++;
        XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
        XtSetArg(args[n], XtNstring, ""); n++;
        text_widget = XtCreateManagedWidget("textInput", asciiTextWidgetClass, form, args, n);

        /* Create data structure to pass to callback */
        TimeJumpDialogData *jump_data = malloc(sizeof(TimeJumpDialogData));
        jump_data->text_widget = text_widget;
        jump_data->dialog_shell = dialog_shell;

        n = 0;
        XtSetArg(args[n], XtNfromVert, text_label); n++;
        XtSetArg(args[n], XtNfromHoriz, text_widget); n++;
        XtSetArg(args[n], XtNlabel, "Go"); n++;
        button = XtCreateManagedWidget("goButton", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_to_typed_callback, (XtPointer)jump_data);

        /* Or quick jump label */
        n = 0;
        XtSetArg(args[n], XtNfromVert, text_widget); n++;
        XtSetArg(args[n], XtNlabel, "Or quick jump:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("orLabel", labelWidgetClass, form, args, n);

        /* Quick jump buttons */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNlabel, "First (1)"); n++;
        button = XtCreateManagedWidget("first", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_to_callback, (XtPointer)0);

        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "1/4"); n++;
        button = XtCreateManagedWidget("quarter", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_to_callback, (XtPointer)3);

        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "Middle"); n++;
        button = XtCreateManagedWidget("middle", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_to_callback, (XtPointer)2);

        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "3/4"); n++;
        button = XtCreateManagedWidget("threequarter", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_to_callback, (XtPointer)4);

        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        snprintf(msg, sizeof(msg), "Last (%d)", n_timesteps);
        XtSetArg(args[n], XtNlabel, msg); n++;
        button = XtCreateManagedWidget("last", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_to_callback, (XtPointer)1);

        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "Close"); n++;
        button = XtCreateManagedWidget("close", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, time_jump_dialog_close_callback, (XtPointer)dialog_shell);

        XtRealizeWidget(dialog_shell);
        XtPopup(dialog_shell, XtGrabExclusive);

        /* Set keyboard focus to text input */
        XtSetKeyboardFocus(dialog_shell, text_widget);
        XSync(display, False);
        Time time = CurrentTime;
        XtCallAcceptFocus(text_widget, &time);
    }
}

/* Time Series button callback */
void time_series_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data && n_timesteps > 1) {
        show_time_series(global_pf);
    }
}

/* Read Header file */
int read_header(PlotfileData *pf) {
    char path[MAX_PATH];
    char line[MAX_LINE];
    FILE *fp;
    int i, idx = 0;
    
    snprintf(path, MAX_PATH, "%s/Header", pf->plotfile_dir);
    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", path);
        return -1;
    }
    
    /* Line 0: version */
    fgets(line, MAX_LINE, fp);
    
    /* Line 1: number of variables */
    fgets(line, MAX_LINE, fp);
    pf->n_vars = atoi(line);
    
    /* Variable names */
    for (i = 0; i < pf->n_vars; i++) {
        fgets(line, MAX_LINE, fp);
        line[strcspn(line, "\n")] = 0;  /* Remove newline */
        strncpy(pf->variables[i], line, 63);
    }
    
    /* Dimensionality */
    fgets(line, MAX_LINE, fp);
    pf->ndim = atoi(line);
    
    /* Time */
    fgets(line, MAX_LINE, fp);
    pf->time = atof(line);
    
    /* Number of levels - read from Header but verify by scanning directories */
    fgets(line, MAX_LINE, fp);
    int header_levels = atoi(line);

    /* Read prob_lo (domain lower bounds) */
    fgets(line, MAX_LINE, fp);
    pf->prob_lo[0] = pf->prob_lo[1] = pf->prob_lo[2] = 0.0;
    if (pf->ndim == 3) sscanf(line, "%lf %lf %lf", &pf->prob_lo[0], &pf->prob_lo[1], &pf->prob_lo[2]);
    else if (pf->ndim == 2) sscanf(line, "%lf %lf", &pf->prob_lo[0], &pf->prob_lo[1]);

    /* Read prob_hi (domain upper bounds) */
    fgets(line, MAX_LINE, fp);
    pf->prob_hi[0] = pf->prob_hi[1] = pf->prob_hi[2] = 0.0;
    if (pf->ndim == 3) sscanf(line, "%lf %lf %lf", &pf->prob_hi[0], &pf->prob_hi[1], &pf->prob_hi[2]);
    else if (pf->ndim == 2) sscanf(line, "%lf %lf", &pf->prob_hi[0], &pf->prob_hi[1]);

    /* Parse refinement ratios - format: "r1 r2 ..." or "(r1,r2,r3) ..." */
    fgets(line, MAX_LINE, fp);
    pf->ref_ratio[0] = 1;  /* Level 0 has no refinement */
    int ref = atoi(line);  /* First number is ref ratio for level 0->1 */
    for (i = 1; i < MAX_LEVELS; i++) {
        pf->ref_ratio[i] = (ref > 0) ? ref : 2;  /* Default to 2 if not parsed */
    }
    /* Initialize overlay mode to off */
    pf->overlay_mode = 0;
    pf->map_mode = 0;
    
    /* Domain box */
    fgets(line, MAX_LINE, fp);
    /* Parse ((lo_x,lo_y,lo_z) (hi_x,hi_y,hi_z) ...) */
    char *p = line;
    int lo[3], hi[3];
    while (*p && (*p == '(' || *p == ' ')) p++;
    for (i = 0; i < pf->ndim; i++) {
        while (*p && !isdigit(*p) && *p != '-') p++;
        lo[i] = atoi(p);
        while (*p && (isdigit(*p) || *p == '-')) p++;
    }
    for (i = 0; i < pf->ndim; i++) {
        while (*p && !isdigit(*p) && *p != '-') p++;
        hi[i] = atoi(p);
        while (*p && (isdigit(*p) || *p == '-')) p++;
    }
    
    for (i = 0; i < pf->ndim; i++) {
        pf->grid_dims[i] = hi[i] - lo[i] + 1;
    }
    /* Ensure non-zero dimensions for 2D plotfiles */
    normalize_plotfile_dims(pf);
    
    fclose(fp);
    
    /* Detect actual levels by scanning directories */
    pf->n_levels = detect_levels(pf);
    
    printf("Loaded: %s\n", pf->plotfile_dir);
    printf("Variables: %d (", pf->n_vars);
    for (i = 0; i < pf->n_vars && i < 5; i++) {
        printf("%s%s", pf->variables[i], i < pf->n_vars-1 ? ", " : "");
    }
    if (pf->n_vars > 5) printf("...");
    printf(")\n");
    printf("Grid: %d x %d x %d\n", pf->grid_dims[0], pf->grid_dims[1], pf->grid_dims[2]);
    printf("Domain: [%.3g, %.3g] x [%.3g, %.3g] x [%.3g, %.3g]\n",
           pf->prob_lo[0], pf->prob_hi[0], pf->prob_lo[1], pf->prob_hi[1],
           pf->prob_lo[2], pf->prob_hi[2]);
    printf("Time: %.3f\n", pf->time);
    printf("Levels: %d (Header says %d)\n", pf->n_levels, header_levels);
    
    return 0;
}

/* Read Cell_H to get box layout and FabOnDisk mapping */
int read_cell_h(PlotfileData *pf) {
    char path[MAX_PATH];
    char line[MAX_LINE];
    FILE *fp;
    int i;
    
    snprintf(path, MAX_PATH, "%s/Level_%d/Cell_H", pf->plotfile_dir, pf->current_level);
    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", path);
        return -1;
    }
    
    /* Reset grid dimensions for this level */
    int level_lo[3] = {0, 0, 0};
    int level_hi[3] = {0, 0, 0};
    int found_domain = 0;
    pf->n_boxes = 0;
    
    /* Skip first few lines until we find box definitions */
    int box_count = 0;
    while (fgets(line, MAX_LINE, fp)) {
        if (strncmp(line, "((", 2) == 0) {
            /* Parse box: ((lo_x,lo_y,lo_z) (hi_x,hi_y,hi_z) ...) */
            char *p = line + 2;
            int lo[3], hi[3];
            for (i = 0; i < pf->ndim; i++) {
                while (*p && !isdigit(*p) && *p != '-') p++;
                lo[i] = atoi(p);
                pf->boxes[box_count].lo[i] = lo[i];
                while (*p && (isdigit(*p) || *p == '-')) p++;
            }
            for (i = 0; i < pf->ndim; i++) {
                while (*p && !isdigit(*p) && *p != '-') p++;
                hi[i] = atoi(p);
                pf->boxes[box_count].hi[i] = hi[i];
                while (*p && (isdigit(*p) || *p == '-')) p++;
            }
            
            /* Track overall domain bounds */
            if (!found_domain) {
                for (i = 0; i < pf->ndim; i++) {
                    level_lo[i] = lo[i];
                    level_hi[i] = hi[i];
                }
                found_domain = 1;
            } else {
                for (i = 0; i < pf->ndim; i++) {
                    if (lo[i] < level_lo[i]) level_lo[i] = lo[i];
                    if (hi[i] > level_hi[i]) level_hi[i] = hi[i];
                }
            }
            
            box_count++;
        } else if (strncmp(line, "FabOnDisk:", 10) == 0) {
            /* Parse FabOnDisk: Cell_D_XXXXX offset */
            char *p = strchr(line, ':');
            if (p) {
                p++;
                while (*p == ' ') p++;
                char *end = strchr(p, ' ');
                if (end) {
                    *end = '\0';
                    end++;
                }
                char *line_end = strchr(p, '\n');
                if (line_end) *line_end = '\0';
                strncpy(pf->boxes[pf->n_boxes].filename, p, 63);
                pf->boxes[pf->n_boxes].filename[63] = '\0';
                pf->boxes[pf->n_boxes].offset = 0;
                if (end) {
                    while (*end == ' ') end++;
                    if (*end) {
                        pf->boxes[pf->n_boxes].offset = strtoll(end, NULL, 10);
                    }
                }
                pf->n_boxes++;
            }
        }
    }
    
    fclose(fp);

    /* Update grid dimensions and level bounds */
    for (i = 0; i < pf->ndim; i++) {
        pf->level_lo[i] = level_lo[i];
        pf->level_hi[i] = level_hi[i];
        pf->grid_dims[i] = level_hi[i] - level_lo[i] + 1;
    }
    /* Ensure non-zero dimensions for 2D plotfiles */
    normalize_plotfile_dims(pf);

    printf("Level %d: Found %d boxes, Grid: %d x %d x %d (lo: %d,%d,%d)\n",
           pf->current_level, pf->n_boxes,
           pf->grid_dims[0], pf->grid_dims[1], pf->grid_dims[2],
           pf->level_lo[0], pf->level_lo[1], pf->level_lo[2]);
    return 0;
}

/* Read variable data from all boxes */
int read_variable_data(PlotfileData *pf, int var_idx) {
    char path[MAX_PATH];
    FILE *fp;
    int box_idx, i, j, k;
    size_t total_size = pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2];
    
    /* Allocate data array (Z, Y, X ordering) */
    if (pf->data) free(pf->data);
    pf->data = (double *)calloc(total_size, sizeof(double));
    
    /* Read each box */
    for (box_idx = 0; box_idx < pf->n_boxes; box_idx++) {
        Box *box = &pf->boxes[box_idx];
        int box_dims[3];
        for (i = 0; i < 3; i++) {
            box_dims[i] = box->hi[i] - box->lo[i] + 1;
        }
        size_t box_size = box_dims[0] * box_dims[1] * box_dims[2];
        
        snprintf(path, MAX_PATH, "%s/Level_%d/%s", pf->plotfile_dir, pf->current_level, box->filename);
        fp = fopen(path, "rb");
        if (!fp) continue;
        
        /* Seek to FAB header offset if provided */
        if (box->offset > 0) {
            fseeko(fp, (off_t)box->offset, SEEK_SET);
        }

        /* Skip FAB header (read until newline) */
        int c;
        while ((c = fgetc(fp)) != EOF && c != '\n');

        /* Skip to variable data */
        fseeko(fp, (off_t)(var_idx * box_size * sizeof(double)), SEEK_CUR);
        
        /* Read box data */
        double *box_data = (double *)malloc(box_size * sizeof(double));
        fread(box_data, sizeof(double), box_size, fp);
        fclose(fp);
        
        /* Insert into global array (Fortran order -> C order) */
        /* Fortran order: X varies fastest */
        /* Use relative indices by subtracting level_lo to handle non-zero level origins */
        size_t idx = 0;
        for (k = 0; k < box_dims[2]; k++) {
            for (j = 0; j < box_dims[1]; j++) {
                for (i = 0; i < box_dims[0]; i++) {
                    int gx = box->lo[0] + i - pf->level_lo[0];
                    int gy = box->lo[1] + j - pf->level_lo[1];
                    int gz = box->lo[2] + k - pf->level_lo[2];
                    /* Global array: data[z][y][x] */
                    size_t gidx = gz * pf->grid_dims[1] * pf->grid_dims[0] +
                                  gy * pf->grid_dims[0] + gx;
                    pf->data[gidx] = box_data[idx++];
                }
            }
        }
        
        free(box_data);
    }
    
    printf("Loaded variable: %s\n", pf->variables[var_idx]);
    return 0;
}

/* ========== Multi-Level Overlay Functions ========== */

/* Read Cell_H for a specific level into LevelData */
int read_cell_h_level(PlotfileData *pf, int level) {
    char path[MAX_PATH];
    char line[MAX_LINE];
    FILE *fp;
    int i;
    LevelData *ld = &pf->levels[level];

    snprintf(path, MAX_PATH, "%s/Level_%d/Cell_H", pf->plotfile_dir, level);
    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", path);
        return -1;
    }

    /* Reset level data */
    int level_lo[3] = {0, 0, 0};
    int level_hi[3] = {0, 0, 0};
    int found_domain = 0;
    ld->n_boxes = 0;

    /* Parse box definitions and FabOnDisk entries */
    int box_count = 0;
    while (fgets(line, MAX_LINE, fp)) {
        if (strncmp(line, "((", 2) == 0) {
            /* Parse box: ((lo_x,lo_y,lo_z) (hi_x,hi_y,hi_z) ...) */
            char *p = line + 2;
            int lo[3], hi[3];
            for (i = 0; i < pf->ndim; i++) {
                while (*p && !isdigit(*p) && *p != '-') p++;
                lo[i] = atoi(p);
                ld->boxes[box_count].lo[i] = lo[i];
                while (*p && (isdigit(*p) || *p == '-')) p++;
            }
            for (i = 0; i < pf->ndim; i++) {
                while (*p && !isdigit(*p) && *p != '-') p++;
                hi[i] = atoi(p);
                ld->boxes[box_count].hi[i] = hi[i];
                while (*p && (isdigit(*p) || *p == '-')) p++;
            }

            /* Track overall domain bounds */
            if (!found_domain) {
                for (i = 0; i < pf->ndim; i++) {
                    level_lo[i] = lo[i];
                    level_hi[i] = hi[i];
                }
                found_domain = 1;
            } else {
                for (i = 0; i < pf->ndim; i++) {
                    if (lo[i] < level_lo[i]) level_lo[i] = lo[i];
                    if (hi[i] > level_hi[i]) level_hi[i] = hi[i];
                }
            }
            box_count++;
        } else if (strncmp(line, "FabOnDisk:", 10) == 0) {
            /* Parse FabOnDisk: Cell_D_XXXXX offset */
            char *p = strchr(line, ':');
            if (p) {
                p++;
                while (*p == ' ') p++;
                char *end = strchr(p, ' ');
                if (end) {
                    *end = '\0';
                    end++;
                }
                char *line_end = strchr(p, '\n');
                if (line_end) *line_end = '\0';
                strncpy(ld->boxes[ld->n_boxes].filename, p, 63);
                ld->boxes[ld->n_boxes].filename[63] = '\0';
                ld->boxes[ld->n_boxes].offset = 0;
                if (end) {
                    while (*end == ' ') end++;
                    if (*end) {
                        ld->boxes[ld->n_boxes].offset = strtoll(end, NULL, 10);
                    }
                }
                ld->n_boxes++;
            }
        }
    }

    fclose(fp);

    /* Store level bounds and grid dimensions */
    for (i = 0; i < pf->ndim; i++) {
        ld->level_lo[i] = level_lo[i];
        ld->level_hi[i] = level_hi[i];
        ld->grid_dims[i] = level_hi[i] - level_lo[i] + 1;
    }
    normalize_level_dims(pf, ld);

    printf("Level %d overlay: Found %d boxes, Grid: %d x %d x %d (lo: %d,%d,%d)\n",
           level, ld->n_boxes,
           ld->grid_dims[0], ld->grid_dims[1], ld->grid_dims[2],
           ld->level_lo[0], ld->level_lo[1], ld->level_lo[2]);
    return 0;
}

/* Read variable data for a specific level into LevelData */
int read_variable_data_level(PlotfileData *pf, int var_idx, int level) {
    char path[MAX_PATH];
    FILE *fp;
    int box_idx, i, j, k;
    LevelData *ld = &pf->levels[level];

    size_t total_size = (size_t)ld->grid_dims[0] * ld->grid_dims[1] * ld->grid_dims[2];

    /* Allocate data array */
    if (ld->data) free(ld->data);
    ld->data = (double *)calloc(total_size, sizeof(double));
    if (!ld->data) {
        fprintf(stderr, "Error: Cannot allocate memory for level %d\n", level);
        return -1;
    }

    /* Read each box */
    for (box_idx = 0; box_idx < ld->n_boxes; box_idx++) {
        Box *box = &ld->boxes[box_idx];
        int box_dims[3];
        for (i = 0; i < 3; i++) {
            box_dims[i] = box->hi[i] - box->lo[i] + 1;
        }
        size_t box_size = box_dims[0] * box_dims[1] * box_dims[2];

        snprintf(path, MAX_PATH, "%s/Level_%d/%s", pf->plotfile_dir, level, box->filename);
        fp = fopen(path, "rb");
        if (!fp) continue;

        /* Seek to FAB header offset if provided */
        if (box->offset > 0) {
            fseeko(fp, (off_t)box->offset, SEEK_SET);
        }

        /* Skip FAB header */
        int c;
        while ((c = fgetc(fp)) != EOF && c != '\n');

        /* Skip to variable data */
        fseeko(fp, (off_t)(var_idx * box_size * sizeof(double)), SEEK_CUR);

        /* Read box data */
        double *box_data = (double *)malloc(box_size * sizeof(double));
        fread(box_data, sizeof(double), box_size, fp);
        fclose(fp);

        /* Insert into level array using relative indices */
        size_t idx = 0;
        for (k = 0; k < box_dims[2]; k++) {
            for (j = 0; j < box_dims[1]; j++) {
                for (i = 0; i < box_dims[0]; i++) {
                    int gx = box->lo[0] + i - ld->level_lo[0];
                    int gy = box->lo[1] + j - ld->level_lo[1];
                    int gz = box->lo[2] + k - ld->level_lo[2];
                    size_t gidx = gz * ld->grid_dims[1] * ld->grid_dims[0] +
                                  gy * ld->grid_dims[0] + gx;
                    ld->data[gidx] = box_data[idx++];
                }
            }
        }

        free(box_data);
    }

    ld->loaded = 1;
    printf("Loaded level %d: %s\n", level, pf->variables[var_idx]);
    return 0;
}

/* Load all levels for overlay rendering */
int load_all_levels(PlotfileData *pf, int var_idx) {
    int level;
    int loaded_count = 0;

    printf("load_all_levels: Loading %d levels for var %d\n", pf->n_levels, var_idx);

    for (level = 0; level < pf->n_levels && level < MAX_LEVELS; level++) {
        /* Always read Cell_H for this level to ensure fresh data */
        if (read_cell_h_level(pf, level) < 0) {
            fprintf(stderr, "Warning: Cannot read Cell_H for level %d\n", level);
            continue;
        }

        /* Read variable data for this level */
        if (read_variable_data_level(pf, var_idx, level) < 0) {
            fprintf(stderr, "Warning: Cannot load variable for level %d\n", level);
            continue;
        }

        loaded_count++;
    }

    printf("Loaded %d of %d levels for overlay\n", loaded_count, pf->n_levels);
    return 0;
}

/* Free all level data */
void free_all_levels(PlotfileData *pf) {
    int level, i;
    for (level = 0; level < MAX_LEVELS; level++) {
        if (pf->levels[level].data) {
            free(pf->levels[level].data);
            pf->levels[level].data = NULL;
        }
        if (pf->levels[level].z_phys_data) {
            free(pf->levels[level].z_phys_data);
            pf->levels[level].z_phys_data = NULL;
        }
        pf->levels[level].loaded = 0;
        pf->levels[level].z_phys_loaded = 0;
        pf->levels[level].n_boxes = 0;
        /* Clear all fields to prevent stale data issues */
        for (i = 0; i < 3; i++) {
            pf->levels[level].grid_dims[i] = 0;
            pf->levels[level].level_lo[i] = 0;
            pf->levels[level].level_hi[i] = 0;
        }
    }
}

/* ========== SDM (Super Droplet Moisture) Functions ========== */

/* Read particle Header from super_droplets_moisture subdirectory */
int read_sdm_header(ParticleData *pd, const char *plotfile_dir) {
    char path[MAX_PATH];
    char line[MAX_LINE];
    FILE *fp;

    snprintf(path, MAX_PATH, "%s/%s/Header", plotfile_dir, SDM_SUBDIR);
    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", path);
        return -1;
    }

    /* Line 1: version string */
    fgets(line, MAX_LINE, fp);
    line[strcspn(line, "\n")] = 0;
    if (strstr(line, "Version_Two") == NULL) {
        fprintf(stderr, "Warning: Unexpected particle version: %s\n", line);
    }

    /* Line 2: ndim */
    fgets(line, MAX_LINE, fp);
    pd->ndim = atoi(line);

    /* Line 3: n_real_comps (excluding x,y,z) */
    fgets(line, MAX_LINE, fp);
    pd->n_real_comps = atoi(line);

    /* Real component names */
    for (int i = 0; i < pd->n_real_comps && i < MAX_SDM_VARS; i++) {
        fgets(line, MAX_LINE, fp);
        line[strcspn(line, "\n")] = 0;
        strncpy(pd->real_comp_names[i], line, 63);
    }

    /* n_int_comps (excluding id, cpu) */
    fgets(line, MAX_LINE, fp);
    pd->n_int_comps = atoi(line);

    /* Int component names */
    for (int i = 0; i < pd->n_int_comps && i < MAX_SDM_VARS; i++) {
        fgets(line, MAX_LINE, fp);
        line[strcspn(line, "\n")] = 0;
        strncpy(pd->int_comp_names[i], line, 63);
    }

    /* Skip: is_checkpoint line */
    fgets(line, MAX_LINE, fp);

    /* Total number of particles */
    fgets(line, MAX_LINE, fp);
    pd->n_particles = atoi(line);

    /* Skip: max_next_id */
    fgets(line, MAX_LINE, fp);

    /* finest_level (should be 0 for single-level) */
    fgets(line, MAX_LINE, fp);
    /* int finest_level = atoi(line); */

    /* Number of grids at level 0 */
    fgets(line, MAX_LINE, fp);
    pd->n_grids = atoi(line);

    /* Per-grid info: file_number count offset */
    for (int i = 0; i < pd->n_grids && i < MAX_BOXES; i++) {
        fgets(line, MAX_LINE, fp);
        sscanf(line, "%d %d %ld", &pd->grid_file_num[i],
               &pd->grid_count[i], &pd->grid_offset[i]);
    }

    fclose(fp);

    /* Find indices for radius, multiplicity, particle_mass in real comp names */
    pd->radius_idx = -1;
    pd->mult_idx = -1;
    pd->mass_idx = -1;

    for (int i = 0; i < pd->n_real_comps; i++) {
        if (strcmp(pd->real_comp_names[i], "radius") == 0)
            pd->radius_idx = i;
        else if (strcmp(pd->real_comp_names[i], "multiplicity") == 0)
            pd->mult_idx = i;
        else if (strcmp(pd->real_comp_names[i], "particle_mass") == 0)
            pd->mass_idx = i;
    }

    if (pd->radius_idx < 0 || pd->mult_idx < 0 || pd->mass_idx < 0) {
        fprintf(stderr, "Error: Missing required particle components (radius=%d, multiplicity=%d, particle_mass=%d)\n",
                pd->radius_idx, pd->mult_idx, pd->mass_idx);
        return -1;
    }

    printf("SDM Header: %d particles, %d real comps, %d int comps, %d grids\n",
           pd->n_particles, pd->n_real_comps, pd->n_int_comps, pd->n_grids);
    printf("  radius_idx=%d, multiplicity_idx=%d, mass_idx=%d\n",
           pd->radius_idx, pd->mult_idx, pd->mass_idx);

    return 0;
}

/* Compute domain volume from main plotfile Header */
double compute_domain_volume(const char *plotfile_dir) {
    char path[MAX_PATH];
    char line[MAX_LINE];
    FILE *fp;

    snprintf(path, MAX_PATH, "%s/Header", plotfile_dir);
    fp = fopen(path, "r");
    if (!fp) return 1.0;

    /* Skip: version */
    fgets(line, MAX_LINE, fp);
    /* Skip: n_vars */
    fgets(line, MAX_LINE, fp);
    int n_vars = atoi(line);
    /* Skip: variable names */
    for (int i = 0; i < n_vars; i++)
        fgets(line, MAX_LINE, fp);
    /* Read: ndim */
    fgets(line, MAX_LINE, fp);
    int ndim = atoi(line);
    /* Skip: time */
    fgets(line, MAX_LINE, fp);
    /* Skip: finest_level */
    fgets(line, MAX_LINE, fp);

    /* Read prob_lo (one value per dimension on one line) */
    fgets(line, MAX_LINE, fp);
    double prob_lo[3] = {0, 0, 0};
    if (ndim == 3) sscanf(line, "%lf %lf %lf", &prob_lo[0], &prob_lo[1], &prob_lo[2]);
    else if (ndim == 2) sscanf(line, "%lf %lf", &prob_lo[0], &prob_lo[1]);

    /* Read prob_hi */
    fgets(line, MAX_LINE, fp);
    double prob_hi[3] = {0, 0, 0};
    if (ndim == 3) sscanf(line, "%lf %lf %lf", &prob_hi[0], &prob_hi[1], &prob_hi[2]);
    else if (ndim == 2) sscanf(line, "%lf %lf", &prob_hi[0], &prob_hi[1]);

    fclose(fp);

    double volume = 1.0;
    for (int d = 0; d < ndim; d++) {
        volume *= (prob_hi[d] - prob_lo[d]);
    }

    printf("Domain volume: %g (bounds: [%g,%g] x [%g,%g] x [%g,%g])\n",
           volume, prob_lo[0], prob_hi[0], prob_lo[1], prob_hi[1], prob_lo[2], prob_hi[2]);
    return volume;
}

/* Read particle binary data from DATA files */
int read_sdm_data(ParticleData *pd, const char *plotfile_dir) {
    char path[MAX_PATH];

    /* Free previous data */
    if (pd->radius) { free(pd->radius); pd->radius = NULL; }
    if (pd->multiplicity) { free(pd->multiplicity); pd->multiplicity = NULL; }
    if (pd->mass) { free(pd->mass); pd->mass = NULL; }

    if (pd->n_particles <= 0) {
        printf("No particles in %s (0 particles)\n", plotfile_dir);
        return 0;  /* Not an error — timestep may simply have no particles yet */
    }

    pd->radius = (double *)malloc(pd->n_particles * sizeof(double));
    pd->multiplicity = (double *)malloc(pd->n_particles * sizeof(double));
    pd->mass = (double *)malloc(pd->n_particles * sizeof(double));

    int ints_per_particle = 2 + pd->n_int_comps;   /* id, cpu, int_comp0, int_comp1, ... */
    int reals_per_particle = pd->ndim + pd->n_real_comps;  /* x,y,z + real comps */

    /* Indices within the real block (0-based, including x,y,z) */
    int real_radius_idx = pd->ndim + pd->radius_idx;
    int real_mult_idx = pd->ndim + pd->mult_idx;
    int real_mass_idx = pd->ndim + pd->mass_idx;

    int particle_offset = 0;  /* Running offset into output arrays */

    for (int g = 0; g < pd->n_grids; g++) {
        int count = pd->grid_count[g];
        if (count <= 0) continue;

        snprintf(path, MAX_PATH, "%s/%s/Level_0/DATA_%05d",
                 plotfile_dir, SDM_SUBDIR, pd->grid_file_num[g]);
        FILE *fp = fopen(path, "rb");
        if (!fp) {
            fprintf(stderr, "Error: Cannot open %s\n", path);
            continue;
        }

        /* Seek to grid offset */
        fseek(fp, pd->grid_offset[g], SEEK_SET);

        /* Skip int block: count * ints_per_particle * sizeof(int32_t) */
        fseek(fp, (long)count * ints_per_particle * sizeof(int), SEEK_CUR);

        /* Read real data for all particles in this grid */
        double *real_buf = (double *)malloc((size_t)count * reals_per_particle * sizeof(double));
        size_t read_count = fread(real_buf, sizeof(double), (size_t)count * reals_per_particle, fp);
        fclose(fp);

        if ((int)read_count != count * reals_per_particle) {
            fprintf(stderr, "Warning: Short read for grid %d: got %zu expected %d\n",
                    g, read_count, count * reals_per_particle);
        }

        /* Extract radius, multiplicity, mass for each particle */
        for (int p = 0; p < count && (particle_offset + p) < pd->n_particles; p++) {
            double *pdata = real_buf + (size_t)p * reals_per_particle;
            pd->radius[particle_offset + p] = pdata[real_radius_idx];
            pd->multiplicity[particle_offset + p] = pdata[real_mult_idx];
            pd->mass[particle_offset + p] = pdata[real_mass_idx];
        }

        particle_offset += count;
        free(real_buf);
    }

    printf("Loaded %d particles from %s\n", particle_offset, plotfile_dir);
    return 0;
}

/* Extract 2D slice from 3D data */
void extract_slice(PlotfileData *pf, double *slice, int axis, int idx) {
    int i, j, k;
    int nx = pf->grid_dims[0];
    int ny = pf->grid_dims[1];
    int nz = pf->grid_dims[2];

    if (axis == 2) {  /* Z slice */
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                slice[j * nx + i] = pf->data[idx * ny * nx + j * nx + i];
            }
        }
    } else if (axis == 1) {  /* Y slice */
        for (k = 0; k < nz; k++) {
            for (i = 0; i < nx; i++) {
                slice[k * nx + i] = pf->data[k * ny * nx + idx * nx + i];
            }
        }
    } else {  /* X slice */
        for (k = 0; k < nz; k++) {
            for (j = 0; j < ny; j++) {
                slice[k * ny + j] = pf->data[k * ny * nx + j * nx + idx];
            }
        }
    }
}

/* Extract slice from a specific level's data */
void extract_slice_level(LevelData *ld, double *slice, int axis, int idx) {
    int i, j, k;
    int nx = ld->grid_dims[0];
    int ny = ld->grid_dims[1];
    int nz = ld->grid_dims[2];

    if (axis == 2) {  /* Z slice */
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                slice[j * nx + i] = ld->data[idx * ny * nx + j * nx + i];
            }
        }
    } else if (axis == 1) {  /* Y slice */
        for (k = 0; k < nz; k++) {
            for (i = 0; i < nx; i++) {
                slice[k * nx + i] = ld->data[k * ny * nx + idx * nx + i];
            }
        }
    } else {  /* X slice */
        for (k = 0; k < nz; k++) {
            for (j = 0; j < ny; j++) {
                slice[k * ny + j] = ld->data[k * ny * nx + j * nx + idx];
            }
        }
    }
}

/* Drop the cached terrain coordinate when the plotfile or AMR level changes. */
void free_z_phys_cache(PlotfileData *pf) {
    if (!pf) return;
    free(pf->z_phys_data);
    pf->z_phys_data = NULL;
    pf->z_phys_level = -1;
}

/* Load the 3D z_phys coordinate without replacing the selected scalar field. */
int ensure_z_phys_data(PlotfileData *pf) {
    if (!pf) return -1;
    int z_idx = find_variable_index(pf, "z_phys");
    if (z_idx < 0) return -1;
    if (pf->z_phys_data && pf->z_phys_level == pf->current_level) return 0;

    free_z_phys_cache(pf);
    size_t count = (size_t)pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2];

    if (z_idx == pf->current_var && pf->data) {
        pf->z_phys_data = (double *)malloc(count * sizeof(double));
        if (!pf->z_phys_data) return -1;
        memcpy(pf->z_phys_data, pf->data, count * sizeof(double));
    } else {
        double *saved_data = pf->data;
        pf->data = NULL;
        if (read_variable_data(pf, z_idx) < 0 || !pf->data) {
            free(pf->data);
            pf->data = saved_data;
            return -1;
        }
        pf->z_phys_data = pf->data;
        pf->data = saved_data;
    }

    pf->z_phys_level = pf->current_level;
    return 0;
}

/* Load z_phys for an overlay level while keeping that level's scalar data intact. */
int ensure_z_phys_level_data(PlotfileData *pf, int level) {
    if (!pf || level < 0 || level >= pf->n_levels || level >= MAX_LEVELS) return -1;
    int z_idx = find_variable_index(pf, "z_phys");
    if (z_idx < 0) return -1;

    LevelData *ld = &pf->levels[level];
    if (ld->z_phys_loaded && ld->z_phys_data) return 0;

    size_t count = (size_t)ld->grid_dims[0] * ld->grid_dims[1] * ld->grid_dims[2];
    if (z_idx == pf->current_var && ld->data) {
        ld->z_phys_data = (double *)malloc(count * sizeof(double));
        if (!ld->z_phys_data) return -1;
        memcpy(ld->z_phys_data, ld->data, count * sizeof(double));
    } else {
        double *saved_data = ld->data;
        int saved_loaded = ld->loaded;
        ld->data = NULL;
        if (read_variable_data_level(pf, z_idx, level) < 0 || !ld->data) {
            free(ld->data);
            ld->data = saved_data;
            ld->loaded = saved_loaded;
            return -1;
        }
        ld->z_phys_data = ld->data;
        ld->data = saved_data;
        ld->loaded = saved_loaded;
    }

    ld->z_phys_loaded = 1;
    return 0;
}

static double z_phys_vertical_interface(const double *z_values, int width, int height,
                                        int column, int boundary) {
    if (!z_values || width <= 0 || height <= 0 || column < 0 || column >= width)
        return NAN;
    if (height == 1) return z_values[column];
    if (boundary <= 0) {
        double z0 = z_values[column];
        double z1 = z_values[width + column];
        return z0 - 0.5 * (z1 - z0);
    }
    if (boundary >= height) {
        double z0 = z_values[(height - 2) * width + column];
        double z1 = z_values[(height - 1) * width + column];
        return z1 + 0.5 * (z1 - z0);
    }
    double below = z_values[(boundary - 1) * width + column];
    double above = z_values[boundary * width + column];
    return 0.5 * (below + above);
}

/* Cell-corner coordinate reconstructed from cell-centered, 3D z_phys values. */
double z_phys_corner(const double *z_values, int width, int height,
                     int x_boundary, int z_boundary) {
    if (x_boundary <= 0)
        return z_phys_vertical_interface(z_values, width, height, 0, z_boundary);
    if (x_boundary >= width)
        return z_phys_vertical_interface(z_values, width, height, width - 1, z_boundary);

    double left = z_phys_vertical_interface(z_values, width, height,
                                            x_boundary - 1, z_boundary);
    double right = z_phys_vertical_interface(z_values, width, height,
                                             x_boundary, z_boundary);
    if (!isfinite(left)) return right;
    if (!isfinite(right)) return left;
    return 0.5 * (left + right);
}

static short z_phys_screen_coord(double value) {
    if (value < -32760.0) return -32760;
    if (value > 32760.0) return 32760;
    return (short)lrint(value);
}

/* Draw an X/Z or Y/Z cross-section in its terrain-following physical geometry. */
void draw_z_phys_cells(const double *z_values, const unsigned long *pixels,
                       const unsigned char *mask, int width, int height,
                       double data_xmin, double data_xmax,
                       double view_xmin, double view_xmax,
                       double view_ymin, double view_ymax,
                       int offset_x, int offset_y, int draw_width, int draw_height) {
    if (!z_values || !pixels || width <= 0 || height <= 0 ||
        view_xmax <= view_xmin || view_ymax <= view_ymin) return;

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            int idx = j * width + i;
            if (mask && !mask[idx]) continue;

            double z00 = z_phys_corner(z_values, width, height, i, j);
            double z10 = z_phys_corner(z_values, width, height, i + 1, j);
            double z11 = z_phys_corner(z_values, width, height, i + 1, j + 1);
            double z01 = z_phys_corner(z_values, width, height, i, j + 1);
            if (!isfinite(z00) || !isfinite(z10) || !isfinite(z11) || !isfinite(z01))
                continue;

            double x0 = data_xmin + (double)i / width * (data_xmax - data_xmin);
            double x1 = data_xmin + (double)(i + 1) / width * (data_xmax - data_xmin);
            double sx0 = offset_x + (x0 - view_xmin) / (view_xmax - view_xmin) * draw_width;
            double sx1 = offset_x + (x1 - view_xmin) / (view_xmax - view_xmin) * draw_width;
            double sy00 = offset_y + (view_ymax - z00) / (view_ymax - view_ymin) * draw_height;
            double sy10 = offset_y + (view_ymax - z10) / (view_ymax - view_ymin) * draw_height;
            double sy11 = offset_y + (view_ymax - z11) / (view_ymax - view_ymin) * draw_height;
            double sy01 = offset_y + (view_ymax - z01) / (view_ymax - view_ymin) * draw_height;

            XPoint points[4] = {
                {z_phys_screen_coord(sx0), z_phys_screen_coord(sy00)},
                {z_phys_screen_coord(sx1), z_phys_screen_coord(sy10)},
                {z_phys_screen_coord(sx1), z_phys_screen_coord(sy11)},
                {z_phys_screen_coord(sx0), z_phys_screen_coord(sy01)}
            };
            XSetForeground(display, gc, pixels[idx]);
            XFillPolygon(display, canvas, gc, points, 4, Convex, CoordModeOrigin);
        }
    }
}

/* Build a fast inverse map from displayed map pixels back to scalar cells. */
int prepare_map_hover_lookup(void) {
    size_t needed = (size_t)canvas_width * canvas_height;
    if (needed == 0) return -1;
    if (needed > map_hover_lookup_size) {
        int *new_lookup = (int *)realloc(map_hover_lookup, needed * sizeof(int));
        if (!new_lookup) {
            map_hover_lookup_active = 0;
            return -1;
        }
        map_hover_lookup = new_lookup;
        map_hover_lookup_size = needed;
    }
    for (size_t p = 0; p < needed; p++) map_hover_lookup[p] = -1;
    map_hover_lookup_active = 1;
    return 0;
}

/* Record the same clipped rectangle that represents a map cell on screen. */
void record_map_hover_rect(int x, int y, int width, int height, int cell_idx) {
    if (!map_hover_lookup_active || !map_hover_lookup || width <= 0 || height <= 0)
        return;

    int x0 = x;
    int y0 = y;
    int x1 = x + width;
    int y1 = y + height;
    if (x0 < vis_area_x) x0 = vis_area_x;
    if (y0 < vis_area_y) y0 = vis_area_y;
    if (x1 > vis_area_x + vis_area_w) x1 = vis_area_x + vis_area_w;
    if (y1 > vis_area_y + vis_area_h) y1 = vis_area_y + vis_area_h;
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 > canvas_width) x1 = canvas_width;
    if (y1 > canvas_height) y1 = canvas_height;
    if (x0 >= x1 || y0 >= y1) return;

    for (int py = y0; py < y1; py++) {
        int *row = map_hover_lookup + (size_t)py * canvas_width;
        for (int px = x0; px < x1; px++) row[px] = cell_idx;
    }
}

/* Convert a canvas location back to a slice cell, including terrain geometry. */
int canvas_to_data_indices(int mouse_x, int mouse_y, int *data_x, int *data_y) {
    if (!data_x || !data_y || render_width <= 0 || render_height <= 0 ||
        slice_width <= 0 || slice_height <= 0) return 0;

    if (map_hover_lookup_active && map_hover_lookup &&
        mouse_x >= 0 && mouse_x < canvas_width && mouse_y >= 0 && mouse_y < canvas_height) {
        int idx = map_hover_lookup[(size_t)mouse_y * canvas_width + mouse_x];
        if (idx < 0 || idx >= slice_width * slice_height) return 0;
        *data_x = idx % slice_width;
        *data_y = idx / slice_width;
        return 1;
    }
    if (global_pf && global_pf->map_mode) return 0;

    int ix = (int)((mouse_x - render_offset_x) * slice_width / (double)render_width);
    if (ix < 0 || ix >= slice_width) return 0;

    if (!render_uses_z_phys || !current_z_phys_slice) {
        int iy = slice_height - 1 -
                 (int)((mouse_y - render_offset_y) * slice_height / (double)render_height);
        if (iy < 0 || iy >= slice_height) return 0;
        *data_x = ix;
        *data_y = iy;
        return 1;
    }

    double z = render_phys_ymax - (mouse_y - render_offset_y) /
               (double)render_height * (render_phys_ymax - render_phys_ymin);
    for (int j = 0; j < slice_height; j++) {
        double z0 = z_phys_vertical_interface(current_z_phys_slice, slice_width,
                                              slice_height, ix, j);
        double z1 = z_phys_vertical_interface(current_z_phys_slice, slice_width,
                                              slice_height, ix, j + 1);
        if (!isfinite(z0) || !isfinite(z1)) continue;
        if (z >= fmin(z0, z1) && z <= fmax(z0, z1)) {
            *data_x = ix;
            *data_y = j;
            return 1;
        }
    }
    return 0;
}

/* Jet colormap */
RGB jet_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    if (t < 0.25) {
        color.r = 0;
        color.g = (unsigned char)(255 * (4 * t));
        color.b = 255;
    } else if (t < 0.5) {
        color.r = 0;
        color.g = 255;
        color.b = (unsigned char)(255 * (1 - 4 * (t - 0.25)));
    } else if (t < 0.75) {
        color.r = (unsigned char)(255 * (4 * (t - 0.5)));
        color.g = 255;
        color.b = 0;
    } else {
        color.r = 255;
        color.g = (unsigned char)(255 * (1 - 4 * (t - 0.75)));
        color.b = 0;
    }
    return color;
}

/* Turbo colormap (approximation) */
RGB turbo_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    double r = t * 0.8 + 0.2;
    double g = sin(t * 3.14159);
    double b = 1.0 - t * 0.9;
    
    color.r = (unsigned char)(255 * r);
    color.g = (unsigned char)(255 * g);
    color.b = (unsigned char)(255 * b);
    return color;
}

/* Plasma colormap */
RGB plasma_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    if (t < 0.5) {
        color.r = (unsigned char)(13 + (177 - 13) * (t / 0.5));
        color.g = (unsigned char)(8 + (42 - 8) * (t / 0.5));
        color.b = (unsigned char)(135 + (127 - 135) * (t / 0.5));
    } else {
        color.r = (unsigned char)(177 + (240 - 177) * ((t - 0.5) / 0.5));
        color.g = (unsigned char)(42 + (249 - 42) * ((t - 0.5) / 0.5));
        color.b = (unsigned char)(127 + (33 - 127) * ((t - 0.5) / 0.5));
    }
    return color;
}

/* Viridis colormap */
RGB viridis_colormap(double t) {
    RGB color;
    /* Simplified viridis approximation */
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    if (t < 0.25) {
        color.r = (unsigned char)(68 + (59 - 68) * (t / 0.25));
        color.g = (unsigned char)(1 + (82 - 1) * (t / 0.25));
        color.b = (unsigned char)(84 + (139 - 84) * (t / 0.25));
    } else if (t < 0.5) {
        color.r = (unsigned char)(59 + (33 - 59) * ((t - 0.25) / 0.25));
        color.g = (unsigned char)(82 + (144 - 82) * ((t - 0.25) / 0.25));
        color.b = (unsigned char)(139 + (140 - 139) * ((t - 0.25) / 0.25));
    } else if (t < 0.75) {
        color.r = (unsigned char)(33 + (93 - 33) * ((t - 0.5) / 0.25));
        color.g = (unsigned char)(144 + (201 - 144) * ((t - 0.5) / 0.25));
        color.b = (unsigned char)(140 + (99 - 140) * ((t - 0.5) / 0.25));
    } else {
        color.r = (unsigned char)(93 + (253 - 93) * ((t - 0.75) / 0.25));
        color.g = (unsigned char)(201 + (231 - 201) * ((t - 0.75) / 0.25));
        color.b = (unsigned char)(99 + (37 - 99) * ((t - 0.75) / 0.25));
    }
    
    return color;
}

/* Hot colormap (black -> red -> yellow -> white) */
RGB hot_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    if (t < 0.33) {
        color.r = (unsigned char)(255 * (t / 0.33));
        color.g = 0;
        color.b = 0;
    } else if (t < 0.67) {
        color.r = 255;
        color.g = (unsigned char)(255 * ((t - 0.33) / 0.34));
        color.b = 0;
    } else {
        color.r = 255;
        color.g = 255;
        color.b = (unsigned char)(255 * ((t - 0.67) / 0.33));
    }
    return color;
}

/* Cool colormap (cyan -> blue -> magenta) */
RGB cool_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    color.r = (unsigned char)(255 * t);
    color.g = (unsigned char)(255 * (1.0 - t));
    color.b = 255;
    return color;
}

/* Gray colormap (black -> white) */
RGB gray_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    unsigned char val = (unsigned char)(255 * t);
    color.r = val;
    color.g = val;
    color.b = val;
    return color;
}

/* Magma colormap (black -> purple -> orange -> white) */
RGB magma_colormap(double t) {
    RGB color;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    
    if (t < 0.25) {
        color.r = (unsigned char)(8 + (72 - 8) * (t / 0.25));
        color.g = (unsigned char)(8 + (22 - 8) * (t / 0.25));
        color.b = (unsigned char)(40 + (84 - 40) * (t / 0.25));
    } else if (t < 0.5) {
        color.r = (unsigned char)(72 + (161 - 72) * ((t - 0.25) / 0.25));
        color.g = (unsigned char)(22 + (51 - 22) * ((t - 0.25) / 0.25));
        color.b = (unsigned char)(84 + (118 - 84) * ((t - 0.25) / 0.25));
    } else if (t < 0.75) {
        color.r = (unsigned char)(161 + (235 - 161) * ((t - 0.5) / 0.25));
        color.g = (unsigned char)(51 + (105 - 51) * ((t - 0.5) / 0.25));
        color.b = (unsigned char)(118 + (81 - 118) * ((t - 0.5) / 0.25));
    } else {
        color.r = (unsigned char)(235 + (252 - 235) * ((t - 0.75) / 0.25));
        color.g = (unsigned char)(105 + (191 - 105) * ((t - 0.75) / 0.25));
        color.b = (unsigned char)(81 + (170 - 81) * ((t - 0.75) / 0.25));
    }
    return color;
}

/* Get RGB for any colormap */
RGB get_colormap_rgb(double t, int cmap_type) {
    switch(cmap_type) {
        case 1: return jet_colormap(t);
        case 2: return turbo_colormap(t);
        case 3: return plasma_colormap(t);
        case 4: return hot_colormap(t);
        case 5: return cool_colormap(t);
        case 6: return gray_colormap(t);
        case 7: return magma_colormap(t);
        default: return viridis_colormap(t);
    }
}

/* Get unit string for a variable based on common naming conventions */
const char *get_variable_unit(const char *varname) {
    if (!varname) return "";

    /* Velocity components */
    if (strstr(varname, "velocity") || strstr(varname, "vel_") ||
        strcmp(varname, "u") == 0 || strcmp(varname, "v") == 0 ||
        strcmp(varname, "w") == 0) {
        return "m/s";
    }
    /* Temperature and potential temperature */
    if (strstr(varname, "temp") || strstr(varname, "theta") ||
        strcmp(varname, "T") == 0) {
        return "K";
    }
    /* Pressure */
    if (strstr(varname, "pressure") || strstr(varname, "pres") ||
        strcmp(varname, "p") == 0 || strcmp(varname, "P") == 0) {
        return "Pa";
    }
    /* Density */
    if (strcmp(varname, "density") == 0 || strcmp(varname, "rho") == 0) {
        return "kg/m^3";
    }
    /* Density times theta */
    if (strstr(varname, "rhotheta")) {
        return "kg K/m^3";
    }
    /* Mixing ratios */
    if (strncmp(varname, "q", 1) == 0 && strlen(varname) <= 6) {
        return "kg/kg";
    }
    /* Relative humidity */
    if (strstr(varname, "humidity") || strstr(varname, "rh") ||
        strcmp(varname, "RH") == 0) {
        return "";  /* fraction or % depending on simulation */
    }
    /* Number density */
    if (strstr(varname, "number_density")) {
        return "1/m^3";
    }
    /* Mass density (for particles) */
    if (strstr(varname, "mass_density")) {
        return "kg/m^3";
    }
    /* Radius (droplets) */
    if (strstr(varname, "radius")) {
        return "m";
    }
    /* Vorticity */
    if (strstr(varname, "vort") || strstr(varname, "omega")) {
        return "1/s";
    }
    /* TKE */
    if (strstr(varname, "tke") || strstr(varname, "TKE")) {
        return "m^2/s^2";
    }

    return "";  /* Unknown - no unit */
}

/* Apply colormap to data with scale_mode support */
void apply_colormap(double *data, int width, int height,
                   unsigned long *pixels, double vmin, double vmax, int cmap_type) {
    int i, j;
    double range = vmax - vmin;
    if (range < 1e-10) range = 1.0;

    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            double val = data[j * width + i];
            double t;

            if (scale_mode == 0) {  /* Linear */
                t = (val - vmin) / range;
            } else if (scale_mode == 1) {  /* Log(+) */
                if (val <= 0) {
                    pixels[j * width + i] = 0xC0C0C0;  /* Gray for non-positive */
                    continue;
                }
                /* Map log(val) to [0,1] using log(vmin) to log(vmax) */
                double log_vmin = (vmin > 0) ? log10(vmin) : log10(1e-10);
                double log_vmax = (vmax > 0) ? log10(vmax) : log10(1e-10);
                if (log_vmax - log_vmin < 1e-10) {
                    t = 0.5;
                } else {
                    t = (log10(val) - log_vmin) / (log_vmax - log_vmin);
                }
            } else {  /* Log(-) for scale_mode == 2 */
                if (val >= 0) {
                    pixels[j * width + i] = 0xC0C0C0;  /* Gray for non-negative */
                    continue;
                }
                /* Use absolute value, flip the range */
                double abs_val = -val;
                double abs_vmax = (vmin < 0) ? -vmin : 1e-10;
                double abs_vmin = (vmax < 0) ? -vmax : 1e-10;
                double log_vmin = log10(abs_vmin > 0 ? abs_vmin : 1e-10);
                double log_vmax = log10(abs_vmax > 0 ? abs_vmax : 1e-10);
                if (log_vmax - log_vmin < 1e-10) {
                    t = 0.5;
                } else {
                    t = (log10(abs_val) - log_vmin) / (log_vmax - log_vmin);
                }
            }

            if (t < 0) t = 0;
            if (t > 1) t = 1;
            RGB color = get_colormap_rgb(t, cmap_type);
            pixels[j * width + i] = (color.r << 16) | (color.g << 8) | color.b;
        }
    }
}

/* Draw colorbar with variable name and units */
void draw_colorbar(double vmin, double vmax, int cmap_type, const char *varname) {
    int height = 256, width = 30;
    int top_margin = 50;   /* Extra space at top for variable name */
    int bottom_margin = 10;

    /* Clear colorbar with white background */
    XSetForeground(display, colorbar_gc, WhitePixel(display, screen));
    XFillRectangle(display, colorbar, colorbar_gc, 0, 0, 100, canvas_height);

    /* Draw variable name at top */
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    if (varname && strlen(varname) > 0) {
        /* Truncate long variable names */
        char short_name[16];
        if (strlen(varname) > 12) {
            strncpy(short_name, varname, 11);
            short_name[11] = '\0';
            strcat(short_name, "..");
        } else {
            strncpy(short_name, varname, 15);
            short_name[15] = '\0';
        }
        XDrawString(display, colorbar, text_gc, 2, 15, short_name, strlen(short_name));

        /* Draw unit below variable name */
        const char *unit = get_variable_unit(varname);
        if (unit && strlen(unit) > 0) {
            char unit_str[20];
            snprintf(unit_str, sizeof(unit_str), "[%s]", unit);
            XDrawString(display, colorbar, text_gc, 2, 30, unit_str, strlen(unit_str));
        }
    }

    /* Draw colorbar as solid rectangles within margins */
    for (int i = 0; i < height; i++) {
        double t = (double)(height - 1 - i) / (height - 1);
        RGB color = get_colormap_rgb(t, cmap_type);
        unsigned long pixel = (color.r << 16) | (color.g << 8) | color.b;

        XSetForeground(display, colorbar_gc, pixel);
        int y = top_margin + (i * (canvas_height - top_margin - bottom_margin)) / height;
        int h = top_margin + ((i + 1) * (canvas_height - top_margin - bottom_margin)) / height - y;
        if (h < 1) h = 1;
        XFillRectangle(display, colorbar, colorbar_gc, 0, y, width, h);
    }

    /* Draw tick marks and labels */
    char text[32];
    XSetForeground(display, text_gc, BlackPixel(display, screen));

    int n_ticks = 11;  /* 11 ticks gives 10 intervals */

    for (int i = 0; i < n_ticks; i++) {
        double fraction = (double)i / (n_ticks - 1);
        double value;

        if (scale_mode == 0) {  /* Linear */
            value = vmin + fraction * (vmax - vmin);
            snprintf(text, sizeof(text), "%.2e", value);
        } else if (scale_mode == 1) {  /* Log(+) */
            /* Logarithmic tick spacing */
            double log_vmin = (vmin > 0) ? log10(vmin) : log10(1e-10);
            double log_vmax = (vmax > 0) ? log10(vmax) : log10(1e-10);
            double log_val = log_vmin + fraction * (log_vmax - log_vmin);
            value = pow(10.0, log_val);
            /* Format as 10^x for cleaner display */
            if (log_val == floor(log_val)) {
                snprintf(text, sizeof(text), "10^%.0f", log_val);
            } else {
                snprintf(text, sizeof(text), "%.1e", value);
            }
        } else {  /* Log(-) */
            /* For negative log scale, show -10^x format */
            double abs_vmax = (vmin < 0) ? -vmin : 1e-10;
            double abs_vmin = (vmax < 0) ? -vmax : 1e-10;
            double log_vmin = log10(abs_vmin > 0 ? abs_vmin : 1e-10);
            double log_vmax = log10(abs_vmax > 0 ? abs_vmax : 1e-10);
            double log_val = log_vmin + fraction * (log_vmax - log_vmin);
            value = -pow(10.0, log_val);
            /* Format as -10^x for cleaner display */
            if (log_val == floor(log_val)) {
                snprintf(text, sizeof(text), "-10^%.0f", log_val);
            } else {
                snprintf(text, sizeof(text), "%.1e", value);
            }
        }

        /* Map to drawable area with margins */
        int y = top_margin + (canvas_height - top_margin - bottom_margin) -
                (int)(fraction * (canvas_height - top_margin - bottom_margin));

        /* Draw tick mark */
        XDrawLine(display, colorbar, text_gc, width, y, width + 5, y);

        /* Draw label with vertical centering adjustment */
        XDrawString(display, colorbar, text_gc, width + 8, y + 4, text, strlen(text));
    }

    XFlush(display);
}

/* Initialize GUI with Athena Widgets */
void init_gui(PlotfileData *pf, int argc, char **argv) {
    Arg args[20];
    int n, i;
    Widget button, label;
    char label_text[64];
    
    global_pf = pf;
    
    toplevel = XtAppInitialize(NULL, "PLTView", NULL, 0, &argc, argv, NULL, NULL, 0);
    display = XtDisplay(toplevel);
    screen = DefaultScreen(display);
    
    /* Load font */
    font = XLoadQueryFont(display, "fixed");
    if (!font) font = XLoadQueryFont(display, "*");
    
    /* Main form container */
    n = 0;
    XtSetArg(args[n], XtNwidth, 1000); n++;
    XtSetArg(args[n], XtNheight, 700); n++;
    form = XtCreateManagedWidget("form", formWidgetClass, toplevel, args, n);
    
    /* Info label at top */
    snprintf(label_text, sizeof(label_text), "PLTView - Loading...");
    n = 0;
    XtSetArg(args[n], XtNlabel, label_text); n++;
    XtSetArg(args[n], XtNwidth, 900); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    info_label = XtCreateManagedWidget("info", labelWidgetClass, form, args, n);
    
    /* Calculate viewport width based on longest variable name */
    int var_viewport_width = 120;  /* Minimum width */
    if (font) {
        int max_text_width = 0;
        for (i = 0; i < pf->n_vars; i++) {
            int w = XTextWidth(font, pf->variables[i], strlen(pf->variables[i]));
            if (w > max_text_width) max_text_width = w;
        }
        /* Add padding: button internalWidth (8) + border (4) + highlight (4) +
           box hSpace (8) + scrollbar (20) + safety margin (36) = 80 */
        var_viewport_width = max_text_width + 80;
        if (var_viewport_width < 120) var_viewport_width = 120;
        if (var_viewport_width > 400) var_viewport_width = 400;
    }

    /* Variable buttons viewport (use internal scrollbar) */
    n = 0;
    XtSetArg(args[n], XtNfromVert, info_label); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNallowVert, True); n++;
    XtSetArg(args[n], XtNallowHoriz, False); n++;
    XtSetArg(args[n], XtNforceBars, True); n++;
    XtSetArg(args[n], XtNheight, canvas_height); n++;
    XtSetArg(args[n], XtNwidth, var_viewport_width); n++;
    var_viewport = XtCreateManagedWidget("varViewport", viewportWidgetClass, form, args, n);
    XtAddEventHandler(var_viewport, ButtonPressMask, False, var_viewport_wheel_handler, NULL);
    XtAddEventHandler(var_viewport, StructureNotifyMask, False, var_viewport_configure_handler, NULL);

    /* Variable buttons box inside viewport */
    n = 0;
    XtSetArg(args[n], XtNorientation, XtorientVertical); n++;
    var_box = XtCreateManagedWidget("varBox", boxWidgetClass, var_viewport, args, n);
    
    /* Add variable buttons (all available variables) */
    for (i = 0; i < pf->n_vars; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, pf->variables[i]); n++;
        button = XtCreateManagedWidget(pf->variables[i], commandWidgetClass, var_box, args, n);
        XtAddCallback(button, XtNcallback, var_button_callback, (XtPointer)(long)i);
    }
    update_var_scrollbar();
    
    /* Canvas drawing area */
    n = 0;
    XtSetArg(args[n], XtNfromVert, info_label); n++;
    XtSetArg(args[n], XtNfromHoriz, var_viewport); n++;
    XtSetArg(args[n], XtNwidth, canvas_width); n++;
    XtSetArg(args[n], XtNheight, canvas_height); n++;
    XtSetArg(args[n], XtNborderWidth, 2); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    canvas_widget = XtCreateManagedWidget("canvas", simpleWidgetClass, form, args, n);
    XtAddCallback(canvas_widget, XtNcallback, canvas_expose_callback, NULL);
    
    /* COLUMN 1: Navigation buttons */
    /* Navigation buttons (+/-) in Column 1, Row 1 */
    n = 0;
    XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
    XtSetArg(args[n], XtNfromHoriz, var_viewport); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    nav_box = XtCreateManagedWidget("navBox", boxWidgetClass, form, args, n);
    
    /* Layer label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Layer"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("layerText", labelWidgetClass, nav_box, args, n);

    /* Navigation buttons (v/^) */
    const char *nav_labels[] = {"v", "^"};
    for (i = 0; i < 2; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, nav_labels[i]); n++;
        button = XtCreateManagedWidget(nav_labels[i], commandWidgetClass, nav_box, args, n);
        XtAddCallback(button, XtNcallback, nav_button_callback, (XtPointer)(long)i);
    }

    /* Layer index display label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "1/1"); n++;
    XtSetArg(args[n], XtNwidth, 60); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    layer_label = XtCreateManagedWidget("layerLabel", labelWidgetClass, nav_box, args, n);
    
    /* Jump button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Jump"); n++;
    button = XtCreateManagedWidget("jump", commandWidgetClass, nav_box, args, n);
    XtAddCallback(button, XtNcallback, jump_button_callback, NULL);

    /* Profile button for slice statistics */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Profile"); n++;
    button = XtCreateManagedWidget("profile", commandWidgetClass, nav_box, args, n);
    XtAddCallback(button, XtNcallback, profile_button_callback, NULL);

    /* COLUMN 2, ROW 1: Axis buttons (X, Y, Z) */
    n = 0;
    XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
    XtSetArg(args[n], XtNfromHoriz, nav_box); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    axis_box = XtCreateManagedWidget("axisBox", boxWidgetClass, form, args, n);

    /* Axis buttons */
    const char *axis_labels[] = {"X", "Y", "Z"};
    for (i = 0; i < 3; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, axis_labels[i]); n++;
        button = XtCreateManagedWidget(axis_labels[i], commandWidgetClass, axis_box, args, n);
        XtAddCallback(button, XtNcallback, axis_button_callback, (XtPointer)(long)i);
    }

    /* Terrain-following physical height option for X/Y cross-sections. */
    n = 0;
    XtSetArg(args[n], XtNlabel, "z_phys: N/A"); n++;
    z_phys_button_widget = XtCreateManagedWidget("z_phys", commandWidgetClass, axis_box, args, n);
    XtAddCallback(z_phys_button_widget, XtNcallback, z_phys_button_callback, NULL);
    update_z_phys_button(pf);
    
    /* Map button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Map: OFF"); n++;
    map_button_widget = XtCreateManagedWidget("Map", commandWidgetClass, axis_box, args, n);
    XtAddCallback(map_button_widget, XtNcallback, map_button_callback, NULL);

    /* COLUMN 2, ROW 2: Tools box (Colorbar, Distrib, Quiver, Zoom) */
    Widget tools_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert, axis_box); n++;
    XtSetArg(args[n], XtNfromHoriz, nav_box); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    tools_box = XtCreateManagedWidget("toolsBox", boxWidgetClass, form, args, n);

    /* Colorbar button - combined dialog for colormap, range, and scale */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Colorbar"); n++;
    button = XtCreateManagedWidget("colorbar", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, colorbar_button_callback, NULL);

    /* Distribution button for current layer histogram */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Distrib"); n++;
    button = XtCreateManagedWidget("distribution", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, distribution_button_callback, NULL);

    /* Quiver button for vector field display */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Quiver"); n++;
    button = XtCreateManagedWidget("quiver", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, quiver_button_callback, NULL);

    /* FFT button for energy spectrum analysis */
    n = 0;
    XtSetArg(args[n], XtNlabel, "FFT"); n++;
    button = XtCreateManagedWidget("fft2d", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, fft2d_button_callback, NULL);

    /* Zoom buttons */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Z+"); n++;
    button = XtCreateManagedWidget("zoomIn", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, zoom_in_callback, NULL);

    n = 0;
    XtSetArg(args[n], XtNlabel, "Z-"); n++;
    button = XtCreateManagedWidget("zoomOut", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, zoom_out_callback, NULL);

    /* COLUMN 3: Level buttons (show if any timestep has multiple levels) */
    /* Use max_levels_all_timesteps for multi-timestep mode, pf->n_levels for single plotfile */
    int total_levels = (n_timesteps > 1) ? max_levels_all_timesteps : pf->n_levels;
    if (total_levels > 1) {
        n = 0;
        XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
        XtSetArg(args[n], XtNfromHoriz, axis_box); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
        XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
        XtSetArg(args[n], XtNleft, XawChainLeft); n++;
        Widget level_box = XtCreateManagedWidget("levelBox", boxWidgetClass, form, args, n);

        /* Add level buttons for all levels across all timesteps (limit to 10) */
        int max_levels = total_levels < 10 ? total_levels : 10;
        for (i = 0; i < max_levels; i++) {
            n = 0;
            snprintf(label_text, sizeof(label_text), "Level %d", i);
            XtSetArg(args[n], XtNlabel, label_text); n++;
            button = XtCreateManagedWidget(label_text, commandWidgetClass, level_box, args, n);
            XtAddCallback(button, XtNcallback, level_button_callback, (XtPointer)(long)i);
        }

        /* Add overlay toggle button if more than one level possible */
        if (total_levels > 1) {
            n = 0;
            XtSetArg(args[n], XtNlabel, "Overlay: OFF"); n++;
            overlay_button = XtCreateManagedWidget("overlay", commandWidgetClass, level_box, args, n);
            XtAddCallback(overlay_button, XtNcallback, overlay_button_callback, NULL);
        }
    }

    /* ROW 2: Time navigation (only if multiple timesteps) */
    if (n_timesteps > 1) {
        Widget time_box;
        n = 0;
        XtSetArg(args[n], XtNfromVert, nav_box); n++;
        XtSetArg(args[n], XtNfromHoriz, var_viewport); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
        XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
        XtSetArg(args[n], XtNleft, XawChainLeft); n++;
        time_box = XtCreateManagedWidget("timeBox", boxWidgetClass, form, args, n);

        /* Time label */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Time"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        XtCreateManagedWidget("timeText", labelWidgetClass, time_box, args, n);

        /* Time navigation buttons (</>)  */
        const char *time_labels[] = {"<", ">"};
        for (i = 0; i < 2; i++) {
            n = 0;
            XtSetArg(args[n], XtNlabel, time_labels[i]); n++;
            button = XtCreateManagedWidget(time_labels[i], commandWidgetClass, time_box, args, n);
            XtAddCallback(button, XtNcallback, time_nav_button_callback, (XtPointer)(long)i);
        }

        /* Time index display label */
        n = 0;
        XtSetArg(args[n], XtNlabel, "1/1"); n++;
        XtSetArg(args[n], XtNwidth, 60); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        time_label = XtCreateManagedWidget("timeLabel", labelWidgetClass, time_box, args, n);

        /* Time Jump button */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Jump"); n++;
        button = XtCreateManagedWidget("timeJump", commandWidgetClass, time_box, args, n);
        XtAddCallback(button, XtNcallback, time_jump_button_callback, NULL);

        /* Time Series button */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Series"); n++;
        button = XtCreateManagedWidget("timeSeries", commandWidgetClass, time_box, args, n);
        XtAddCallback(button, XtNcallback, time_series_button_callback, NULL);

        n = 0;
        XtSetArg(args[n], XtNlabel, "Video"); n++;
        button = XtCreateManagedWidget("video", commandWidgetClass, time_box, args, n);
        XtAddCallback(button, XtNcallback, video_button_callback, NULL);
    }

    /* Colorbar widget */
    n = 0;
    XtSetArg(args[n], XtNfromVert, info_label); n++;
    XtSetArg(args[n], XtNfromHoriz, canvas_widget); n++;
    XtSetArg(args[n], XtNwidth, 100); n++;
    XtSetArg(args[n], XtNheight, canvas_height); n++;
    XtSetArg(args[n], XtNborderWidth, 2); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    colorbar_widget = XtCreateManagedWidget("colorbar", simpleWidgetClass, form, args, n);
    
    XtRealizeWidget(toplevel);
    
    /* Get canvas window and create GC */
    canvas = XtWindow(canvas_widget);
    gc = XCreateGC(display, canvas, 0, NULL);
    XSetForeground(display, gc, BlackPixel(display, screen));
    XSetFillStyle(display, gc, FillSolid);
    XSetFunction(display, gc, GXcopy);
    
    /* Get colorbar window */
    colorbar = XtWindow(colorbar_widget);
    colorbar_gc = XCreateGC(display, colorbar, 0, NULL);
    XSetFillStyle(display, colorbar_gc, FillSolid);
    XSetFunction(display, colorbar_gc, GXcopy);
    
    /* Create text GC for overlay */
    text_gc = XCreateGC(display, canvas, 0, NULL);
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    XSetBackground(display, text_gc, WhitePixel(display, screen));
    if (font) XSetFont(display, text_gc, font->fid);
    
    /* Allocate pixel buffer */
    pixel_data_size = (size_t)canvas_width * (size_t)canvas_height;
    pixel_data = (unsigned long *)malloc(pixel_data_size * sizeof(unsigned long));
    pixmap = XCreatePixmap(display, canvas, canvas_width, canvas_height, 
                          DefaultDepth(display, screen));
    colorbar_pixmap = XCreatePixmap(display, colorbar, 100, 256,
                                   DefaultDepth(display, screen));
    
    /* Add event handlers */
    XSelectInput(display, canvas, ExposureMask | KeyPressMask | PointerMotionMask | ButtonPressMask);
    XSelectInput(display, colorbar, ExposureMask);
    
    /* Add mouse event handlers - use raw event handler for proper event handling */
    XtAddRawEventHandler(canvas_widget, PointerMotionMask, False, canvas_motion_handler, NULL);
    XtAddRawEventHandler(canvas_widget, ButtonPressMask, False, canvas_button_handler, NULL);
    XtAddRawEventHandler(canvas_widget, ButtonReleaseMask, False, canvas_button_release_handler, NULL);
}

static void get_plotfile_basename(const char *path, char *name, size_t name_size) {
    if (!name || name_size == 0) return;
    name[0] = '\0';
    if (!path || !path[0]) {
        snprintf(name, name_size, "unknown");
        return;
    }

    size_t end = strlen(path);
    while (end > 1 && path[end - 1] == '/') end--;
    size_t start = end;
    while (start > 0 && path[start - 1] != '/') start--;
    size_t length = end - start;
    if (length >= name_size) length = name_size - 1;
    memcpy(name, path + start, length);
    name[length] = '\0';
}

/* Update info label */
void update_info_label(PlotfileData *pf) {
    char text[512];
    char plotfile_name[128];
    const char *axis_names[] = {"X", "Y", "Z"};
    int max_idx = pf->grid_dims[pf->slice_axis];
    get_plotfile_basename(pf->plotfile_dir, plotfile_name, sizeof(plotfile_name));
    
    if (hover_value_text[0] != '\0') {
        if (pf->n_levels > 1) {
            snprintf(text, sizeof(text), 
                     "%s | File: %s | Level: %d | Axis: %s | Layer: %d/%d | Time: %.3f | %s",
                     pf->variables[pf->current_var],
                     plotfile_name,
                     pf->current_level,
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time,
                     hover_value_text);
        } else {
            snprintf(text, sizeof(text), 
                     "%s | File: %s | Axis: %s | Layer: %d/%d | Time: %.3f | %s",
                     pf->variables[pf->current_var],
                     plotfile_name,
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time,
                     hover_value_text);
        }
    } else {
        if (pf->n_levels > 1) {
            snprintf(text, sizeof(text), 
                     "%s | File: %s | Level: %d | Axis: %s | Layer: %d/%d | Time: %.3f",
                     pf->variables[pf->current_var],
                     plotfile_name,
                     pf->current_level,
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time);
        } else {
            snprintf(text, sizeof(text), 
                     "%s | File: %s | Axis: %s | Layer: %d/%d | Time: %.3f",
                     pf->variables[pf->current_var],
                     plotfile_name,
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time);
        }
    }
    
    Arg args[1];
    XtSetArg(args[0], XtNlabel, text);
    XtSetValues(info_label, args, 1);
}

/* Variable button callback */
void var_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int var = (int)(long)client_data;
    if (global_pf && var < global_pf->n_vars) {
        global_pf->current_var = var;
        read_variable_data(global_pf, var);

        /* If overlay mode is on, reload all overlay levels with the new variable */
        if (global_pf->overlay_mode) {
            load_all_levels(global_pf, var);
        }

        zoom_reset();
        update_info_label(global_pf);
        render_slice(global_pf);
        update_distribution_histogram(-1);  /* Auto-update distribution popup */
    }
}

void update_z_phys_button(PlotfileData *pf) {
    if (!z_phys_button_widget || !pf) return;

    int available = find_variable_index(pf, "z_phys") >= 0;
    if (!available) {
        XtVaSetValues(z_phys_button_widget, XtNlabel, "z_phys: unavailable", NULL);
        XtSetSensitive(z_phys_button_widget, False);
    } else if (pf->slice_axis == 2) {
        XtVaSetValues(z_phys_button_widget, XtNlabel, "z_phys: N/A", NULL);
        XtSetSensitive(z_phys_button_widget, False);
    } else {
        XtVaSetValues(z_phys_button_widget, XtNlabel,
                      pf->use_z_phys ? "z_phys: ON" : "z_phys: OFF", NULL);
        XtSetSensitive(z_phys_button_widget, True);
    }
}

void z_phys_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w;
    (void)client_data;
    (void)call_data;
    if (!global_pf || global_pf->slice_axis == 2 ||
        find_variable_index(global_pf, "z_phys") < 0) return;

    global_pf->use_z_phys = !global_pf->use_z_phys;
    zoom_reset();
    update_z_phys_button(global_pf);
    render_slice(global_pf);
}

/* Axis button callback */
void axis_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int axis = (int)(long)client_data;
    if (global_pf) {
        global_pf->slice_axis = axis;
        global_pf->slice_idx = 0;  /* Start at first layer */
        zoom_reset();

        /* Update quiver components to match new axis */
        if (quiver_data.enabled) {
            char default_x[64], default_y[64];
            get_default_quiver_components(global_pf, default_x, default_y);
            quiver_data.x_comp_index = find_variable_index(global_pf, default_x);
            quiver_data.y_comp_index = find_variable_index(global_pf, default_y);
            /* Update dialog labels if dialog is open */
            if (quiver_data.x_comp_text && quiver_data.x_comp_index >= 0)
                XtVaSetValues(quiver_data.x_comp_text, XtNlabel, global_pf->variables[quiver_data.x_comp_index], NULL);
            if (quiver_data.y_comp_text && quiver_data.y_comp_index >= 0)
                XtVaSetValues(quiver_data.y_comp_text, XtNlabel, global_pf->variables[quiver_data.y_comp_index], NULL);
        }

        update_layer_label(global_pf);
        update_z_phys_button(global_pf);
        update_info_label(global_pf);
        render_slice(global_pf);
        update_distribution_histogram(-1);  /* Auto-update distribution popup */
    }
}

/* Map button callback */
void map_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        /* Check if lon_m and lat_m variables exist */
        int lon_idx = find_variable_index(global_pf, "lon_m");
        int lat_idx = find_variable_index(global_pf, "lat_m");
        
        if (lon_idx >= 0 && lat_idx >= 0) {
            if (!map_layers_available) {
                show_map_unavailable_dialog("Map layers not installed. To enable Map function, include `PLTVIEW_WITH_MAP=1` before pip install or upgrade.");
                return;
            }
            /* Toggle map mode */
            global_pf->map_mode = !global_pf->map_mode;

            /* Update button label */
            if (global_pf->map_mode) {
                XtVaSetValues(map_button_widget, XtNlabel, "Map: ON", NULL);
                printf("Map mode enabled: using lon_m and lat_m for coordinates\n");
            } else {
                XtVaSetValues(map_button_widget, XtNlabel, "Map: OFF", NULL);
                printf("Map mode disabled: using physical coordinates\n");
                if (map_dialog_shell) XtPopdown(map_dialog_shell);
            }
            
            update_info_label(global_pf);
            render_slice(global_pf);
            if (global_pf->map_mode) {
                show_map_settings_dialog(global_pf);
            }
        } else {
            printf("Map mode requires lon_m and lat_m variables (not found)\n");
            show_map_unavailable_dialog("lat_m and lon_m are not available");
        }
    }
}

static unsigned long get_named_color_pixel(const char *name, unsigned long fallback) {
    XColor color, exact;
    Colormap cmap = DefaultColormap(display, screen);
    if (XAllocNamedColor(display, cmap, name, &color, &exact)) {
        return color.pixel;
    }
    return fallback;
}

static void update_map_color_pixel(void) {
    switch (map_color_option) {
        case 1:
            map_color_pixel = get_named_color_pixel("red", BlackPixel(display, screen));
            break;
        case 2:
            map_color_pixel = get_named_color_pixel("gray", BlackPixel(display, screen));
            break;
        case 3:
            map_color_pixel = WhitePixel(display, screen);
            break;
        case 0:
        default:
            map_color_pixel = BlackPixel(display, screen);
            break;
    }
}

static int dir_exists(const char *path) {
    struct stat st;
    return (path && path[0] && stat(path, &st) == 0 && S_ISDIR(st.st_mode));
}

static void init_map_layers_dir(const char *argv0) {
    const char *env_dir = getenv("PLTVIEW_MAP_LAYERS");
    if (dir_exists(env_dir)) {
        strncpy(map_layers_dir, env_dir, sizeof(map_layers_dir) - 1);
        map_layers_dir[sizeof(map_layers_dir) - 1] = '\0';
    } else if (dir_exists("map_layers")) {
        strncpy(map_layers_dir, "map_layers", sizeof(map_layers_dir) - 1);
        map_layers_dir[sizeof(map_layers_dir) - 1] = '\0';
    } else if (argv0 && argv0[0]) {
        char resolved[MAX_PATH];
        if (realpath(argv0, resolved)) {
            char *slash = strrchr(resolved, '/');
            if (slash) {
                *slash = '\0';
                char candidate[MAX_PATH];
                snprintf(candidate, sizeof(candidate), "%s/map_layers", resolved);
                if (dir_exists(candidate)) {
                    strncpy(map_layers_dir, candidate, sizeof(map_layers_dir) - 1);
                    map_layers_dir[sizeof(map_layers_dir) - 1] = '\0';
                }
            }
        }
    }

    map_layers_available = dir_exists(map_layers_dir);
}

static void scan_coastline_files(void) {
    if (n_coastlines > 0) return;
    if (!map_layers_available) return;

    DIR *dir = opendir(map_layers_dir);
    if (!dir) return;

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] == '.') continue;

        const char *name = entry->d_name;
        size_t len = strlen(name);
        if (len < 5) continue;
        int is_json = (len >= 5 && strcmp(name + len - 5, ".json") == 0);
        int is_geojson = (len >= 8 && strcmp(name + len - 8, ".geojson") == 0);
        if (!is_json && !is_geojson) continue;

        if (n_coastlines >= MAX_COASTLINES) break;

        CoastlineEntry *ce = &coastlines[n_coastlines++];
        snprintf(ce->filename, sizeof(ce->filename), "%s/%s", map_layers_dir, name);

        /* Build label without extension */
        strncpy(ce->label, name, sizeof(ce->label) - 1);
        ce->label[sizeof(ce->label) - 1] = '\0';
        char *dot = strrchr(ce->label, '.');
        if (dot) *dot = '\0';

        ce->enabled = 0;
        ce->bbox_loaded = 0;
        ce->lon_min = ce->lat_min = 1e30;
        ce->lon_max = ce->lat_max = -1e30;
        ce->button = NULL;
    }

    closedir(dir);
}

static int compute_geojson_bbox(const char *path, double *lon_min, double *lon_max,
                                double *lat_min, double *lat_max) {
    FILE *fp = fopen(path, "r");
    if (!fp) return 0;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    if (fsize <= 0) {
        fclose(fp);
        return 0;
    }
    fseek(fp, 0, SEEK_SET);

    char *buf = (char *)malloc((size_t)fsize + 1);
    if (!buf) {
        fclose(fp);
        return 0;
    }
    size_t nread = fread(buf, 1, (size_t)fsize, fp);
    buf[nread] = '\0';
    fclose(fp);

    int in_coords = 0;
    int coords_pending = 0;
    int depth = 0;
    int coords_depth = -1;

    double point_vals[2];
    int nums_in_point = 0;

    *lon_min = 1e30;
    *lon_max = -1e30;
    *lat_min = 1e30;
    *lat_max = -1e30;

    for (char *p = buf; *p; p++) {
        if (!in_coords) {
            if (*p == 'c' && strncmp(p, "coordinates", 11) == 0) {
                coords_pending = 1;
                p += 10;
                continue;
            }
        }

        if (*p == '[') {
            depth++;
            if (coords_pending && !in_coords) {
                in_coords = 1;
                coords_pending = 0;
                coords_depth = depth;
                nums_in_point = 0;
            }
            continue;
        }

        if (*p == ']') {
            depth--;
            if (in_coords && coords_depth >= 0 && depth < coords_depth) {
                in_coords = 0;
                coords_depth = -1;
                nums_in_point = 0;
            }
            continue;
        }

        if (in_coords && (*p == '-' || (*p >= '0' && *p <= '9'))) {
            char *endptr = NULL;
            double val = strtod(p, &endptr);
            if (endptr && endptr != p) {
                point_vals[nums_in_point++] = val;
                if (nums_in_point == 2) {
                    double lon = point_vals[0];
                    double lat = point_vals[1];

                    if (lon < *lon_min) *lon_min = lon;
                    if (lon > *lon_max) *lon_max = lon;
                    if (lat < *lat_min) *lat_min = lat;
                    if (lat > *lat_max) *lat_max = lat;

                    nums_in_point = 0;
                }
                p = endptr - 1;
            }
        }
    }

    free(buf);
    return (*lon_min <= *lon_max && *lat_min <= *lat_max);
}

static void update_coastline_button_label(CoastlineEntry *ce) {
    if (!ce->button) return;
    char label[160];
    snprintf(label, sizeof(label), "%s: %s", ce->label, ce->enabled ? "ON" : "OFF");
    XtVaSetValues(ce->button, XtNlabel, label, NULL);
}

static void auto_detect_coastlines(void) {
    if (!map_has_bounds) return;
    for (int i = 0; i < n_coastlines; i++) {
        CoastlineEntry *ce = &coastlines[i];
        if (!ce->bbox_loaded) {
            if (compute_geojson_bbox(ce->filename, &ce->lon_min, &ce->lon_max, &ce->lat_min, &ce->lat_max)) {
                ce->bbox_loaded = 1;
            }
        }
        if (ce->bbox_loaded) {
            int overlap = !(ce->lon_max < map_last_lon_min || ce->lon_min > map_last_lon_max ||
                            ce->lat_max < map_last_lat_min || ce->lat_min > map_last_lat_max);
            ce->enabled = overlap ? 1 : 0;
        }
    }
    map_auto_detected = 1;
}

void map_color_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    map_color_option = (int)(long)client_data;
    update_map_color_pixel();
    if (global_pf && global_pf->map_mode) {
        render_slice(global_pf);
    }
}

void map_remove_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w;
    (void)client_data;
    (void)call_data;
    if (!global_pf) return;
    global_pf->map_mode = 0;
    XtVaSetValues(map_button_widget, XtNlabel, "Map: OFF", NULL);
    if (map_dialog_shell) {
        XtPopdown(map_dialog_shell);
    }
    update_info_label(global_pf);
    render_slice(global_pf);
}

void show_map_settings_dialog(PlotfileData *pf) {
    (void)pf;
    if (map_dialog_shell) {
        XtPopup(map_dialog_shell, XtGrabNone);
        return;
    }

    update_map_color_pixel();
    scan_coastline_files();
    if (!map_auto_detected) {
        auto_detect_coastlines();
    }

    map_dialog_shell = XtVaCreatePopupShell(
        "mapSettings", transientShellWidgetClass, toplevel,
        XtNtitle, "Map Properties",
        NULL);

    Widget map_form = XtVaCreateManagedWidget("mapForm", formWidgetClass, map_dialog_shell, NULL);
    Widget title = XtVaCreateManagedWidget("mapTitle", labelWidgetClass, map_form,
                                           XtNlabel, "Map Properties",
                                           XtNborderWidth, 0,
                                           NULL);

    Widget color_label = XtVaCreateManagedWidget("mapColorLabel", labelWidgetClass, map_form,
                                                 XtNlabel, "Coastline Color:",
                                                 XtNfromVert, title,
                                                 XtNborderWidth, 0,
                                                 NULL);

    Widget color_box = XtVaCreateManagedWidget("mapColorBox", boxWidgetClass, map_form,
                                               XtNfromVert, color_label,
                                               XtNorientation, XtorientHorizontal,
                                               NULL);

    Widget black_btn = XtVaCreateManagedWidget("mapColorBlack", commandWidgetClass, color_box,
                                               XtNlabel, "Black", NULL);
    Widget red_btn = XtVaCreateManagedWidget("mapColorRed", commandWidgetClass, color_box,
                                             XtNlabel, "Red", NULL);
    Widget gray_btn = XtVaCreateManagedWidget("mapColorGray", commandWidgetClass, color_box,
                                              XtNlabel, "Gray", NULL);
    Widget white_btn = XtVaCreateManagedWidget("mapColorWhite", commandWidgetClass, color_box,
                                               XtNlabel, "White", NULL);

    XtAddCallback(black_btn, XtNcallback, map_color_callback, (XtPointer)0);
    XtAddCallback(red_btn, XtNcallback, map_color_callback, (XtPointer)1);
    XtAddCallback(gray_btn, XtNcallback, map_color_callback, (XtPointer)2);
    XtAddCallback(white_btn, XtNcallback, map_color_callback, (XtPointer)3);

    Widget list_label = XtVaCreateManagedWidget("mapListLabel", labelWidgetClass, map_form,
                                                XtNlabel, "Map Layers:",
                                                XtNfromVert, color_box,
                                                XtNborderWidth, 0,
                                                NULL);

    Widget list_box = XtVaCreateManagedWidget("mapListBox", boxWidgetClass, map_form,
                                              XtNfromVert, list_label,
                                              XtNorientation, XtorientVertical,
                                              NULL);

    for (int i = 0; i < n_coastlines; i++) {
        CoastlineEntry *ce = &coastlines[i];
        ce->button = XtVaCreateManagedWidget("mapCoastline", commandWidgetClass, list_box,
                                             XtNlabel, ce->label,
                                             NULL);
        XtAddCallback(ce->button, XtNcallback, map_coastline_toggle_callback, (XtPointer)(long)i);
        update_coastline_button_label(ce);
    }

    Widget action_box = XtVaCreateManagedWidget("mapActionBox", boxWidgetClass, map_form,
                                                XtNfromVert, list_box,
                                                XtNorientation, XtorientHorizontal,
                                                NULL);

    Widget remove_btn = XtVaCreateManagedWidget("mapRemove", commandWidgetClass, action_box,
                                                XtNlabel, "Remove", NULL);
    XtAddCallback(remove_btn, XtNcallback, map_remove_callback, NULL);

    XtPopup(map_dialog_shell, XtGrabNone);
}

void map_coastline_toggle_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w;
    (void)call_data;
    int idx = (int)(long)client_data;
    if (idx < 0 || idx >= n_coastlines) return;

    CoastlineEntry *ce = &coastlines[idx];
    ce->enabled = !ce->enabled;
    update_coastline_button_label(ce);

    if (ce->enabled) {
        map_coastlines_enabled = 1;
    }

    if (global_pf && global_pf->map_mode) {
        render_slice(global_pf);
    }
}

void map_unavailable_ok_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (map_unavailable_shell) {
        XtPopdown(map_unavailable_shell);
    }
}

void show_map_unavailable_dialog(const char *message) {
    const char *label_text = (message && message[0]) ? message : "Map is unavailable";

    if (map_unavailable_shell) {
        if (map_unavailable_label) {
            XtVaSetValues(map_unavailable_label, XtNlabel, label_text, NULL);
        }
        XtPopup(map_unavailable_shell, XtGrabNone);
        return;
    }

    map_unavailable_shell = XtVaCreatePopupShell(
        "mapUnavailable", transientShellWidgetClass, toplevel,
        XtNtitle, "Map",
        NULL);

    Widget msg_form = XtVaCreateManagedWidget("mapUnavailableForm", formWidgetClass,
                                              map_unavailable_shell, NULL);
    map_unavailable_label = XtVaCreateManagedWidget("mapUnavailableLabel", labelWidgetClass,
                                                    msg_form,
                                                    XtNlabel, label_text,
                                                    XtNborderWidth, 0,
                                                    NULL);
    Widget ok_btn = XtVaCreateManagedWidget("mapUnavailableOk", commandWidgetClass, msg_form,
                                            XtNlabel, "OK",
                                            XtNfromVert, map_unavailable_label,
                                            NULL);

    XtAddCallback(ok_btn, XtNcallback, map_unavailable_ok_callback, NULL);

    XtPopup(map_unavailable_shell, XtGrabNone);
}

/* Level button callback */
void level_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int level = (int)(long)client_data;
    if (!global_pf) return;

    /* Check if the requested level is available at the current timestep */
    if (level >= global_pf->n_levels) {
        show_level_warning(level);
        return;
    }

    zoom_reset();
    free_z_phys_cache(global_pf);
    global_pf->current_level = level;

    /* Reload data for new level */
    global_pf->n_boxes = 0;
    read_cell_h(global_pf);
    read_variable_data(global_pf, global_pf->current_var);

    /* Clamp slice_idx if new level has fewer layers */
    int max_idx = global_pf->grid_dims[global_pf->slice_axis] - 1;
    if (global_pf->slice_idx > max_idx) {
        global_pf->slice_idx = max_idx;
    }

    update_layer_label(global_pf);
    update_z_phys_button(global_pf);
    update_info_label(global_pf);
    render_slice(global_pf);
    update_distribution_histogram(-1);  /* Auto-update distribution popup */
}

/* Overlay toggle button callback */
void overlay_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        global_pf->overlay_mode = !global_pf->overlay_mode;

        if (global_pf->overlay_mode) {
            /* Load all levels data for overlay */
            printf("Enabling overlay mode - loading all levels...\n");
            if (global_pf->n_levels > 1) {
                load_all_levels(global_pf, global_pf->current_var);
            }

            /* Update button label */
            Arg args[1];
            if (global_pf->n_levels > 1) {
                XtSetArg(args[0], XtNlabel, "Overlay: ON");
            } else {
                XtSetArg(args[0], XtNlabel, "Overlay: ON (no L1)");
            }
            XtSetValues(overlay_button, args, 1);
        } else {
            /* Disable overlay mode */
            printf("Disabling overlay mode\n");
            free_all_levels(global_pf);

            /* Update button label */
            Arg args[1];
            XtSetArg(args[0], XtNlabel, "Overlay: OFF");
            XtSetValues(overlay_button, args, 1);
        }

        render_slice(global_pf);
    }
}

/* Navigation button callback (+/-) */
void nav_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int dir = (int)(long)client_data;  /* 0 = minus, 1 = plus */
    if (global_pf) {
        int max_idx = global_pf->grid_dims[global_pf->slice_axis] - 1;

        if (dir == 1) {
            /* Plus: go to next layer, wrap to 0 if at end */
            global_pf->slice_idx++;
            if (global_pf->slice_idx > max_idx) {
                global_pf->slice_idx = 0;
            }
        } else {
            /* Minus: go to previous layer, wrap to end if at 0 */
            global_pf->slice_idx--;
            if (global_pf->slice_idx < 0) {
                global_pf->slice_idx = max_idx;
            }
        }

        update_layer_label(global_pf);
        update_info_label(global_pf);
        render_slice(global_pf);
        update_distribution_histogram(-1);  /* Auto-update distribution popup */
    }
}

/* Jump to specific layer positions - button-based for X11 forwarding reliability */
void jump_to_layer_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    long layer_type = (long)client_data;
    
    if (global_pf) {
        int max_idx = global_pf->grid_dims[global_pf->slice_axis];
        int new_idx = global_pf->slice_idx;
        
        switch (layer_type) {
            case 0: new_idx = 0; break;                    /* First */
            case 1: new_idx = max_idx - 1; break;          /* Last */
            case 2: new_idx = max_idx / 2; break;          /* Middle */
            case 3: new_idx = max_idx / 4; break;          /* 1/4 */
            case 4: new_idx = 3 * max_idx / 4; break;      /* 3/4 */
        }
        
        if (new_idx >= 0 && new_idx < max_idx) {
            global_pf->slice_idx = new_idx;
            update_layer_label(global_pf);
            update_info_label(global_pf);
            render_slice(global_pf);
            update_distribution_histogram(-1);  /* Auto-update distribution popup */
        }
    }

    /* Close the dialog */
    Widget shell = XtParent(XtParent(w));
    XtPopdown(shell);
    XtDestroyWidget(shell);
    dialog_active = 0;
    active_text_widget = NULL;
}

/* Structure to pass both text widget and shell to callback */
typedef struct {
    Widget text_widget;
    Widget dialog_shell;
} JumpDialogData;

/* Jump to typed layer number */
void jump_to_typed_layer_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    JumpDialogData *data = (JumpDialogData *)client_data;
    
    if (global_pf && data) {
        String value;
        Arg args[1];
        XtSetArg(args[0], XtNstring, &value);
        XtGetValues(data->text_widget, args, 1);
        
        if (value && strlen(value) > 0) {
            int layer = atoi(value);
            int max_idx = global_pf->grid_dims[global_pf->slice_axis];
            
            /* Convert from 1-indexed to 0-indexed and clamp */
            layer = layer - 1;
            if (layer < 0) layer = 0;
            if (layer >= max_idx) layer = max_idx - 1;
            
            global_pf->slice_idx = layer;
            update_layer_label(global_pf);
            update_info_label(global_pf);
            render_slice(global_pf);
            update_distribution_histogram(-1);  /* Auto-update distribution popup */
        }

        /* Close the dialog */
        XtPopdown(data->dialog_shell);
        XtDestroyWidget(data->dialog_shell);
        free(data);
        dialog_active = 0;
        active_text_widget = NULL;
    }
}

/* Close jump dialog */
void jump_dialog_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Widget shell = (Widget)client_data;
    XtPopdown(shell);
    XtDestroyWidget(shell);
    dialog_active = 0;
    active_text_widget = NULL;
}

/* Jump button callback - hybrid dialog with both text input and quick-jump buttons */
void jump_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        Arg args[10];
        int n;
        Widget dialog_shell, form, label, button, text_widget, text_label;
        char msg[128];
        int max_idx = global_pf->grid_dims[global_pf->slice_axis];
        
        snprintf(msg, sizeof(msg), "Jump to layer (1-%d)", max_idx);
        
        n = 0;
        XtSetArg(args[n], XtNtitle, "Jump to Layer"); n++;
        dialog_shell = XtCreatePopupShell("jumpDialog", transientShellWidgetClass, toplevel, args, n);
        
        n = 0;
        form = XtCreateManagedWidget("form", formWidgetClass, dialog_shell, args, n);
        
        /* Title label */
        n = 0;
        XtSetArg(args[n], XtNlabel, msg); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("label", labelWidgetClass, form, args, n);
        
        /* Text input section */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNlabel, "Type layer:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        text_label = XtCreateManagedWidget("textLabel", labelWidgetClass, form, args, n);
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, text_label); n++;
        XtSetArg(args[n], XtNwidth, 100); n++;
        XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
        XtSetArg(args[n], XtNstring, ""); n++;
        text_widget = XtCreateManagedWidget("textInput", asciiTextWidgetClass, form, args, n);
        
        /* Create data structure to pass to callback */
        JumpDialogData *jump_data = malloc(sizeof(JumpDialogData));
        jump_data->text_widget = text_widget;
        jump_data->dialog_shell = dialog_shell;
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, text_label); n++;
        XtSetArg(args[n], XtNfromHoriz, text_widget); n++;
        XtSetArg(args[n], XtNlabel, "Go"); n++;
        button = XtCreateManagedWidget("goButton", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_to_typed_layer_callback, (XtPointer)jump_data);
        
        /* Or quick jump label */
        n = 0;
        XtSetArg(args[n], XtNfromVert, text_widget); n++;
        XtSetArg(args[n], XtNlabel, "Or quick jump:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("orLabel", labelWidgetClass, form, args, n);
        
        /* Quick jump buttons */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNlabel, "First (1)"); n++;
        button = XtCreateManagedWidget("first", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_to_layer_callback, (XtPointer)0);
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "1/4"); n++;
        button = XtCreateManagedWidget("quarter", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_to_layer_callback, (XtPointer)3);
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "Middle"); n++;
        button = XtCreateManagedWidget("middle", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_to_layer_callback, (XtPointer)2);
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "3/4"); n++;
        button = XtCreateManagedWidget("threequarter", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_to_layer_callback, (XtPointer)4);
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        snprintf(msg, sizeof(msg), "Last (%d)", max_idx);
        XtSetArg(args[n], XtNlabel, msg); n++;
        button = XtCreateManagedWidget("last", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_to_layer_callback, (XtPointer)1);
        
        n = 0;
        XtSetArg(args[n], XtNfromVert, button); n++;
        XtSetArg(args[n], XtNlabel, "Close"); n++;
        button = XtCreateManagedWidget("close", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, jump_dialog_close_callback, (XtPointer)dialog_shell);
        
        XtRealizeWidget(dialog_shell);
        XtPopup(dialog_shell, XtGrabExclusive);

        /* Set keyboard focus to text input - needed for remote X11 */
        XtSetKeyboardFocus(dialog_shell, text_widget);

        /* Force the text widget to accept focus */
        XSync(display, False);
        Time time = CurrentTime;
        XtCallAcceptFocus(text_widget, &time);

        dialog_active = 1;
        active_text_widget = text_widget;
    }
}

/* Structure for colorbar dialog (combines range + colormap + scale) */
typedef struct {
    Widget max_text;      /* Max input (shown first) */
    Widget min_text;      /* Min input (shown second) */
    Widget dialog_shell;
    Widget scale_buttons[3];  /* Linear, Log(+), Log(-) buttons */
} ColorbarDialogData;

/* Global pointer to colorbar dialog data for keyboard input */
ColorbarDialogData *active_colorbar_dialog = NULL;
int active_field = 0;  /* 0 = max, 1 = min */

/* Update scale button labels to show current selection */
void update_scale_button_labels(ColorbarDialogData *data) {
    if (!data) return;
    Arg args[1];
    const char *labels[] = {"Linear", "Log(+)", "Log(-)"};
    for (int i = 0; i < 3; i++) {
        char label[16];
        if (i == scale_mode) {
            snprintf(label, sizeof(label), "[%s]", labels[i]);
        } else {
            snprintf(label, sizeof(label), "%s", labels[i]);
        }
        XtSetArg(args[0], XtNlabel, label);
        XtSetValues(data->scale_buttons[i], args, 1);
    }
}

/* Scale mode toggle callback */
void colorbar_scale_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int new_mode = (int)(long)client_data;
    scale_mode = new_mode;

    /* Update button labels to show selection */
    if (active_colorbar_dialog) {
        update_scale_button_labels(active_colorbar_dialog);
    }

    /* Re-render with new scale mode */
    if (global_pf) {
        render_slice(global_pf);
        draw_colorbar(current_vmin, current_vmax, global_pf->colormap,
                      global_pf->variables[global_pf->current_var]);
    }
}

/* Apply colorbar settings callback */
void colorbar_apply_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    ColorbarDialogData *data = (ColorbarDialogData *)client_data;

    if (data) {
        String min_str, max_str;
        Arg args[1];

        XtSetArg(args[0], XtNstring, &max_str);
        XtGetValues(data->max_text, args, 1);
        XtSetArg(args[0], XtNstring, &min_str);
        XtGetValues(data->min_text, args, 1);

        if (min_str && strlen(min_str) > 0 && max_str && strlen(max_str) > 0) {
            double new_min = atof(min_str);
            double new_max = atof(max_str);

            if (new_min < new_max) {
                custom_vmin = new_min;
                custom_vmax = new_max;
                use_custom_range = 1;

                /* Re-render with new range */
                if (global_pf) {
                    render_slice(global_pf);
                }
            }
        }

        /* Close the dialog */
        XtPopdown(data->dialog_shell);
        XtDestroyWidget(data->dialog_shell);
        free(data);
        dialog_active = 0;
        active_text_widget = NULL;
        active_colorbar_dialog = NULL;
    }
}

/* Auto range callback - reset to data-driven min/max */
void colorbar_auto_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    ColorbarDialogData *data = (ColorbarDialogData *)client_data;

    use_custom_range = 0;

    /* Re-render with auto range */
    if (global_pf) {
        render_slice(global_pf);
    }

    if (data) {
        /* Close the dialog */
        XtPopdown(data->dialog_shell);
        XtDestroyWidget(data->dialog_shell);
        free(data);
        dialog_active = 0;
        active_text_widget = NULL;
        active_colorbar_dialog = NULL;
    }
}

/* Close colorbar dialog callback */
void colorbar_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    ColorbarDialogData *data = (ColorbarDialogData *)client_data;

    if (data) {
        XtPopdown(data->dialog_shell);
        XtDestroyWidget(data->dialog_shell);
        free(data);
        dialog_active = 0;
        active_text_widget = NULL;
        active_colorbar_dialog = NULL;
    }
}

/* Ensure focus follows mouse clicks in colorbar dialog text fields */
void colorbar_text_click_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    ColorbarDialogData *data = (ColorbarDialogData *)client_data;
    if (!data || event->type != ButtonPress) return;

    XtSetKeyboardFocus(data->dialog_shell, w);
    Time time = CurrentTime;
    XtCallAcceptFocus(w, &time);

    active_text_widget = w;
    active_field = (w == data->min_text) ? 1 : 0;  /* 0=max, 1=min */
}

/* Colormap selection callback for colorbar dialog */
void colorbar_cmap_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int cmap = (int)(long)client_data;
    if (global_pf) {
        global_pf->colormap = cmap;
        render_slice(global_pf);
        draw_colorbar(current_vmin, current_vmax, cmap,
                      global_pf->variables[global_pf->current_var]);
    }
}

/* Colorbar button callback - combined dialog for range, scale, and colormap */
void colorbar_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        Arg args[10];
        int n, i;
        Widget dialog_shell, form, label, button, min_text, max_text;
        char min_str[64], max_str[64];

        /* Format current values */
        snprintf(max_str, sizeof(max_str), "%.6e", use_custom_range ? custom_vmax : current_vmax);
        snprintf(min_str, sizeof(min_str), "%.6e", use_custom_range ? custom_vmin : current_vmin);

        n = 0;
        XtSetArg(args[n], XtNtitle, "Colorbar Settings"); n++;
        dialog_shell = XtCreatePopupShell("colorbarDialog", transientShellWidgetClass, toplevel, args, n);

        n = 0;
        form = XtCreateManagedWidget("form", formWidgetClass, dialog_shell, args, n);

        /* --- Range Section --- */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Range:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("rangeLabel", labelWidgetClass, form, args, n);

        /* Max label (shown first per user request) */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNlabel, "Max:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        Widget max_label = XtCreateManagedWidget("maxLabel", labelWidgetClass, form, args, n);

        /* Max text input */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNfromHoriz, max_label); n++;
        XtSetArg(args[n], XtNwidth, 150); n++;
        XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
        XtSetArg(args[n], XtNstring, max_str); n++;
        max_text = XtCreateManagedWidget("maxInput", asciiTextWidgetClass, form, args, n);

        /* Min label */
        n = 0;
        XtSetArg(args[n], XtNfromVert, max_label); n++;
        XtSetArg(args[n], XtNlabel, "Min:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        Widget min_label = XtCreateManagedWidget("minLabel", labelWidgetClass, form, args, n);

        /* Min text input */
        n = 0;
        XtSetArg(args[n], XtNfromVert, max_label); n++;
        XtSetArg(args[n], XtNfromHoriz, min_label); n++;
        XtSetArg(args[n], XtNwidth, 150); n++;
        XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
        XtSetArg(args[n], XtNstring, min_str); n++;
        min_text = XtCreateManagedWidget("minInput", asciiTextWidgetClass, form, args, n);

        /* --- Scale Section --- */
        n = 0;
        XtSetArg(args[n], XtNfromVert, min_label); n++;
        XtSetArg(args[n], XtNlabel, "Scale:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        Widget scale_label = XtCreateManagedWidget("scaleLabel", labelWidgetClass, form, args, n);

        /* Create data structure early so scale buttons can use it */
        ColorbarDialogData *colorbar_data = malloc(sizeof(ColorbarDialogData));
        colorbar_data->max_text = max_text;
        colorbar_data->min_text = min_text;
        colorbar_data->dialog_shell = dialog_shell;

        /* Scale toggle buttons: Linear, Log(+), Log(-) */
        const char *scale_labels[] = {"Linear", "Log(+)", "Log(-)"};
        Widget prev_scale = scale_label;
        for (i = 0; i < 3; i++) {
            char btn_label[16];
            if (i == scale_mode) {
                snprintf(btn_label, sizeof(btn_label), "[%s]", scale_labels[i]);
            } else {
                snprintf(btn_label, sizeof(btn_label), "%s", scale_labels[i]);
            }
            n = 0;
            XtSetArg(args[n], XtNfromVert, min_label); n++;
            XtSetArg(args[n], XtNfromHoriz, prev_scale); n++;
            XtSetArg(args[n], XtNlabel, btn_label); n++;
            colorbar_data->scale_buttons[i] = XtCreateManagedWidget(scale_labels[i], commandWidgetClass, form, args, n);
            XtAddCallback(colorbar_data->scale_buttons[i], XtNcallback, colorbar_scale_callback, (XtPointer)(long)i);
            prev_scale = colorbar_data->scale_buttons[i];
        }

        /* --- Colormap Section --- */
        n = 0;
        XtSetArg(args[n], XtNfromVert, scale_label); n++;
        XtSetArg(args[n], XtNlabel, "Colormap:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        Widget cmap_label = XtCreateManagedWidget("cmapLabel", labelWidgetClass, form, args, n);

        /* Colormap buttons in a row */
        const char *cmap_names[] = {"viridis", "jet", "turbo", "plasma", "hot", "cool", "gray", "magma"};
        Widget prev_cmap = cmap_label;
        for (i = 0; i < 8; i++) {
            n = 0;
            XtSetArg(args[n], XtNfromVert, scale_label); n++;
            XtSetArg(args[n], XtNfromHoriz, prev_cmap); n++;
            XtSetArg(args[n], XtNlabel, cmap_names[i]); n++;
            button = XtCreateManagedWidget(cmap_names[i], commandWidgetClass, form, args, n);
            XtAddCallback(button, XtNcallback, colorbar_cmap_callback, (XtPointer)(long)i);
            prev_cmap = button;
        }

        /* Click-to-focus handlers for text fields */
        XtAddEventHandler(max_text, ButtonPressMask, False, colorbar_text_click_handler, (XtPointer)colorbar_data);
        XtAddEventHandler(min_text, ButtonPressMask, False, colorbar_text_click_handler, (XtPointer)colorbar_data);

        /* --- Action Buttons --- */
        n = 0;
        XtSetArg(args[n], XtNfromVert, cmap_label); n++;
        XtSetArg(args[n], XtNlabel, "Apply"); n++;
        button = XtCreateManagedWidget("apply", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, colorbar_apply_callback, (XtPointer)colorbar_data);

        n = 0;
        XtSetArg(args[n], XtNfromVert, cmap_label); n++;
        XtSetArg(args[n], XtNfromHoriz, button); n++;
        XtSetArg(args[n], XtNlabel, "Auto"); n++;
        Widget auto_button = XtCreateManagedWidget("auto", commandWidgetClass, form, args, n);
        XtAddCallback(auto_button, XtNcallback, colorbar_auto_callback, (XtPointer)colorbar_data);

        n = 0;
        XtSetArg(args[n], XtNfromVert, cmap_label); n++;
        XtSetArg(args[n], XtNfromHoriz, auto_button); n++;
        XtSetArg(args[n], XtNlabel, "Close"); n++;
        button = XtCreateManagedWidget("close", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, colorbar_close_callback, (XtPointer)colorbar_data);

        XtRealizeWidget(dialog_shell);
        XtPopup(dialog_shell, XtGrabExclusive);

        /* Set keyboard focus to max text input (shown first) */
        XtSetKeyboardFocus(dialog_shell, max_text);
        XSync(display, False);
        Time time = CurrentTime;
        XtCallAcceptFocus(max_text, &time);

        dialog_active = 1;
        active_text_widget = max_text;
        active_colorbar_dialog = colorbar_data;
        active_field = 0;  /* 0=max */
    }
}

/* Update layer display label */
void update_layer_label(PlotfileData *pf) {
    char text[32];
    int max_idx = pf->grid_dims[pf->slice_axis];
    snprintf(text, sizeof(text), "%d/%d", pf->slice_idx + 1, max_idx);
    
    Arg args[1];
    XtSetArg(args[0], XtNlabel, text);
    XtSetValues(layer_label, args, 1);
}

/* Colormap button callback (used in popup) */
void cmap_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int cmap = (int)(long)client_data;
    if (global_pf) {
        global_pf->colormap = cmap;
        render_slice(global_pf);
        draw_colorbar(current_vmin, current_vmax, cmap,
                      global_pf->variables[global_pf->current_var]);
    }
}

/* Close colormap dialog and apply selection */
void cmap_select_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int cmap = (int)(long)client_data;
    if (global_pf) {
        global_pf->colormap = cmap;
        render_slice(global_pf);
        draw_colorbar(current_vmin, current_vmax, cmap,
                      global_pf->variables[global_pf->current_var]);
    }

    /* Close the dialog */
    Widget shell = XtParent(XtParent(w));
    XtPopdown(shell);
    XtDestroyWidget(shell);
}

/* Close colormap dialog without changing */
void cmap_dialog_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Widget shell = (Widget)client_data;
    XtPopdown(shell);
    XtDestroyWidget(shell);
}

/* Colormap button callback - opens popup with colormap options */
void colormap_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        Arg args[10];
        int n, i;
        Widget dialog_shell, form, label, button;

        n = 0;
        XtSetArg(args[n], XtNtitle, "Select Colormap"); n++;
        dialog_shell = XtCreatePopupShell("cmapDialog", transientShellWidgetClass, toplevel, args, n);

        n = 0;
        form = XtCreateManagedWidget("form", formWidgetClass, dialog_shell, args, n);

        /* Title label */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Choose colormap:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("label", labelWidgetClass, form, args, n);

        /* Colormap buttons with numbered labels */
        const char *cmap_names[] = {"viridis", "jet", "turbo", "plasma", "hot", "cool", "gray", "magma"};
        Widget prev_button = label;
        char cmap_label[32];

        for (i = 0; i < 8; i++) {
            n = 0;
            XtSetArg(args[n], XtNfromVert, prev_button); n++;
            snprintf(cmap_label, sizeof(cmap_label), "%d. %s", i + 1, cmap_names[i]);
            XtSetArg(args[n], XtNlabel, cmap_label); n++;
            XtSetArg(args[n], XtNwidth, 100); n++;
            button = XtCreateManagedWidget(cmap_names[i], commandWidgetClass, form, args, n);
            XtAddCallback(button, XtNcallback, cmap_select_callback, (XtPointer)(long)i);
            prev_button = button;
        }

        /* Close button */
        n = 0;
        XtSetArg(args[n], XtNfromVert, prev_button); n++;
        XtSetArg(args[n], XtNlabel, "Close"); n++;
        button = XtCreateManagedWidget("close", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, cmap_dialog_close_callback, (XtPointer)dialog_shell);

        XtRealizeWidget(dialog_shell);
        XtPopup(dialog_shell, XtGrabNone);
    }
}

/* Canvas expose callback */
void canvas_expose_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data) {
        render_slice(global_pf);
    }
}

/* Colorbar expose callback */
void colorbar_expose_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        draw_colorbar(current_vmin, current_vmax, global_pf->colormap,
                      global_pf->variables[global_pf->current_var]);
    }
}

static void draw_y_label_ccw(const char *label, int x, int y_center) {
    if (!label || !label[0] || !font) return;

    int len = (int)strlen(label);
    int text_w = XTextWidth(font, label, len);
    int text_h = font->ascent + font->descent;
    if (text_w <= 0 || text_h <= 0) return;

    Pixmap tmp = XCreatePixmap(display, canvas, text_w, text_h,
                               DefaultDepth(display, screen));
    if (!tmp) return;

    GC tmp_gc = XCreateGC(display, tmp, 0, NULL);
    if (!tmp_gc) {
        XFreePixmap(display, tmp);
        return;
    }

    XSetForeground(display, tmp_gc, WhitePixel(display, screen));
    XFillRectangle(display, tmp, tmp_gc, 0, 0, text_w, text_h);
    XSetForeground(display, tmp_gc, BlackPixel(display, screen));
    if (font) {
        XSetFont(display, tmp_gc, font->fid);
    }
    XDrawString(display, tmp, tmp_gc, 0, font->ascent, label, len);

    XImage *src = XGetImage(display, tmp, 0, 0, text_w, text_h, AllPlanes, ZPixmap);
    if (!src) {
        XFreeGC(display, tmp_gc);
        XFreePixmap(display, tmp);
        return;
    }

    int dest_w = text_h;
    int dest_h = text_w;
    XImage *rot = XCreateImage(display, DefaultVisual(display, screen),
                               DefaultDepth(display, screen), ZPixmap, 0, NULL,
                               dest_w, dest_h, 32, 0);
    if (!rot) {
        XDestroyImage(src);
        XFreeGC(display, tmp_gc);
        XFreePixmap(display, tmp);
        return;
    }
    rot->data = (char *)malloc((size_t)rot->bytes_per_line * (size_t)dest_h);
    if (!rot->data) {
        XDestroyImage(rot);
        XDestroyImage(src);
        XFreeGC(display, tmp_gc);
        XFreePixmap(display, tmp);
        return;
    }

    for (int y = 0; y < text_h; y++) {
        for (int xx = 0; xx < text_w; xx++) {
            unsigned long pixel = XGetPixel(src, xx, y);
            int rx = y;
            int ry = dest_h - 1 - xx;
            XPutPixel(rot, rx, ry, pixel);
        }
    }

    int draw_y = y_center - dest_h / 2;
    XPutImage(display, canvas, text_gc, rot, 0, 0, x, draw_y, dest_w, dest_h);

    XDestroyImage(rot);
    XDestroyImage(src);
    XFreeGC(display, tmp_gc);
    XFreePixmap(display, tmp);
}
void render_slice(PlotfileData *pf) {
    int width, height;
    double *slice;
    double *z_phys_slice = NULL;
    int use_z_phys_coords = 0;
    double vmin = 1e30, vmax = -1e30;
    double vsum = 0.0;
    int vcount = 0;
    int i, j;
    char stats_text[128];

    /* A map render rebuilds this lookup from projected cell positions below. */
    map_hover_lookup_active = 0;

    /* Axis margin sizes */
    int left_margin = 60;    /* Space for Y-axis labels */
    int bottom_margin = 40;  /* Space for X-axis labels */
    int top_margin = 10;     /* Small top margin */
    int right_margin = 10;   /* Small right margin */

    /* Determine slice dimensions and physical coordinates */
    int x_axis, y_axis;  /* Which physical dimensions map to screen x,y */
    if (pf->slice_axis == 2) {       /* Z-slice: X horizontal, Y vertical */
        width = pf->grid_dims[0];
        height = pf->grid_dims[1];
        x_axis = 0;  /* X */
        y_axis = 1;  /* Y */
    } else if (pf->slice_axis == 1) { /* Y-slice: X horizontal, Z vertical */
        width = pf->grid_dims[0];
        height = pf->grid_dims[2];
        x_axis = 0;  /* X */
        y_axis = 2;  /* Z */
    } else {                          /* X-slice: Y horizontal, Z vertical */
        width = pf->grid_dims[1];
        height = pf->grid_dims[2];
        x_axis = 1;  /* Y */
        y_axis = 2;  /* Z */
    }

    slice = (double *)malloc(width * height * sizeof(double));
    {
        size_t needed_pixels = (size_t)width * (size_t)height;
        if (needed_pixels > pixel_data_size) {
            unsigned long *new_pixels = (unsigned long *)realloc(pixel_data, needed_pixels * sizeof(unsigned long));
            if (new_pixels) {
                pixel_data = new_pixels;
                pixel_data_size = needed_pixels;
            }
        }
    }
    extract_slice(pf, slice, pf->slice_axis, pf->slice_idx);

    /* z_phys is a 3D cell-centered coordinate.  Extract the same X/Z or Y/Z
     * cross-section as the scalar field so every displayed cell gets its own
     * terrain-following physical height. */
    if (pf->use_z_phys && pf->slice_axis != 2 &&
        find_variable_index(pf, "z_phys") >= 0 && ensure_z_phys_data(pf) == 0) {
        z_phys_slice = (double *)malloc((size_t)width * height * sizeof(double));
        if (z_phys_slice) {
            extract_slice_from_data(pf->z_phys_data, pf, z_phys_slice,
                                    pf->slice_axis, pf->slice_idx);
            use_z_phys_coords = 1;
        }
    }

    /* Physical coordinate ranges for axes */
    double phys_xmin = pf->prob_lo[x_axis];
    double phys_xmax = pf->prob_hi[x_axis];
    double phys_ymin = pf->prob_lo[y_axis];
    double phys_ymax = pf->prob_hi[y_axis];
    
    if (pf->map_mode) {
        /* Map mode: phys_xmin/max and phys_ymin/max will be set in the map rendering section */
    } else {
        /* Normal mode: use physical coordinates */
        phys_xmin = pf->prob_lo[x_axis];
        phys_xmax = pf->prob_hi[x_axis];
        if (use_z_phys_coords) {
            phys_ymin = INFINITY;
            phys_ymax = -INFINITY;
            for (int zb = 0; zb <= height; zb++) {
                for (int xb = 0; xb <= width; xb++) {
                    double z = z_phys_corner(z_phys_slice, width, height, xb, zb);
                    if (!isfinite(z)) continue;
                    if (z < phys_ymin) phys_ymin = z;
                    if (z > phys_ymax) phys_ymax = z;
                }
            }
            if (!isfinite(phys_ymin) || !isfinite(phys_ymax)) {
                use_z_phys_coords = 0;
                phys_ymin = pf->prob_lo[y_axis];
                phys_ymax = pf->prob_hi[y_axis];
            } else if (phys_ymax <= phys_ymin) {
                double pad = fmax(fabs(phys_ymin) * 1.0e-6, 1.0e-6);
                phys_ymin -= pad;
                phys_ymax += pad;
            }
        } else {
            phys_ymin = pf->prob_lo[y_axis];
            phys_ymax = pf->prob_hi[y_axis];
        }
    }

    /* Store current slice for mouse interaction */
    if (current_slice_data) free(current_slice_data);
    current_slice_data = (double *)malloc(width * height * sizeof(double));
    memcpy(current_slice_data, slice, width * height * sizeof(double));
    if (current_z_phys_slice) free(current_z_phys_slice);
    current_z_phys_slice = NULL;
    if (use_z_phys_coords) {
        current_z_phys_slice = (double *)malloc((size_t)width * height * sizeof(double));
        if (current_z_phys_slice)
            memcpy(current_z_phys_slice, z_phys_slice,
                   (size_t)width * height * sizeof(double));
    }
    slice_width = width;
    slice_height = height;

    /* Build mask for cells inside actual boxes (needed when level has
     * non-contiguous boxes with zero-filled gaps in between) */
    unsigned char *base_in_box = NULL;
    if (pf->current_level > 0 && pf->n_boxes > 1) {
        base_in_box = (unsigned char *)calloc(width * height, 1);
        int base_slice_coord = pf->slice_idx + pf->level_lo[pf->slice_axis];
        for (int bi = 0; bi < pf->n_boxes; bi++) {
            Box *box = &pf->boxes[bi];
            if (base_slice_coord < box->lo[pf->slice_axis] || base_slice_coord > box->hi[pf->slice_axis])
                continue;
            int dim_x, dim_y;
            if (pf->slice_axis == 2) { dim_x = 0; dim_y = 1; }
            else if (pf->slice_axis == 1) { dim_x = 0; dim_y = 2; }
            else { dim_x = 1; dim_y = 2; }
            int mi_lo = box->lo[dim_x] - pf->level_lo[dim_x];
            int mi_hi = box->hi[dim_x] - pf->level_lo[dim_x];
            int mj_lo = box->lo[dim_y] - pf->level_lo[dim_y];
            int mj_hi = box->hi[dim_y] - pf->level_lo[dim_y];
            if (mi_lo < 0) mi_lo = 0;
            if (mj_lo < 0) mj_lo = 0;
            if (mi_hi >= width) mi_hi = width - 1;
            if (mj_hi >= height) mj_hi = height - 1;
            for (int mj = mj_lo; mj <= mj_hi; mj++)
                for (int mi = mi_lo; mi <= mi_hi; mi++)
                    base_in_box[mj * width + mi] = 1;
        }
    }

    /* Find data min/max/mean, skipping gap cells when mask is active */
    for (i = 0; i < width * height; i++) {
        if (base_in_box && !base_in_box[i]) continue;
        if (slice[i] < vmin) vmin = slice[i];
        if (slice[i] > vmax) vmax = slice[i];
        vsum += slice[i];
        vcount++;
    }

    /* When overlay mode is on, include all overlay levels in min/max for consistent colorbar */
    if (pf->overlay_mode && pf->n_levels > 1) {
        LevelData *ld0 = &pf->levels[0];
        int level0_dims[3];
        double dx0_overlay[3];
        for (i = 0; i < 3; i++) {
            level0_dims[i] = (ld0->loaded && ld0->grid_dims[i] > 0) ? ld0->grid_dims[i] : pf->grid_dims[i];
            dx0_overlay[i] = (pf->prob_hi[i] - pf->prob_lo[i]) / level0_dims[i];
        }

        /* Compute current level's cell size for proper physical position */
        double dx_current_mm[3];
        for (i = 0; i < 3; i++) {
            if (pf->level_lo[i] == 0 && pf->grid_dims[i] == level0_dims[i]) {
                dx_current_mm[i] = dx0_overlay[i];
            } else {
                int curr_apparent = pf->level_lo[i] + pf->grid_dims[i];
                int curr_estimated = level0_dims[i] * pf->ref_ratio[pf->current_level > 0 ? pf->current_level : 1];
                if (curr_apparent < curr_estimated) curr_apparent = curr_estimated;
                dx_current_mm[i] = (pf->prob_hi[i] - pf->prob_lo[i]) / curr_apparent;
            }
        }

        /* Compute physical position of current slice (same formula as overlay rendering) */
        double phys_slice = pf->prob_lo[pf->slice_axis] +
                            (pf->level_lo[pf->slice_axis] + pf->slice_idx + 0.5) * dx_current_mm[pf->slice_axis];

        for (int level = pf->current_level + 1; level < pf->n_levels && level < MAX_LEVELS; level++) {
            LevelData *ld = &pf->levels[level];
            if (!ld->loaded || !ld->data) continue;

            /* Compute cell size for this level (same logic as overlay rendering) */
            double dx_lev[3];
            for (i = 0; i < 3; i++) {
                if (ld->level_lo[i] == 0 && ld->grid_dims[i] == level0_dims[i]) {
                    dx_lev[i] = dx0_overlay[i];
                } else {
                    int apparent_full_res = ld->level_lo[i] + ld->grid_dims[i];
                    int estimated_full_res = level0_dims[i] * pf->ref_ratio[level];
                    if (apparent_full_res < estimated_full_res) {
                        apparent_full_res = estimated_full_res;
                    }
                    dx_lev[i] = (pf->prob_hi[i] - pf->prob_lo[i]) / apparent_full_res;
                }
            }

            /* Compute slice index for this level using same physical position */
            int lev_slice_idx = (int)((phys_slice - pf->prob_lo[pf->slice_axis]) / dx_lev[pf->slice_axis]);
            lev_slice_idx -= ld->level_lo[pf->slice_axis];

            if (lev_slice_idx < 0 || lev_slice_idx >= ld->grid_dims[pf->slice_axis]) continue;

            /* Extract slice and find min/max */
            int lw, lh;
            if (pf->slice_axis == 2) { lw = ld->grid_dims[0]; lh = ld->grid_dims[1]; }
            else if (pf->slice_axis == 1) { lw = ld->grid_dims[0]; lh = ld->grid_dims[2]; }
            else { lw = ld->grid_dims[1]; lh = ld->grid_dims[2]; }

            double *lev_slice = (double *)malloc(lw * lh * sizeof(double));
            extract_slice_level(ld, lev_slice, pf->slice_axis, lev_slice_idx);

            /* Build mask so we only consider cells inside actual boxes, not zero-filled gaps */
            unsigned char *mm_in_box = (unsigned char *)calloc(lw * lh, 1);
            int mm_slice_coord = lev_slice_idx + ld->level_lo[pf->slice_axis];
            for (int bi = 0; bi < ld->n_boxes; bi++) {
                Box *box = &ld->boxes[bi];
                if (mm_slice_coord < box->lo[pf->slice_axis] || mm_slice_coord > box->hi[pf->slice_axis])
                    continue;
                int dim_x, dim_y;
                if (pf->slice_axis == 2) { dim_x = 0; dim_y = 1; }
                else if (pf->slice_axis == 1) { dim_x = 0; dim_y = 2; }
                else { dim_x = 1; dim_y = 2; }
                int mi_lo = box->lo[dim_x] - ld->level_lo[dim_x];
                int mi_hi = box->hi[dim_x] - ld->level_lo[dim_x];
                int mj_lo = box->lo[dim_y] - ld->level_lo[dim_y];
                int mj_hi = box->hi[dim_y] - ld->level_lo[dim_y];
                if (mi_lo < 0) mi_lo = 0;
                if (mj_lo < 0) mj_lo = 0;
                if (mi_hi >= lw) mi_hi = lw - 1;
                if (mj_hi >= lh) mj_hi = lh - 1;
                for (int mj = mj_lo; mj <= mj_hi; mj++)
                    for (int mi = mi_lo; mi <= mi_hi; mi++)
                        mm_in_box[mj * lw + mi] = 1;
            }

            for (int j = 0; j < lw * lh; j++) {
                if (!mm_in_box[j]) continue;
                if (lev_slice[j] < vmin) vmin = lev_slice[j];
                if (lev_slice[j] > vmax) vmax = lev_slice[j];
            }
            free(mm_in_box);
            free(lev_slice);
        }
    }

    /* Use custom range if set, otherwise use data min/max */
    double display_vmin, display_vmax;
    if (use_custom_range) {
        display_vmin = custom_vmin;
        display_vmax = custom_vmax;
    } else {
        display_vmin = vmin;
        display_vmax = vmax;
    }

    /* Store current vmin/vmax for colorbar */
    current_vmin = display_vmin;
    current_vmax = display_vmax;

    /* Compute mean for stats display */
    double vmean = (vcount > 0) ? vsum / vcount : 0.0;

    /* Clear canvas with white background */
    XSetForeground(display, gc, WhitePixel(display, screen));
    XFillRectangle(display, canvas, gc, 0, 0, canvas_width, canvas_height);

    /* Declare rendering variables */
    int offset_x, offset_y, local_render_width, local_render_height;

    if (pf->map_mode) {
        /* Map mode: Use appropriate geographic coordinate based on slice axis */
        
        int lon_idx = find_variable_index(pf, "lon_m");
        int lat_idx = find_variable_index(pf, "lat_m");
        
        if (lon_idx >= 0 && lat_idx >= 0) {
            /* Determine which geographic coordinate to use based on slice axis */
            double *x_geo_slice, *y_coord_slice;  /* Geographic x-axis and actual y-coordinate */
            double *x_geo_extent, *y_coord_extent;
            const char *x_label, *y_label;
            
            if (pf->slice_axis == 2) {
                /* Z-slice: longitude as x, latitude as y (normal map view) */
                x_geo_slice = (double *)malloc(width * height * sizeof(double));
                y_coord_slice = (double *)malloc(width * height * sizeof(double));
                x_geo_extent = x_geo_slice; 
                y_coord_extent = y_coord_slice;
                x_label = "lon_m"; y_label = "lat_m";
                
                int prev_var = pf->current_var;
                read_variable_data(pf, lon_idx);
                extract_slice_from_data(pf->data, pf, x_geo_slice, pf->slice_axis, pf->slice_idx);
                read_variable_data(pf, lat_idx);
                extract_slice_from_data(pf->data, pf, y_coord_slice, pf->slice_axis, pf->slice_idx);
                read_variable_data(pf, prev_var);
            } else if (pf->slice_axis == 1) {
                /* Y-slice: longitude as x, Z as y */
                x_geo_slice = (double *)malloc(width * height * sizeof(double));
                y_coord_slice = (double *)malloc(width * height * sizeof(double));
                x_geo_extent = x_geo_slice;
                y_coord_extent = y_coord_slice;
                x_label = "lon_m"; y_label = "Z";
                
                int prev_var = pf->current_var;
                read_variable_data(pf, lon_idx);
                extract_slice_from_data(pf->data, pf, x_geo_slice, pf->slice_axis, pf->slice_idx);
                read_variable_data(pf, prev_var);
                
                if (use_z_phys_coords) {
                    memcpy(y_coord_slice, z_phys_slice,
                           (size_t)width * height * sizeof(double));
                } else {
                    /* Generate uniform Z coordinates for this slice. */
                    for (j = 0; j < height; j++) {
                        for (i = 0; i < width; i++) {
                            int idx = j * width + i;
                            double z_coord = pf->prob_lo[2] + (j + 0.5) *
                                (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                            y_coord_slice[idx] = z_coord;
                        }
                    }
                }
            } else {
                /* X-slice: latitude as x, Z as y */
                x_geo_slice = (double *)malloc(width * height * sizeof(double));
                y_coord_slice = (double *)malloc(width * height * sizeof(double));
                x_geo_extent = x_geo_slice;
                y_coord_extent = y_coord_slice;
                x_label = "lat_m"; y_label = "Z";
                
                int prev_var = pf->current_var;
                read_variable_data(pf, lat_idx);
                extract_slice_from_data(pf->data, pf, x_geo_slice, pf->slice_axis, pf->slice_idx);
                read_variable_data(pf, prev_var);
                
                if (use_z_phys_coords) {
                    memcpy(y_coord_slice, z_phys_slice,
                           (size_t)width * height * sizeof(double));
                } else {
                    /* Generate uniform Z coordinates for this slice. */
                    for (j = 0; j < height; j++) {
                        for (i = 0; i < width; i++) {
                            int idx = j * width + i;
                            double z_coord = pf->prob_lo[2] + (j + 0.5) *
                                (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                            y_coord_slice[idx] = z_coord;
                        }
                    }
                }
            }
            
            /* Find actual data extent */
            double data_x_min = x_geo_extent[0], data_x_max = x_geo_extent[0];
            double data_y_min = y_coord_extent[0], data_y_max = y_coord_extent[0];
            
            for (i = 0; i < width * height; i++) {
                if (x_geo_extent[i] < data_x_min) data_x_min = x_geo_extent[i];
                if (x_geo_extent[i] > data_x_max) data_x_max = x_geo_extent[i];
                if (y_coord_extent[i] < data_y_min) data_y_min = y_coord_extent[i];
                if (y_coord_extent[i] > data_y_max) data_y_max = y_coord_extent[i];
            }
            
            /* Add small padding around data */
            double x_range = data_x_max - data_x_min;
            double y_range = data_y_max - data_y_min;
            phys_xmin = data_x_min - 0.1 * x_range;
            phys_xmax = data_x_max + 0.1 * x_range;
            phys_ymin = data_y_min - 0.1 * y_range;
            phys_ymax = data_y_max + 0.1 * y_range;

            map_last_lon_min = phys_xmin;
            map_last_lon_max = phys_xmax;
            map_last_lat_min = phys_ymin;
            map_last_lat_max = phys_ymax;
            map_has_bounds = 1;

            if (!map_auto_detected) {
                scan_coastline_files();
                auto_detect_coastlines();
                if (map_dialog_shell) {
                    for (int ci = 0; ci < n_coastlines; ci++) {
                        update_coastline_button_label(&coastlines[ci]);
                    }
                }
            }
            
            /* Set map rendering area */
            int avail_width = canvas_width - left_margin - right_margin;
            int avail_height = canvas_height - top_margin - bottom_margin;

            local_render_width = avail_width;
            local_render_height = avail_height;
            offset_x = left_margin;
            offset_y = top_margin;

            /* Save visible area for zoom interaction */
            vis_area_x = offset_x;
            vis_area_y = offset_y;
            vis_area_w = local_render_width;
            vis_area_h = local_render_height;

            /* Apply zoom */
            int zoomed_rw = (int)(local_render_width * zoom_level);
            int zoomed_rh = (int)(local_render_height * zoom_level);
            int max_sx = zoomed_rw - local_render_width;
            int max_sy = zoomed_rh - local_render_height;
            if (max_sx < 0) max_sx = 0;
            if (max_sy < 0) max_sy = 0;
            if (zoom_scroll_x > max_sx) zoom_scroll_x = max_sx;
            if (zoom_scroll_y > max_sy) zoom_scroll_y = max_sy;
            if (zoom_scroll_x < 0) zoom_scroll_x = 0;
            if (zoom_scroll_y < 0) zoom_scroll_y = 0;
            int zoom_base_x = offset_x - zoom_scroll_x;
            int zoom_base_y = offset_y - zoom_scroll_y;

            if (zoom_level > 1.0) {
                XRectangle clip = {offset_x, offset_y, local_render_width, local_render_height};
                XSetClipRectangles(display, gc, 0, 0, &clip, 1, Unsorted);
            }

            /* Create pixel data for individual points */
            unsigned long *point_pixels = (unsigned long *)malloc(width * height * sizeof(unsigned long));
            apply_colormap(slice, width, height, point_pixels, display_vmin, display_vmax, pf->colormap);

            /* Compute dot size from actual screen-space cell size so dots tile
             * without gaps regardless of grid resolution or slice axis.
             * +2 ensures overlap at boundaries. */
            int dot_w = (int)ceil((double)zoomed_rw / width) + 2;
            int dot_h = (int)ceil((double)zoomed_rh / height) + 2;
            if (dot_w < 1) dot_w = 1;
            if (dot_h < 1) dot_h = 1;

            /* Hover hit-testing must use the projected cell positions, not the
             * cell's relative i/j position in the rectangular frame. */
            prepare_map_hover_lookup();

            /* Render each data point at its coordinate */
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    int idx = j * width + i;
                    if (base_in_box && !base_in_box[idx]) continue;
                    double x_coord = x_geo_extent[idx];
                    double y_coord = y_coord_extent[idx];

                    /* Map coordinates to screen coordinates */
                    if (x_coord >= phys_xmin && x_coord <= phys_xmax && y_coord >= phys_ymin && y_coord <= phys_ymax) {
                        int screen_x = zoom_base_x + (int)((x_coord - phys_xmin) / (phys_xmax - phys_xmin) * zoomed_rw);
                        int screen_y = zoom_base_y + (int)((phys_ymax - y_coord) / (phys_ymax - phys_ymin) * zoomed_rh);

                        /* Draw a rectangle for each data point */
                        XSetForeground(display, gc, point_pixels[idx]);
                        XFillRectangle(display, canvas, gc, screen_x, screen_y, dot_w, dot_h);
                        record_map_hover_rect(screen_x, screen_y, dot_w, dot_h, idx);
                    }
                }
            }

            if (zoom_level > 1.0) XSetClipMask(display, gc, None);
            
            free(x_geo_slice);
            free(y_coord_slice);
            free(point_pixels);
        } else {
            /* Fallback to normal rendering if lon/lat not available */
            phys_xmin = pf->prob_lo[x_axis];
            phys_xmax = pf->prob_hi[x_axis];
            phys_ymin = pf->prob_lo[y_axis];
            phys_ymax = pf->prob_hi[y_axis];

            /* Use normal rendering code */
            apply_colormap(slice, width, height, pixel_data, display_vmin, display_vmax, pf->colormap);

            int avail_width = canvas_width - left_margin - right_margin;
            int avail_height = canvas_height - top_margin - bottom_margin;

            double data_aspect = (double)width / height;
            double avail_aspect = (double)avail_width / avail_height;

            if (data_aspect > avail_aspect) {
                local_render_width = avail_width;
                local_render_height = (int)(avail_width / data_aspect);
                offset_x = left_margin;
                offset_y = top_margin + (avail_height - local_render_height) / 2;
            } else {
                local_render_width = (int)(avail_height * data_aspect);
                local_render_height = avail_height;
                offset_x = left_margin + (avail_width - local_render_width) / 2;
                offset_y = top_margin;
            }

            /* Save visible area for zoom interaction */
            vis_area_x = offset_x;
            vis_area_y = offset_y;
            vis_area_w = local_render_width;
            vis_area_h = local_render_height;

            /* Apply zoom */
            int zoomed_rw = (int)(local_render_width * zoom_level);
            int zoomed_rh = (int)(local_render_height * zoom_level);
            int max_sx = zoomed_rw - local_render_width;
            int max_sy = zoomed_rh - local_render_height;
            if (max_sx < 0) max_sx = 0;
            if (max_sy < 0) max_sy = 0;
            if (zoom_scroll_x > max_sx) zoom_scroll_x = max_sx;
            if (zoom_scroll_y > max_sy) zoom_scroll_y = max_sy;
            if (zoom_scroll_x < 0) zoom_scroll_x = 0;
            if (zoom_scroll_y < 0) zoom_scroll_y = 0;
            int zoom_base_x = offset_x - zoom_scroll_x;
            int zoom_base_y = offset_y - zoom_scroll_y;

            if (zoom_level > 1.0) {
                XRectangle clip = {offset_x, offset_y, local_render_width, local_render_height};
                XSetClipRectangles(display, gc, 0, 0, &clip, 1, Unsorted);
            }

            double pixel_width = (double)zoomed_rw / width;
            double pixel_height = (double)zoomed_rh / height;

            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    if (base_in_box && !base_in_box[j * width + i]) continue;

                    unsigned long pixel = pixel_data[j * width + i];
                    XSetForeground(display, gc, pixel);

                    int x = zoom_base_x + (int)(i * pixel_width);
                    int flipped_j = height - 1 - j;
                    int y = zoom_base_y + (int)(flipped_j * pixel_height);
                    int w = (int)((i + 1) * pixel_width) - (int)(i * pixel_width);
                    int h = (int)((flipped_j + 1) * pixel_height) - (int)(flipped_j * pixel_height);
                    if (w < 1) w = 1;
                    if (h < 1) h = 1;

                    XFillRectangle(display, canvas, gc, x, y, w, h);
                }
            }

            if (zoom_level > 1.0) XSetClipMask(display, gc, None);
        }
    } else {
        /* Normal mode: apply colormap and render as regular grid */
        apply_colormap(slice, width, height, pixel_data, display_vmin, display_vmax, pf->colormap);

        /* Available area for data (excluding margins) */
        int avail_width = canvas_width - left_margin - right_margin;
        int avail_height = canvas_height - top_margin - bottom_margin;

        /* Calculate scaling to maintain aspect ratio within available area */
        double data_aspect = (double)width / height;
        double avail_aspect = (double)avail_width / avail_height;

        if (data_aspect > avail_aspect) {
            /* Width-limited */
            local_render_width = avail_width;
            local_render_height = (int)(avail_width / data_aspect);
            offset_x = left_margin;
            offset_y = top_margin + (avail_height - local_render_height) / 2;
        } else {
            /* Height-limited */
            local_render_width = (int)(avail_height * data_aspect);
            local_render_height = avail_height;
            offset_x = left_margin + (avail_width - local_render_width) / 2;
            offset_y = top_margin;
        }

        /* Save visible area for zoom interaction */
        vis_area_x = offset_x;
        vis_area_y = offset_y;
        vis_area_w = local_render_width;
        vis_area_h = local_render_height;

        /* Apply zoom */
        int zoomed_rw = (int)(local_render_width * zoom_level);
        int zoomed_rh = (int)(local_render_height * zoom_level);

        /* Clamp scroll offsets */
        int max_sx = zoomed_rw - local_render_width;
        int max_sy = zoomed_rh - local_render_height;
        if (max_sx < 0) max_sx = 0;
        if (max_sy < 0) max_sy = 0;
        if (zoom_scroll_x > max_sx) zoom_scroll_x = max_sx;
        if (zoom_scroll_y > max_sy) zoom_scroll_y = max_sy;
        if (zoom_scroll_x < 0) zoom_scroll_x = 0;
        if (zoom_scroll_y < 0) zoom_scroll_y = 0;

        /* Compute zoomed rendering base (may be off-screen) */
        int zoom_base_x = offset_x - zoom_scroll_x;
        int zoom_base_y = offset_y - zoom_scroll_y;

        /* Set clip region to visible data area when zoomed */
        if (zoom_level > 1.0) {
            XRectangle clip = {offset_x, offset_y, local_render_width, local_render_height};
            XSetClipRectangles(display, gc, 0, 0, &clip, 1, Unsorted);
        }

        /* Draw pixels as filled rectangles with correct aspect ratio */
        double pixel_width = (double)zoomed_rw / width;
        double pixel_height = (double)zoomed_rh / height;

        /* Compute visible cell range for performance optimization */
        int i_start = 0, i_end = width, j_start = 0, j_end = height;
        if (zoom_level > 1.0) {
            i_start = (int)((double)zoom_scroll_x / pixel_width);
            i_end = (int)((double)(zoom_scroll_x + local_render_width) / pixel_width) + 2;
            j_start = height - 1 - (int)((double)(zoom_scroll_y + local_render_height) / pixel_height) - 1;
            j_end = height - 1 - (int)((double)zoom_scroll_y / pixel_height) + 2;
            if (i_start < 0) i_start = 0;
            if (i_end > width) i_end = width;
            if (j_start < 0) j_start = 0;
            if (j_end > height) j_end = height;
        }

        if (use_z_phys_coords) {
            draw_z_phys_cells(z_phys_slice, pixel_data, base_in_box, width, height,
                              pf->prob_lo[x_axis], pf->prob_hi[x_axis],
                              phys_xmin, phys_xmax, phys_ymin, phys_ymax,
                              zoom_base_x, zoom_base_y, zoomed_rw, zoomed_rh);
        } else {
            for (j = j_start; j < j_end; j++) {
                for (i = i_start; i < i_end; i++) {
                    if (base_in_box && !base_in_box[j * width + i]) continue;

                    unsigned long pixel = pixel_data[j * width + i];
                    XSetForeground(display, gc, pixel);

                    int x = zoom_base_x + (int)(i * pixel_width);
                    /* Flip y-axis: higher j (higher y in data) should be at top of screen */
                    int flipped_j = height - 1 - j;
                    int y = zoom_base_y + (int)(flipped_j * pixel_height);
                    int w = (int)((i + 1) * pixel_width) - (int)(i * pixel_width);
                    int h = (int)((flipped_j + 1) * pixel_height) - (int)(flipped_j * pixel_height);
                    if (w < 1) w = 1;
                    if (h < 1) h = 1;

                    XFillRectangle(display, canvas, gc, x, y, w, h);
                }
            }
        }

        /* Reset clip for non-data elements */
        if (zoom_level > 1.0) XSetClipMask(display, gc, None);
    }

    /* Store rendering parameters for mouse interaction (use zoomed values) */
    {
        int zoomed_rw = (int)(local_render_width * zoom_level);
        int zoomed_rh = (int)(local_render_height * zoom_level);
        int zoom_base_x = offset_x - zoom_scroll_x;
        int zoom_base_y = offset_y - zoom_scroll_y;
        render_offset_x = zoom_base_x;
        render_offset_y = zoom_base_y;
        render_width = zoomed_rw;
        render_height = zoomed_rh;
        render_phys_ymin = phys_ymin;
        render_phys_ymax = phys_ymax;
        render_uses_z_phys = use_z_phys_coords && !pf->map_mode &&
                             current_z_phys_slice != NULL;
    }

    /* For overlay and subsequent rendering, use zoomed coordinate system */
    if (zoom_level > 1.0) {
        offset_x = render_offset_x;
        offset_y = render_offset_y;
        local_render_width = render_width;
        local_render_height = render_height;
        XRectangle clip = {vis_area_x, vis_area_y, vis_area_w, vis_area_h};
        XSetClipRectangles(display, gc, 0, 0, &clip, 1, Unsorted);
    }

    /* Overlay higher levels if overlay_mode is enabled */
    printf("render_slice: overlay_mode=%d, n_levels=%d, current_level=%d\n",
           pf->overlay_mode, pf->n_levels, pf->current_level);
    if (pf->overlay_mode && pf->n_levels > 1) {
        /* Level 0 cell size in physical units - MUST use Level 0's grid dimensions */
        LevelData *ld0 = &pf->levels[0];
        double dx0[3], dx_level[3];
        int level0_dims[3];
        for (i = 0; i < 3; i++) {
            /* Use Level 0 grid dims, not current level */
            level0_dims[i] = (ld0->loaded && ld0->grid_dims[i] > 0) ? ld0->grid_dims[i] : pf->grid_dims[i];
            dx0[i] = (pf->prob_hi[i] - pf->prob_lo[i]) / level0_dims[i];
        }

        /* Only overlay levels HIGHER than the current level being displayed */
        int start_level = pf->current_level + 1;
        printf("render_slice: Overlay loop from level %d to %d\n", start_level, pf->n_levels - 1);

        for (int level = start_level; level < pf->n_levels && level < MAX_LEVELS; level++) {
            LevelData *ld = &pf->levels[level];
            printf("render_slice: Level %d: loaded=%d, data=%p\n", level, ld->loaded, (void*)ld->data);
            if (!ld->loaded || !ld->data) continue;

            /* Detect per-dimension refinement by comparing grid sizes
             * If a dimension has same grid size and starts at 0, it's not refined */
            for (i = 0; i < 3; i++) {
                if (ld->level_lo[i] == 0 && ld->grid_dims[i] == level0_dims[i]) {
                    /* No refinement in this dimension - same cell size as Level 0 */
                    dx_level[i] = dx0[i];
                } else {
                    /* Refined dimension - compute actual cell size from level indices */
                    /* The full domain at this level's resolution has (level_hi_max+1) cells */
                    /* For refined dims: dx = domain_size / (level0_dims * ref_ratio) */
                    /* Approximate ref_ratio from the level bounds */
                    int apparent_full_res = ld->level_hi[i] + 1;  /* Assuming level covers to near edge */
                    if (ld->level_lo[i] > 0) {
                        /* Level doesn't start at 0, estimate full resolution */
                        apparent_full_res = ld->level_lo[i] + ld->grid_dims[i];
                    }
                    /* Use apparent full resolution to compute cell size */
                    double domain_size = pf->prob_hi[i] - pf->prob_lo[i];
                    /* Estimate based on ref_ratio if apparent_full_res seems too small */
                    int estimated_full_res = level0_dims[i] * pf->ref_ratio[level];
                    if (apparent_full_res < estimated_full_res) {
                        apparent_full_res = estimated_full_res;
                    }
                    dx_level[i] = domain_size / apparent_full_res;
                }
            }

            /* Physical bounds of this level */
            double level_phys_lo[3], level_phys_hi[3];
            for (i = 0; i < 3; i++) {
                level_phys_lo[i] = pf->prob_lo[i] + ld->level_lo[i] * dx_level[i];
                level_phys_hi[i] = pf->prob_lo[i] + (ld->level_hi[i] + 1) * dx_level[i];
            }

            /* Compute the current view level's cell size for proper slice mapping */
            double dx_current[3];
            for (i = 0; i < 3; i++) {
                /* Current level cell size: use pf->level_lo to detect if refined */
                if (pf->level_lo[i] == 0 && pf->grid_dims[i] == level0_dims[i]) {
                    dx_current[i] = dx0[i];  /* Not refined */
                } else {
                    /* Refined - compute from current level's apparent resolution */
                    int curr_apparent = pf->level_lo[i] + pf->grid_dims[i];
                    int curr_estimated = level0_dims[i] * pf->ref_ratio[pf->current_level > 0 ? pf->current_level : 1];
                    if (curr_apparent < curr_estimated) curr_apparent = curr_estimated;
                    dx_current[i] = (pf->prob_hi[i] - pf->prob_lo[i]) / curr_apparent;
                }
            }

            /* Determine level slice dimensions */
            int lwidth, lheight;
            int level_slice_idx;
            double level_x_lo, level_x_hi, level_y_lo, level_y_hi;

            /* Compute physical position of current slice using current level's cell size */
            double phys_slice_pos = pf->prob_lo[pf->slice_axis] +
                                    (pf->level_lo[pf->slice_axis] + pf->slice_idx + 0.5) * dx_current[pf->slice_axis];

            if (pf->slice_axis == 2) {  /* Z slice */
                lwidth = ld->grid_dims[0];
                lheight = ld->grid_dims[1];
                level_x_lo = level_phys_lo[0];
                level_x_hi = level_phys_hi[0];
                level_y_lo = level_phys_lo[1];
                level_y_hi = level_phys_hi[1];
                /* Map physical position to overlay level's slice index */
                level_slice_idx = (int)((phys_slice_pos - pf->prob_lo[2]) / dx_level[2]);
                level_slice_idx = level_slice_idx - ld->level_lo[2];
            } else if (pf->slice_axis == 1) {  /* Y slice */
                lwidth = ld->grid_dims[0];
                lheight = ld->grid_dims[2];
                level_x_lo = level_phys_lo[0];
                level_x_hi = level_phys_hi[0];
                level_y_lo = level_phys_lo[2];
                level_y_hi = level_phys_hi[2];
                level_slice_idx = (int)((phys_slice_pos - pf->prob_lo[1]) / dx_level[1]);
                level_slice_idx = level_slice_idx - ld->level_lo[1];
            } else {  /* X slice */
                lwidth = ld->grid_dims[1];
                lheight = ld->grid_dims[2];
                level_x_lo = level_phys_lo[1];
                level_x_hi = level_phys_hi[1];
                level_y_lo = level_phys_lo[2];
                level_y_hi = level_phys_hi[2];
                level_slice_idx = (int)((phys_slice_pos - pf->prob_lo[0]) / dx_level[0]);
                level_slice_idx = level_slice_idx - ld->level_lo[0];
            }

            /* Check if slice is within this level's bounds */
            if (level_slice_idx < 0 || level_slice_idx >= ld->grid_dims[pf->slice_axis]) {
                continue;  /* Slice not in this level */
            }

            /* Map mode overlay: draw only box boundaries using geographic coordinates
             * so they curve to follow the map projection instead of appearing as
             * straight rectangles. */
            if (pf->map_mode) {
                int lon_idx_m = find_variable_index(pf, "lon_m");
                int lat_idx_m = find_variable_index(pf, "lat_m");

                /* Determine which geo variable(s) we need for each screen axis */
                int geo_x_var = -1, geo_y_var = -1;
                int need_geo_y = 0;  /* Whether y-axis needs a geo variable (vs physical Z) */

                if (pf->slice_axis == 2) {
                    /* Z-slice: lon as x, lat as y */
                    geo_x_var = lon_idx_m;
                    geo_y_var = lat_idx_m;
                    need_geo_y = 1;
                } else if (pf->slice_axis == 1) {
                    /* Y-slice: lon as x, Z as y */
                    geo_x_var = lon_idx_m;
                } else {
                    /* X-slice: lat as x, Z as y */
                    geo_x_var = lat_idx_m;
                }

                if (geo_x_var < 0) goto skip_map_overlay;

                /* Read geo coordinate data for this level by temporarily using ld->data.
                 * Save and restore the original variable data pointer. */
                double *saved_data = ld->data;

                ld->data = NULL;
                read_variable_data_level(pf, geo_x_var, level);
                double *geo_x_3d = ld->data;

                double *geo_y_3d = NULL;
                if (need_geo_y && geo_y_var >= 0) {
                    ld->data = NULL;
                    read_variable_data_level(pf, geo_y_var, level);
                    geo_y_3d = ld->data;
                }

                ld->data = saved_data;

                if (!geo_x_3d) goto skip_map_overlay;

                /* Extract 2D geo coordinate slices */
                double *geo_x_slice = (double *)malloc(lwidth * lheight * sizeof(double));
                double *geo_y_slice = (double *)malloc(lwidth * lheight * sizeof(double));

                /* Extract x-axis geographic coordinates using temporary data swap */
                ld->data = geo_x_3d;
                extract_slice_level(ld, geo_x_slice, pf->slice_axis, level_slice_idx);

                /* Extract or compute y-axis coordinates */
                ld->data = saved_data;
                if (need_geo_y && geo_y_3d) {
                    ld->data = geo_y_3d;
                    extract_slice_level(ld, geo_y_slice, pf->slice_axis, level_slice_idx);
                } else if (use_z_phys_coords && ensure_z_phys_level_data(pf, level) == 0) {
                    ld->data = ld->z_phys_data;
                    extract_slice_level(ld, geo_y_slice, pf->slice_axis, level_slice_idx);
                } else {
                    /* Y-axis is physical Z coordinate */
                    for (int lj = 0; lj < lheight; lj++) {
                        double z_phys = pf->prob_lo[2] +
                            (ld->level_lo[2] + lj + 0.5) * dx_level[2];
                        for (int li = 0; li < lwidth; li++) {
                            geo_y_slice[lj * lwidth + li] = z_phys;
                        }
                    }
                }

                ld->data = saved_data;

                /* Compute slice_coord for box intersection test */
                int map_slice_coord = level_slice_idx + ld->level_lo[pf->slice_axis];

                /* Build in_box mask for this level */
                unsigned char *map_in_box = (unsigned char *)calloc(lwidth * lheight, 1);
                for (int bi = 0; bi < ld->n_boxes; bi++) {
                    Box *box = &ld->boxes[bi];
                    if (map_slice_coord < box->lo[pf->slice_axis] ||
                        map_slice_coord > box->hi[pf->slice_axis])
                        continue;
                    int dim_x, dim_y;
                    if (pf->slice_axis == 2) { dim_x = 0; dim_y = 1; }
                    else if (pf->slice_axis == 1) { dim_x = 0; dim_y = 2; }
                    else { dim_x = 1; dim_y = 2; }
                    int mli_lo = box->lo[dim_x] - ld->level_lo[dim_x];
                    int mli_hi = box->hi[dim_x] - ld->level_lo[dim_x];
                    int mlj_lo = box->lo[dim_y] - ld->level_lo[dim_y];
                    int mlj_hi = box->hi[dim_y] - ld->level_lo[dim_y];
                    if (mli_lo < 0) mli_lo = 0;
                    if (mlj_lo < 0) mlj_lo = 0;
                    if (mli_hi >= lwidth) mli_hi = lwidth - 1;
                    if (mlj_hi >= lheight) mlj_hi = lheight - 1;
                    for (int mj = mlj_lo; mj <= mlj_hi; mj++)
                        for (int mi = mli_lo; mi <= mli_hi; mi++)
                            map_in_box[mj * lwidth + mi] = 1;
                }

                /* Extract data slice and apply colormap */
                double *map_level_slice = (double *)malloc(lwidth * lheight * sizeof(double));
                extract_slice_level(ld, map_level_slice, pf->slice_axis, level_slice_idx);
                unsigned long *map_level_pixels = (unsigned long *)malloc(lwidth * lheight * sizeof(unsigned long));
                apply_colormap(map_level_slice, lwidth, lheight, map_level_pixels, display_vmin, display_vmax, pf->colormap);

                /* Render data points at geo coordinates (finer grid than base level) */
                for (int lj = 0; lj < lheight; lj++) {
                    for (int li = 0; li < lwidth; li++) {
                        if (!map_in_box[lj * lwidth + li]) continue;
                        int idx = lj * lwidth + li;
                        double x_coord = geo_x_slice[idx];
                        double y_coord = geo_y_slice[idx];
                        if (x_coord >= phys_xmin && x_coord <= phys_xmax &&
                            y_coord >= phys_ymin && y_coord <= phys_ymax) {
                            int sx = offset_x + (int)((x_coord - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                            int sy = offset_y + (int)((phys_ymax - y_coord) / (phys_ymax - phys_ymin) * local_render_height);
                            XSetForeground(display, gc, map_level_pixels[idx]);
                            XFillRectangle(display, canvas, gc, sx - 1, sy - 1, 3, 3);
                        }
                    }
                }

                /* Draw box boundaries as line segments following the geo coordinate grid */
                XSetForeground(display, gc, 0xFF0000);  /* Red */
                XSetLineAttributes(display, gc, 2, LineSolid, CapButt, JoinMiter);

                for (int bi = 0; bi < ld->n_boxes; bi++) {
                    Box *box = &ld->boxes[bi];
                    if (map_slice_coord < box->lo[pf->slice_axis] ||
                        map_slice_coord > box->hi[pf->slice_axis])
                        continue;

                    int dim_x, dim_y;
                    if (pf->slice_axis == 2) { dim_x = 0; dim_y = 1; }
                    else if (pf->slice_axis == 1) { dim_x = 0; dim_y = 2; }
                    else { dim_x = 1; dim_y = 2; }

                    int li_lo = box->lo[dim_x] - ld->level_lo[dim_x];
                    int li_hi = box->hi[dim_x] - ld->level_lo[dim_x];
                    int lj_lo = box->lo[dim_y] - ld->level_lo[dim_y];
                    int lj_hi = box->hi[dim_y] - ld->level_lo[dim_y];

                    /* Clamp to grid bounds */
                    if (li_lo < 0) li_lo = 0;
                    if (lj_lo < 0) lj_lo = 0;
                    if (li_hi >= lwidth) li_hi = lwidth - 1;
                    if (lj_hi >= lheight) lj_hi = lheight - 1;

                    /* Bottom edge: j = lj_lo, i varies */
                    for (int li = li_lo; li < li_hi; li++) {
                        int idx0 = lj_lo * lwidth + li;
                        int idx1 = lj_lo * lwidth + (li + 1);
                        int sx0 = offset_x + (int)((geo_x_slice[idx0] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy0 = offset_y + (int)((phys_ymax - geo_y_slice[idx0]) / (phys_ymax - phys_ymin) * local_render_height);
                        int sx1 = offset_x + (int)((geo_x_slice[idx1] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy1 = offset_y + (int)((phys_ymax - geo_y_slice[idx1]) / (phys_ymax - phys_ymin) * local_render_height);
                        XDrawLine(display, canvas, gc, sx0, sy0, sx1, sy1);
                    }

                    /* Top edge: j = lj_hi, i varies */
                    for (int li = li_lo; li < li_hi; li++) {
                        int idx0 = lj_hi * lwidth + li;
                        int idx1 = lj_hi * lwidth + (li + 1);
                        int sx0 = offset_x + (int)((geo_x_slice[idx0] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy0 = offset_y + (int)((phys_ymax - geo_y_slice[idx0]) / (phys_ymax - phys_ymin) * local_render_height);
                        int sx1 = offset_x + (int)((geo_x_slice[idx1] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy1 = offset_y + (int)((phys_ymax - geo_y_slice[idx1]) / (phys_ymax - phys_ymin) * local_render_height);
                        XDrawLine(display, canvas, gc, sx0, sy0, sx1, sy1);
                    }

                    /* Left edge: i = li_lo, j varies */
                    for (int lj = lj_lo; lj < lj_hi; lj++) {
                        int idx0 = lj * lwidth + li_lo;
                        int idx1 = (lj + 1) * lwidth + li_lo;
                        int sx0 = offset_x + (int)((geo_x_slice[idx0] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy0 = offset_y + (int)((phys_ymax - geo_y_slice[idx0]) / (phys_ymax - phys_ymin) * local_render_height);
                        int sx1 = offset_x + (int)((geo_x_slice[idx1] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy1 = offset_y + (int)((phys_ymax - geo_y_slice[idx1]) / (phys_ymax - phys_ymin) * local_render_height);
                        XDrawLine(display, canvas, gc, sx0, sy0, sx1, sy1);
                    }

                    /* Right edge: i = li_hi, j varies */
                    for (int lj = lj_lo; lj < lj_hi; lj++) {
                        int idx0 = lj * lwidth + li_hi;
                        int idx1 = (lj + 1) * lwidth + li_hi;
                        int sx0 = offset_x + (int)((geo_x_slice[idx0] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy0 = offset_y + (int)((phys_ymax - geo_y_slice[idx0]) / (phys_ymax - phys_ymin) * local_render_height);
                        int sx1 = offset_x + (int)((geo_x_slice[idx1] - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int sy1 = offset_y + (int)((phys_ymax - geo_y_slice[idx1]) / (phys_ymax - phys_ymin) * local_render_height);
                        XDrawLine(display, canvas, gc, sx0, sy0, sx1, sy1);
                    }
                }

                XSetLineAttributes(display, gc, 0, LineSolid, CapButt, JoinMiter);

                free(map_in_box);
                free(map_level_slice);
                free(map_level_pixels);
                free(geo_x_slice);
                free(geo_y_slice);
                free(geo_x_3d);
                if (geo_y_3d) free(geo_y_3d);

                printf("Overlay level %d (map mode): rendered data + boundaries, slice %d\n",
                       level, level_slice_idx);

            skip_map_overlay:
                continue;  /* Skip normal overlay rendering for this level */
            }

            /* Extract slice from this level */
            double *level_slice = (double *)malloc(lwidth * lheight * sizeof(double));
            extract_slice_level(ld, level_slice, pf->slice_axis, level_slice_idx);
            double *level_z_phys_slice = NULL;
            if (use_z_phys_coords && ensure_z_phys_level_data(pf, level) == 0) {
                level_z_phys_slice = (double *)malloc((size_t)lwidth * lheight * sizeof(double));
                if (level_z_phys_slice) {
                    double *saved_level_data = ld->data;
                    ld->data = ld->z_phys_data;
                    extract_slice_level(ld, level_z_phys_slice, pf->slice_axis, level_slice_idx);
                    ld->data = saved_level_data;
                }
            }

            /* Build mask: only render cells that fall inside an actual box.
             * Gaps between non-contiguous boxes are left unmasked (0) so the
             * underlying coarser level shows through. */
            unsigned char *in_box = (unsigned char *)calloc(lwidth * lheight, 1);
            int slice_coord = level_slice_idx + ld->level_lo[pf->slice_axis];
            for (int bi = 0; bi < ld->n_boxes; bi++) {
                Box *box = &ld->boxes[bi];
                /* Check if this box intersects the current slice */
                if (slice_coord < box->lo[pf->slice_axis] || slice_coord > box->hi[pf->slice_axis])
                    continue;
                /* Determine the 2D range this box covers in the slice plane */
                int dim_x, dim_y;  /* which 3D dims map to li, lj */
                if (pf->slice_axis == 2) { dim_x = 0; dim_y = 1; }
                else if (pf->slice_axis == 1) { dim_x = 0; dim_y = 2; }
                else { dim_x = 1; dim_y = 2; }
                int li_lo = box->lo[dim_x] - ld->level_lo[dim_x];
                int li_hi = box->hi[dim_x] - ld->level_lo[dim_x];
                int lj_lo = box->lo[dim_y] - ld->level_lo[dim_y];
                int lj_hi = box->hi[dim_y] - ld->level_lo[dim_y];
                /* Clamp to grid bounds */
                if (li_lo < 0) li_lo = 0;
                if (lj_lo < 0) lj_lo = 0;
                if (li_hi >= lwidth) li_hi = lwidth - 1;
                if (lj_hi >= lheight) lj_hi = lheight - 1;
                for (int mj = lj_lo; mj <= lj_hi; mj++) {
                    for (int mi = li_lo; mi <= li_hi; mi++) {
                        in_box[mj * lwidth + mi] = 1;
                    }
                }
            }

            /* Apply colormap to level slice */
            unsigned long *level_pixels = (unsigned long *)malloc(lwidth * lheight * sizeof(unsigned long));
            apply_colormap(level_slice, lwidth, lheight, level_pixels, display_vmin, display_vmax, pf->colormap);

            /* Map level physical bounds to screen coordinates */
            double frac_x_lo = (level_x_lo - phys_xmin) / (phys_xmax - phys_xmin);
            double frac_x_hi = (level_x_hi - phys_xmin) / (phys_xmax - phys_xmin);
            double frac_y_lo = (level_y_lo - phys_ymin) / (phys_ymax - phys_ymin);
            double frac_y_hi = (level_y_hi - phys_ymin) / (phys_ymax - phys_ymin);

            int screen_x0 = offset_x + (int)(frac_x_lo * local_render_width);
            int screen_x1 = offset_x + (int)(frac_x_hi * local_render_width);
            int screen_y0 = offset_y + local_render_height - (int)(frac_y_hi * local_render_height);
            int screen_y1 = offset_y + local_render_height - (int)(frac_y_lo * local_render_height);

            double lpixel_width = (double)(screen_x1 - screen_x0) / lwidth;
            double lpixel_height = (double)(screen_y1 - screen_y0) / lheight;

            /* Draw level pixels, skipping cells not inside any box. */
            if (level_z_phys_slice) {
                draw_z_phys_cells(level_z_phys_slice, level_pixels, in_box, lwidth, lheight,
                                  level_x_lo, level_x_hi,
                                  phys_xmin, phys_xmax, phys_ymin, phys_ymax,
                                  offset_x, offset_y, local_render_width, local_render_height);
            } else {
                for (int lj = 0; lj < lheight; lj++) {
                    for (int li = 0; li < lwidth; li++) {
                        if (!in_box[lj * lwidth + li]) continue;

                        unsigned long pixel = level_pixels[lj * lwidth + li];
                        XSetForeground(display, gc, pixel);

                        int lx = screen_x0 + (int)(li * lpixel_width);
                        int flipped_lj = lheight - 1 - lj;
                        int ly = screen_y0 + (int)(flipped_lj * lpixel_height);
                        int lw = (int)((li + 1) * lpixel_width) - (int)(li * lpixel_width);
                        int lh = (int)((flipped_lj + 1) * lpixel_height) - (int)(flipped_lj * lpixel_height);
                        if (lw < 1) lw = 1;
                        if (lh < 1) lh = 1;

                        XFillRectangle(display, canvas, gc, lx, ly, lw, lh);
                    }
                }
            }

            /* Draw box outlines for each actual box at this level */
            XSetForeground(display, gc, 0xFF0000);  /* Red */
            for (int bi = 0; bi < ld->n_boxes; bi++) {
                Box *box = &ld->boxes[bi];
                if (slice_coord < box->lo[pf->slice_axis] || slice_coord > box->hi[pf->slice_axis])
                    continue;
                int dim_x, dim_y;
                if (pf->slice_axis == 2) { dim_x = 0; dim_y = 1; }
                else if (pf->slice_axis == 1) { dim_x = 0; dim_y = 2; }
                else { dim_x = 1; dim_y = 2; }
                double box_x_lo = pf->prob_lo[dim_x] + box->lo[dim_x] * dx_level[dim_x];
                double box_x_hi = pf->prob_lo[dim_x] + (box->hi[dim_x] + 1) * dx_level[dim_x];
                double box_y_lo = pf->prob_lo[dim_y] + box->lo[dim_y] * dx_level[dim_y];
                double box_y_hi = pf->prob_lo[dim_y] + (box->hi[dim_y] + 1) * dx_level[dim_y];
                double bfx_lo = (box_x_lo - phys_xmin) / (phys_xmax - phys_xmin);
                double bfx_hi = (box_x_hi - phys_xmin) / (phys_xmax - phys_xmin);
                double bfy_lo = (box_y_lo - phys_ymin) / (phys_ymax - phys_ymin);
                double bfy_hi = (box_y_hi - phys_ymin) / (phys_ymax - phys_ymin);
                int bsx0 = offset_x + (int)(bfx_lo * render_width);
                int bsx1 = offset_x + (int)(bfx_hi * render_width);
                int bsy0 = offset_y + local_render_height - (int)(bfy_hi * local_render_height);
                int bsy1 = offset_y + local_render_height - (int)(bfy_lo * local_render_height);
                XDrawRectangle(display, canvas, gc, bsx0, bsy0, bsx1 - bsx0, bsy1 - bsy0);
            }

            free(in_box);

            free(level_slice);
            free(level_z_phys_slice);
            free(level_pixels);

            printf("Overlay level %d: slice %d, screen [%d,%d]-[%d,%d]\n",
                   level, level_slice_idx, screen_x0, screen_y0, screen_x1, screen_y1);
        }
    }

    /* Reset clip region before drawing axes and labels */
    if (zoom_level > 1.0) XSetClipMask(display, gc, None);

    /* Use visible area coordinates for axis frame/labels */
    int axis_ox = vis_area_x;
    int axis_oy = vis_area_y;
    int axis_w = vis_area_w;
    int axis_h = vis_area_h;

    /* Compute visible physical range when zoomed */
    double vis_phys_xmin = phys_xmin, vis_phys_xmax = phys_xmax;
    double vis_phys_ymin = phys_ymin, vis_phys_ymax = phys_ymax;
    if (zoom_level > 1.0) {
        int zoomed_rw = render_width;
        int zoomed_rh = render_height;
        double vis_frac_lo_x = (double)zoom_scroll_x / zoomed_rw;
        double vis_frac_hi_x = (double)(zoom_scroll_x + vis_area_w) / zoomed_rw;
        double vis_frac_lo_y = (double)zoom_scroll_y / zoomed_rh;
        double vis_frac_hi_y = (double)(zoom_scroll_y + vis_area_h) / zoomed_rh;
        vis_phys_xmin = phys_xmin + vis_frac_lo_x * (phys_xmax - phys_xmin);
        vis_phys_xmax = phys_xmin + vis_frac_hi_x * (phys_xmax - phys_xmin);
        /* Y is flipped: top of screen = high y, bottom = low y */
        vis_phys_ymin = phys_ymin + (1.0 - vis_frac_hi_y) * (phys_ymax - phys_ymin);
        vis_phys_ymax = phys_ymin + (1.0 - vis_frac_lo_y) * (phys_ymax - phys_ymin);
    }

    /* Draw axis frame (border around data) */
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    XDrawRectangle(display, canvas, text_gc, axis_ox, axis_oy, axis_w, axis_h);

    /* Draw X-axis ticks and labels */
    int n_xticks = 5;
    char label[32];
    for (i = 0; i <= n_xticks; i++) {
        double frac = (double)i / n_xticks;
        int tick_x = axis_ox + (int)(frac * axis_w);
        double phys_val = vis_phys_xmin + frac * (vis_phys_xmax - vis_phys_xmin);

        /* Draw tick mark */
        XDrawLine(display, canvas, text_gc, tick_x, axis_oy + axis_h,
              tick_x, axis_oy + axis_h + 5);

        /* Draw label */
        snprintf(label, sizeof(label), "%.3g", phys_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(display, canvas, text_gc, tick_x - label_width / 2,
                axis_oy + axis_h + 18, label, strlen(label));
    }

    /* Draw Y-axis ticks and labels */
    int n_yticks = 5;
    for (i = 0; i <= n_yticks; i++) {
        double frac = (double)i / n_yticks;
        int tick_y = axis_oy + axis_h - (int)(frac * axis_h);
        double phys_val = vis_phys_ymin + frac * (vis_phys_ymax - vis_phys_ymin);

        /* Draw tick mark */
        XDrawLine(display, canvas, text_gc, axis_ox - 5, tick_y, axis_ox, tick_y);

        /* Draw label */
        snprintf(label, sizeof(label), "%.3g", phys_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(display, canvas, text_gc, axis_ox - label_width - 8,
                    tick_y + 4, label, strlen(label));
    }

    /* Draw axis labels with units */
    const char *axis_names[] = {"X", "Y", "Z"};
    char x_label[32], y_label[32];

    if (pf->map_mode) {
        /* Map mode: axis labels depend on slice axis */
        if (pf->slice_axis == 2) {
            /* Z-slice: lon on X, lat on Y */
            strcpy(x_label, "Longitude (deg)");
            strcpy(y_label, "Latitude (deg)");
        } else if (pf->slice_axis == 1) {
            /* Y-slice: lon on X, height on Y */
            strcpy(x_label, "Longitude (deg)");
            strcpy(y_label, use_z_phys_coords ? "z_phys (m)" : "Height (m)");
        } else {
            /* X-slice: lat on X, height on Y */
            strcpy(x_label, "Latitude (deg)");
            strcpy(y_label, use_z_phys_coords ? "z_phys (m)" : "Height (m)");
        }
    } else {
        /* Normal mode: use physical coordinates with units */
        const char *unit_str = "(m)";
        snprintf(x_label, sizeof(x_label), "%s %s", axis_names[x_axis], unit_str);
        if (use_z_phys_coords)
            snprintf(y_label, sizeof(y_label), "z_phys %s", unit_str);
        else
            snprintf(y_label, sizeof(y_label), "%s %s", axis_names[y_axis], unit_str);
    }

    /* X-axis label (centered below ticks) */
    int xlabel_width = XTextWidth(font, x_label, strlen(x_label));
    XDrawString(display, canvas, text_gc,
                axis_ox + axis_w / 2 - xlabel_width / 2,
                axis_oy + axis_h + 35, x_label, strlen(x_label));

    /* Y-axis label (rotated 90° CCW) */
    int y_label_x = axis_ox - left_margin;
    draw_y_label_ccw(y_label, y_label_x, axis_oy + axis_h / 2);

    /* Draw text overlay - show display range (custom if set) */
    if (use_custom_range) {
        snprintf(stats_text, sizeof(stats_text), "range: %.3e to %.3e (custom)  mean: %.3e", display_vmin, display_vmax, vmean);
    } else {
        snprintf(stats_text, sizeof(stats_text), "min: %.3e  max: %.3e  mean: %.3e", display_vmin, display_vmax, vmean);
    }
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    XSetBackground(display, text_gc, WhitePixel(display, screen));
    XDrawImageString(display, canvas, text_gc, left_margin, canvas_height - 5,
                    stats_text, strlen(stats_text));

    /* Draw colorbar */
    draw_colorbar(display_vmin, display_vmax, pf->colormap,
                  pf->variables[pf->current_var]);

    /* Set clip region for quiver/map overlays when zoomed */
    if (zoom_level > 1.0) {
        XRectangle clip = {vis_area_x, vis_area_y, vis_area_w, vis_area_h};
        XSetClipRectangles(display, gc, 0, 0, &clip, 1, Unsorted);
    }

    /* Draw quiver overlay if enabled */
    if (quiver_data.enabled) {
        render_quiver_overlay(pf);
    }

    /* Draw map overlay (coastlines etc.) only for Z-slice in map mode */
    if (pf->map_mode && pf->slice_axis == 2) {
        render_map_overlay(pf, phys_xmin, phys_xmax, phys_ymin, phys_ymax);
    }

    /* Reset clip after overlays */
    if (zoom_level > 1.0) XSetClipMask(display, gc, None);

    XFlush(display);
    
    printf("Rendered: %s, slice %d/%d (%.3e to %.3e)\n", 
           pf->variables[pf->current_var], pf->slice_idx + 1,
           pf->grid_dims[pf->slice_axis], vmin, vmax);
    
    free(slice);
    free(z_phys_slice);
    if (base_in_box) free(base_in_box);
}

/* Clamp zoom scroll offsets to valid range */
void clamp_zoom_scroll(void) {
    int zoomed_w = (int)(vis_area_w * zoom_level);
    int zoomed_h = (int)(vis_area_h * zoom_level);
    int max_sx = zoomed_w - vis_area_w;
    int max_sy = zoomed_h - vis_area_h;
    if (max_sx < 0) max_sx = 0;
    if (max_sy < 0) max_sy = 0;
    if (zoom_scroll_x > max_sx) zoom_scroll_x = max_sx;
    if (zoom_scroll_y > max_sy) zoom_scroll_y = max_sy;
    if (zoom_scroll_x < 0) zoom_scroll_x = 0;
    if (zoom_scroll_y < 0) zoom_scroll_y = 0;
}

/* Reset zoom to 1.0 */
void zoom_reset(void) {
    zoom_level = 1.0;
    zoom_scroll_x = 0;
    zoom_scroll_y = 0;
    zoom_dragging = 0;
}

/* Zoom in button callback */
void zoom_in_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pf) return;
    if (zoom_level >= 50.0) return;
    double old_zoom = zoom_level;
    zoom_level *= 1.5;
    if (zoom_level > 50.0) zoom_level = 50.0;
    /* Keep center of visible area centered */
    double cx = (zoom_scroll_x + vis_area_w / 2.0) / (vis_area_w * old_zoom);
    double cy = (zoom_scroll_y + vis_area_h / 2.0) / (vis_area_h * old_zoom);
    zoom_scroll_x = (int)(cx * vis_area_w * zoom_level - vis_area_w / 2.0);
    zoom_scroll_y = (int)(cy * vis_area_h * zoom_level - vis_area_h / 2.0);
    clamp_zoom_scroll();
    render_slice(global_pf);
}

/* Zoom out button callback */
void zoom_out_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pf) return;
    if (zoom_level <= 1.0) return;
    double old_zoom = zoom_level;
    zoom_level /= 1.5;
    if (zoom_level < 1.0) zoom_level = 1.0;
    /* Keep center of visible area centered */
    double cx = (zoom_scroll_x + vis_area_w / 2.0) / (vis_area_w * old_zoom);
    double cy = (zoom_scroll_y + vis_area_h / 2.0) / (vis_area_h * old_zoom);
    zoom_scroll_x = (int)(cx * vis_area_w * zoom_level - vis_area_w / 2.0);
    zoom_scroll_y = (int)(cy * vis_area_h * zoom_level - vis_area_h / 2.0);
    clamp_zoom_scroll();
    render_slice(global_pf);
}

/* Mouse motion handler - show value at cursor */
void canvas_motion_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (!global_pf || !current_slice_data) return;

    int mouse_x = event->xmotion.x;
    int mouse_y = event->xmotion.y;

    /* Handle drag / click detection */
    if (zoom_dragging) {
        /* Button released — detect click vs drag */
        if (!(event->xmotion.state & Button1Mask)) {
            int dx = mouse_x - zoom_drag_start_x;
            int dy = mouse_y - zoom_drag_start_y;
            zoom_dragging = 0;

            /* Small movement = click: show line profiles */
            if (abs(dx) < 5 && abs(dy) < 5) {
                if (mouse_x >= vis_area_x && mouse_x < vis_area_x + vis_area_w &&
                    mouse_y >= vis_area_y && mouse_y < vis_area_y + vis_area_h) {
                    int data_x, data_y;
                    if (canvas_to_data_indices(mouse_x, mouse_y, &data_x, &data_y)) {
                        show_line_profiles(global_pf, data_x, data_y);
                    }
                }
            }
        } else if (zoom_level > 1.0) {
            /* Button still held and zoomed: pan */
            int dx = mouse_x - zoom_drag_start_x;
            int dy = mouse_y - zoom_drag_start_y;
            zoom_scroll_x = zoom_drag_scroll_x0 - dx;
            zoom_scroll_y = zoom_drag_scroll_y0 - dy;
            clamp_zoom_scroll();
            render_slice(global_pf);
            return;
        }
    }

    /* Check if cursor is within visible data area */
    if (mouse_x < vis_area_x || mouse_x >= vis_area_x + vis_area_w ||
        mouse_y < vis_area_y || mouse_y >= vis_area_y + vis_area_h) {
        /* Outside data region - clear hover text */
        if (hover_value_text[0] != '\0') {
            hover_value_text[0] = '\0';
            update_info_label(global_pf);
        }
        return;
    }

    /* Convert mouse coordinates to the regular or terrain-following cell. */
    int data_x, data_y;
    if (canvas_to_data_indices(mouse_x, mouse_y, &data_x, &data_y)) {
        double value = current_slice_data[data_y * slice_width + data_x];

        /* Update hover value text and info label */
        snprintf(hover_value_text, sizeof(hover_value_text), "[%d,%d]: %.6e", data_x, data_y, value);
        update_info_label(global_pf);
    } else if (hover_value_text[0] != '\0') {
        /* In map mode, blank areas inside the rectangular frame may not
         * correspond to any projected scalar cell. */
        hover_value_text[0] = '\0';
        update_info_label(global_pf);
    }
}

/* Mouse button handler - show line profiles through clicked point */
void canvas_button_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (!global_pf || !current_slice_data) return;

    /* Only process events on the canvas window */
    if (event->xbutton.window != canvas) return;

    /* Set keyboard focus to canvas - needed for remote X11 forwarding */
    XSetInputFocus(display, canvas, RevertToParent, CurrentTime);

    int mouse_x = event->xbutton.x;
    int mouse_y = event->xbutton.y;

    /* Mouse wheel zoom */
    if (event->xbutton.button == Button4 || event->xbutton.button == Button5) {
        /* Only zoom if cursor is within visible data area */
        if (mouse_x < vis_area_x || mouse_x >= vis_area_x + vis_area_w ||
            mouse_y < vis_area_y || mouse_y >= vis_area_y + vis_area_h) {
            return;
        }

        double old_zoom = zoom_level;
        double factor = (event->xbutton.button == Button4) ? 1.25 : (1.0 / 1.25);
        double new_zoom = zoom_level * factor;
        if (new_zoom < 1.0) new_zoom = 1.0;
        if (new_zoom > 50.0) new_zoom = 50.0;
        if (new_zoom == old_zoom) return;
        zoom_level = new_zoom;

        /* Adjust scroll so the point under cursor stays fixed */
        double fx = (mouse_x - vis_area_x + zoom_scroll_x) / (vis_area_w * old_zoom);
        double fy = (mouse_y - vis_area_y + zoom_scroll_y) / (vis_area_h * old_zoom);
        zoom_scroll_x = (int)(fx * vis_area_w * zoom_level - (mouse_x - vis_area_x));
        zoom_scroll_y = (int)(fy * vis_area_h * zoom_level - (mouse_y - vis_area_y));
        clamp_zoom_scroll();
        render_slice(global_pf);
        return;
    }

    if (event->xbutton.button != Button1) return;

    /* Always record drag start (used for both zoomed drag-pan and click detection) */
    zoom_dragging = 1;
    zoom_drag_start_x = mouse_x;
    zoom_drag_start_y = mouse_y;
    zoom_drag_scroll_x0 = zoom_scroll_x;
    zoom_drag_scroll_y0 = zoom_scroll_y;
}

/* Mouse button release handler - end drag pan, detect click */
void canvas_button_release_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (!global_pf || event->xbutton.button != Button1) return;
    if (event->xbutton.window != canvas) return;

    if (zoom_dragging) {
        int mouse_x = event->xbutton.x;
        int mouse_y = event->xbutton.y;
        int dx = mouse_x - zoom_drag_start_x;
        int dy = mouse_y - zoom_drag_start_y;
        zoom_dragging = 0;

        /* If barely moved, treat as click for line profile */
        if (abs(dx) < 5 && abs(dy) < 5) {
            if (mouse_x < vis_area_x || mouse_x >= vis_area_x + vis_area_w ||
                mouse_y < vis_area_y || mouse_y >= vis_area_y + vis_area_h) {
                return;
            }
            int data_x, data_y;
            if (canvas_to_data_indices(mouse_x, mouse_y, &data_x, &data_y)) {
                show_line_profiles(global_pf, data_x, data_y);
            }
        }
    }
}

/* Draw a line plot on a window */
void draw_line_plot(Display *dpy, Window win, GC plot_gc, double *data, double *x_values,
                   int n_points, int width, int height, double vmin, double vmax, 
                   double xmin, double xmax, const char *title, const char *xlabel) {
    /* Clear background */
    XSetForeground(dpy, plot_gc, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, plot_gc, 0, 0, width, height);
    
    /* Draw border */
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, plot_gc, 0, 0, width - 1, height - 1);
    
    /* Draw title */
    if (font) {
        XSetFont(dpy, plot_gc, font->fid);
        XDrawString(dpy, win, plot_gc, 10, 20, title, strlen(title));
    }
    
    /* Plot area */
    int plot_left = 50;
    int plot_right = width - 20;
    int plot_top = 40;
    int plot_bottom = height - 45;  /* More space for x-axis labels */
    int plot_width = plot_right - plot_left;
    int plot_height = plot_bottom - plot_top;
    
    if (plot_width <= 0 || plot_height <= 0 || n_points < 2) return;
    
    /* Draw axes */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_bottom, plot_right, plot_bottom);  /* x-axis */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_top, plot_left, plot_bottom);      /* y-axis */
    
    /* Draw y-axis ticks and labels */
    char label[64];
    int num_y_ticks = 4;
    for (int i = 0; i <= num_y_ticks; i++) {
        double y_val = vmin + (vmax - vmin) * i / num_y_ticks;
        int y_pos = plot_bottom - (int)(plot_height * i / num_y_ticks);
        
        /* Draw tick mark */
        XDrawLine(dpy, win, plot_gc, plot_left - 3, y_pos, plot_left, y_pos);
        
        /* Draw label */
        snprintf(label, sizeof(label), "%.2e", y_val);
        XDrawString(dpy, win, plot_gc, 5, y_pos + 4, label, strlen(label));
    }
    
    /* Draw x-axis ticks and labels */
    int num_x_ticks = 10;
    for (int i = 0; i <= num_x_ticks; i++) {
        double x_val = xmin + (xmax - xmin) * i / num_x_ticks;
        int x_pos = plot_left + (int)(plot_width * i / num_x_ticks);
        
        /* Draw tick mark */
        XDrawLine(dpy, win, plot_gc, x_pos, plot_bottom, x_pos, plot_bottom + 3);
        
        /* Draw label */
        snprintf(label, sizeof(label), "%.0f", x_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(dpy, win, plot_gc, x_pos - label_width / 2, plot_bottom + 14, label, strlen(label));
    }
    
    /* Draw x-axis label */
    if (xlabel && xlabel[0]) {
        int xlabel_width = XTextWidth(font, xlabel, strlen(xlabel));
        XDrawString(dpy, win, plot_gc, plot_left + (plot_width - xlabel_width) / 2, 
                   plot_bottom + 28, xlabel, strlen(xlabel));
    }
    
    /* Draw line plot */
    XSetForeground(dpy, plot_gc, 0x0000FF);  /* Blue */
    double range = vmax - vmin;
    if (range == 0) range = 1;
    double xrange = xmax - xmin;
    if (xrange == 0) xrange = 1;
    
    for (int i = 0; i < n_points - 1; i++) {
        int x1 = plot_left + (int)((x_values[i] - xmin) / xrange * plot_width);
        int x2 = plot_left + (int)((x_values[i + 1] - xmin) / xrange * plot_width);
        int y1 = plot_bottom - (int)((data[i] - vmin) / range * plot_height);
        int y2 = plot_bottom - (int)((data[i + 1] - vmin) / range * plot_height);
        
        /* Clamp to plot area */
        if (y1 < plot_top) y1 = plot_top;
        if (y1 > plot_bottom) y1 = plot_bottom;
        if (y2 < plot_top) y2 = plot_top;
        if (y2 > plot_bottom) y2 = plot_bottom;
        
        XDrawLine(dpy, win, plot_gc, x1, y1, x2, y2);
    }

    XFlush(dpy);
}

/* Draw a horizontal line plot (layer on Y axis, values on X axis) */
void draw_horizontal_plot(Display *dpy, Window win, GC plot_gc, double *data, double *y_values,
                          int n_points, int width, int height, double vmin, double vmax,
                          double ymin, double ymax, const char *title, const char *ylabel,
                          const char *vlabel) {
    /* Clear background */
    XSetForeground(dpy, plot_gc, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, plot_gc, 0, 0, width, height);

    /* Draw border */
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, plot_gc, 0, 0, width - 1, height - 1);

    /* Draw title */
    if (font) {
        XSetFont(dpy, plot_gc, font->fid);
        XDrawString(dpy, win, plot_gc, 10, 20, title, strlen(title));
    }

    /* Plot area - more space on left for Y axis labels, bottom for X label */
    int plot_left = 60;
    int plot_right = width - 20;
    int plot_top = 40;
    int plot_bottom = height - 55;
    int plot_width = plot_right - plot_left;
    int plot_height = plot_bottom - plot_top;

    if (plot_width <= 0 || plot_height <= 0 || n_points < 2) return;

    /* Draw axes */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_bottom, plot_right, plot_bottom);  /* x-axis */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_top, plot_left, plot_bottom);      /* y-axis */

    /* Draw x-axis (value) ticks and labels */
    char label[64];
    int num_x_ticks = 4;
    for (int i = 0; i <= num_x_ticks; i++) {
        double x_val = vmin + (vmax - vmin) * i / num_x_ticks;
        int x_pos = plot_left + (int)(plot_width * i / num_x_ticks);

        /* Draw tick mark */
        XDrawLine(dpy, win, plot_gc, x_pos, plot_bottom, x_pos, plot_bottom + 3);

        /* Draw label */
        snprintf(label, sizeof(label), "%.2e", x_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(dpy, win, plot_gc, x_pos - label_width / 2, plot_bottom + 14, label, strlen(label));
    }

    /* Draw y-axis (layer) ticks and labels */
    int num_y_ticks = 5;
    if (n_points < num_y_ticks) num_y_ticks = n_points - 1;
    for (int i = 0; i <= num_y_ticks; i++) {
        double y_val = ymin + (ymax - ymin) * i / num_y_ticks;
        /* Y axis is inverted: higher layer values at top */
        int y_pos = plot_bottom - (int)(plot_height * i / num_y_ticks);

        /* Draw tick mark */
        XDrawLine(dpy, win, plot_gc, plot_left - 3, y_pos, plot_left, y_pos);

        /* Draw label — %.6g avoids scientific notation for values up to ~999999 */
        snprintf(label, sizeof(label), "%.6g", y_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(dpy, win, plot_gc, plot_left - label_width - 5, y_pos + 4, label, strlen(label));
    }

    /* Draw y-axis label (rotated would be ideal but just put at top) */
    if (ylabel && ylabel[0]) {
        XDrawString(dpy, win, plot_gc, 5, plot_top - 5, ylabel, strlen(ylabel));
    }

    /* Draw x-axis label (value label) centered below the plot */
    if (vlabel && vlabel[0]) {
        int vlabel_width = XTextWidth(font, vlabel, strlen(vlabel));
        XDrawString(dpy, win, plot_gc, plot_left + (plot_width - vlabel_width) / 2,
                    plot_bottom + 30, vlabel, strlen(vlabel));
    }

    /* Draw horizontal line plot (layer on Y, value on X) */
    XSetForeground(dpy, plot_gc, 0x0000FF);  /* Blue */
    double xrange = vmax - vmin;
    if (xrange == 0) xrange = 1;
    double yrange = ymax - ymin;
    if (yrange == 0) yrange = 1;

    for (int i = 0; i < n_points - 1; i++) {
        /* X position based on data value */
        int x1 = plot_left + (int)((data[i] - vmin) / xrange * plot_width);
        int x2 = plot_left + (int)((data[i + 1] - vmin) / xrange * plot_width);
        /* Y position based on layer (inverted: higher layer at top) */
        int y1 = plot_bottom - (int)((y_values[i] - ymin) / yrange * plot_height);
        int y2 = plot_bottom - (int)((y_values[i + 1] - ymin) / yrange * plot_height);

        /* Clamp to plot area */
        if (x1 < plot_left) x1 = plot_left;
        if (x1 > plot_right) x1 = plot_right;
        if (x2 < plot_left) x2 = plot_left;
        if (x2 > plot_right) x2 = plot_right;

        XDrawLine(dpy, win, plot_gc, x1, y1, x2, y2);
    }

    XFlush(dpy);
}

/* Expose event handler for horizontal plot canvas */
void horizontal_plot_expose_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != Expose) return;

    PlotData *plot_data = (PlotData *)client_data;
    if (!plot_data || !plot_data->data) return;

    Window win = XtWindow(w);
    if (!win) return;

    Dimension width, height;
    XtVaGetValues(w, XtNwidth, &width, XtNheight, &height, NULL);

    GC plot_gc = XCreateGC(display, win, 0, NULL);
    /* Note: For horizontal plot, x_values are actually layer indices (Y axis) */
    /* vmin/vmax are value range (X axis), xmin/xmax are layer range (Y axis) */
    draw_horizontal_plot(display, win, plot_gc, plot_data->data, plot_data->x_values,
                         plot_data->n_points, width, height, plot_data->vmin, plot_data->vmax,
                         plot_data->xmin, plot_data->xmax, plot_data->title, plot_data->xlabel,
                         plot_data->vlabel);
    XFreeGC(display, plot_gc);
}

/* Expose event handler for plot canvas */
void plot_expose_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != Expose) return;
    
    PlotData *plot_data = (PlotData *)client_data;
    if (!plot_data || !plot_data->data) return;
    
    Window win = XtWindow(w);
    if (!win) return;
    
    Dimension width, height;
    XtVaGetValues(w, XtNwidth, &width, XtNheight, &height, NULL);
    
    GC plot_gc = XCreateGC(display, win, 0, NULL);
    draw_line_plot(display, win, plot_gc, plot_data->data, plot_data->x_values,
                   plot_data->n_points, width, height, plot_data->vmin, plot_data->vmax,
                   plot_data->xmin, plot_data->xmax, plot_data->title, plot_data->xlabel);
    XFreeGC(display, plot_gc);
}

/* Line-profile popup data (3 standard line plots: X, Y, Z profiles) */
typedef struct {
    Widget shell;
    Widget canvases[3];     /* x_canvas, y_canvas, z_canvas */
    PlotData *plots[3];     /* x_plot, y_plot, z_plot */
    double *phys_values[3]; /* physical coordinate arrays (m) per axis */
    double *layer_values[3];/* 0-indexed layer arrays per axis */
    double phys_min[3], phys_max[3];
    double layer_min[3], layer_max[3];
    int show_layer;         /* 0 = physical (default), 1 = layer index */
    char phys_labels[3][32];
    char layer_labels[3][32];
} LineProfilePopupData;

/* Close callback for line-profile popup */
void close_line_profile_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    LineProfilePopupData *pd = (LineProfilePopupData *)client_data;
    if (pd) {
        for (int i = 0; i < 3; i++) {
            if (pd->plots[i]) {
                if (pd->plots[i]->data) free(pd->plots[i]->data);
                /* x_values is owned by phys_values/layer_values, not the plot */
                free(pd->plots[i]);
            }
            if (pd->phys_values[i])  free(pd->phys_values[i]);
            if (pd->layer_values[i]) free(pd->layer_values[i]);
        }
        XtDestroyWidget(pd->shell);
        free(pd);
    }
}

/* Toggle physical/layer axis on all three line-profile plots */
void toggle_line_layer_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    LineProfilePopupData *pd = (LineProfilePopupData *)client_data;
    pd->show_layer = !pd->show_layer;

    for (int i = 0; i < 3; i++) {
        pd->plots[i]->x_values = pd->show_layer ? pd->layer_values[i] : pd->phys_values[i];
        pd->plots[i]->xmin     = pd->show_layer ? pd->layer_min[i]    : pd->phys_min[i];
        pd->plots[i]->xmax     = pd->show_layer ? pd->layer_max[i]    : pd->phys_max[i];
        snprintf(pd->plots[i]->xlabel, sizeof(pd->plots[i]->xlabel), "%s",
                 pd->show_layer ? pd->layer_labels[i] : pd->phys_labels[i]);
        if (XtWindow(pd->canvases[i]))
            XClearArea(XtDisplay(pd->canvases[i]), XtWindow(pd->canvases[i]), 0, 0, 0, 0, True);
    }
    XtVaSetValues(w, XtNlabel, pd->show_layer ? "m" : "Layer", NULL);
}

/* Show 1D line profiles through clicked point along x, y, z */
void show_line_profiles(PlotfileData *pf, int data_x, int data_y) {
    /* Get 3D coordinates based on current slice */
    int x_coord, y_coord, z_coord;

    if (pf->slice_axis == 2) {
        x_coord = data_x; y_coord = data_y; z_coord = pf->slice_idx;
    } else if (pf->slice_axis == 1) {
        x_coord = data_x; y_coord = pf->slice_idx; z_coord = data_y;
    } else {
        x_coord = pf->slice_idx; y_coord = data_x; z_coord = data_y;
    }

    /* Compute physical coordinate arrays for all three axes */
    int dims[3] = { pf->grid_dims[0], pf->grid_dims[1], pf->grid_dims[2] };
    double *phys[3], *layer[3];
    for (int a = 0; a < 3; a++) {
        phys[a]  = (double *)malloc(dims[a] * sizeof(double));
        layer[a] = (double *)malloc(dims[a] * sizeof(double));
        double dphys = (pf->prob_hi[a] - pf->prob_lo[a]) / dims[a];
        for (int n = 0; n < dims[a]; n++) {
            phys[a][n]  = pf->prob_lo[a] + (n + 0.5) * dphys;
            layer[a][n] = n;
        }
    }
    if (pf->use_z_phys && ensure_z_phys_data(pf) == 0) {
        for (int k = 0; k < dims[2]; k++) {
            int idx = k * dims[0] * dims[1] + y_coord * dims[0] + x_coord;
            phys[2][k] = pf->z_phys_data[idx];
        }
    }

    /* Create plot data — default x axis = physical (m) */
    PlotData *x_plot_data = (PlotData *)malloc(sizeof(PlotData));
    x_plot_data->n_points = dims[0];
    x_plot_data->data     = (double *)malloc(dims[0] * sizeof(double));
    x_plot_data->x_values = phys[0];
    x_plot_data->vmin = 1e30; x_plot_data->vmax = -1e30;
    for (int i = 0; i < dims[0]; i++) {
        int idx = z_coord * dims[0] * dims[1] + y_coord * dims[0] + i;
        x_plot_data->data[i] = pf->data[idx];
        if (x_plot_data->data[i] < x_plot_data->vmin) x_plot_data->vmin = x_plot_data->data[i];
        if (x_plot_data->data[i] > x_plot_data->vmax) x_plot_data->vmax = x_plot_data->data[i];
    }
    x_plot_data->xmin = phys[0][0];
    x_plot_data->xmax = phys[0][dims[0]-1];
    snprintf(x_plot_data->title, sizeof(x_plot_data->title), "%s along X (Y=%d, Z=%d)",
             pf->variables[pf->current_var], y_coord, z_coord);
    snprintf(x_plot_data->xlabel, sizeof(x_plot_data->xlabel), "X (m)");

    PlotData *y_plot_data = (PlotData *)malloc(sizeof(PlotData));
    y_plot_data->n_points = dims[1];
    y_plot_data->data     = (double *)malloc(dims[1] * sizeof(double));
    y_plot_data->x_values = phys[1];
    y_plot_data->vmin = 1e30; y_plot_data->vmax = -1e30;
    for (int j = 0; j < dims[1]; j++) {
        int idx = z_coord * dims[0] * dims[1] + j * dims[0] + x_coord;
        y_plot_data->data[j] = pf->data[idx];
        if (y_plot_data->data[j] < y_plot_data->vmin) y_plot_data->vmin = y_plot_data->data[j];
        if (y_plot_data->data[j] > y_plot_data->vmax) y_plot_data->vmax = y_plot_data->data[j];
    }
    y_plot_data->xmin = phys[1][0];
    y_plot_data->xmax = phys[1][dims[1]-1];
    snprintf(y_plot_data->title, sizeof(y_plot_data->title), "%s along Y (X=%d, Z=%d)",
             pf->variables[pf->current_var], x_coord, z_coord);
    snprintf(y_plot_data->xlabel, sizeof(y_plot_data->xlabel), "Y (m)");

    PlotData *z_plot_data = (PlotData *)malloc(sizeof(PlotData));
    z_plot_data->n_points = dims[2];
    z_plot_data->data     = (double *)malloc(dims[2] * sizeof(double));
    z_plot_data->x_values = phys[2];
    z_plot_data->vmin = 1e30; z_plot_data->vmax = -1e30;
    for (int k = 0; k < dims[2]; k++) {
        int idx = k * dims[0] * dims[1] + y_coord * dims[0] + x_coord;
        z_plot_data->data[k] = pf->data[idx];
        if (z_plot_data->data[k] < z_plot_data->vmin) z_plot_data->vmin = z_plot_data->data[k];
        if (z_plot_data->data[k] > z_plot_data->vmax) z_plot_data->vmax = z_plot_data->data[k];
    }
    z_plot_data->xmin = phys[2][0];
    z_plot_data->xmax = phys[2][dims[2]-1];
    snprintf(z_plot_data->title, sizeof(z_plot_data->title), "%s along Z (X=%d, Y=%d)",
             pf->variables[pf->current_var], x_coord, y_coord);
    snprintf(z_plot_data->xlabel, sizeof(z_plot_data->xlabel), "%s",
             pf->use_z_phys && pf->z_phys_data ? "z_phys (m)" : "Height (m)");

    /* Build LineProfilePopupData */
    LineProfilePopupData *pd = (LineProfilePopupData *)malloc(sizeof(LineProfilePopupData));
    pd->plots[0] = x_plot_data;
    pd->plots[1] = y_plot_data;
    pd->plots[2] = z_plot_data;
    for (int a = 0; a < 3; a++) {
        pd->phys_values[a]  = phys[a];
        pd->layer_values[a] = layer[a];
        pd->phys_min[a]     = phys[a][0];
        pd->phys_max[a]     = phys[a][dims[a]-1];
        pd->layer_min[a]    = 0;
        pd->layer_max[a]    = dims[a] - 1;
    }
    pd->show_layer = 0;
    snprintf(pd->phys_labels[0], sizeof(pd->phys_labels[0]), "X (m)");
    snprintf(pd->phys_labels[1], sizeof(pd->phys_labels[1]), "Y (m)");
    snprintf(pd->phys_labels[2], sizeof(pd->phys_labels[2]), "%s",
             pf->use_z_phys && pf->z_phys_data ? "z_phys (m)" : "Height (m)");
    snprintf(pd->layer_labels[0], sizeof(pd->layer_labels[0]), "X Layer");
    snprintf(pd->layer_labels[1], sizeof(pd->layer_labels[1]), "Y Layer");
    snprintf(pd->layer_labels[2], sizeof(pd->layer_labels[2]), "Z Layer");

    /* Create popup shell */
    Widget popup_shell = XtVaCreatePopupShell("Line Profiles",
        transientShellWidgetClass, toplevel,
        XtNwidth, 900,
        XtNheight, 700,
        NULL);
    pd->shell = popup_shell;

    Widget popup_form = XtVaCreateManagedWidget("form",
        formWidgetClass, popup_shell,
        NULL);

    /* Title label */
    char title_text[256];
    snprintf(title_text, sizeof(title_text),
             "Line profiles through [%d,%d] at 3D position [%d,%d,%d]",
             data_x, data_y, x_coord, y_coord, z_coord);
    Widget title_label = XtVaCreateManagedWidget("title",
        labelWidgetClass, popup_form,
        XtNlabel, title_text,
        XtNwidth, 880,
        NULL);

    Widget x_canvas = XtVaCreateManagedWidget("x_plot",
        simpleWidgetClass, popup_form,
        XtNfromVert, title_label,
        XtNwidth, 880, XtNheight, 180, XtNborderWidth, 1,
        NULL);
    Widget y_canvas = XtVaCreateManagedWidget("y_plot",
        simpleWidgetClass, popup_form,
        XtNfromVert, x_canvas,
        XtNwidth, 880, XtNheight, 180, XtNborderWidth, 1,
        NULL);
    Widget z_canvas = XtVaCreateManagedWidget("z_plot",
        simpleWidgetClass, popup_form,
        XtNfromVert, y_canvas,
        XtNwidth, 880, XtNheight, 180, XtNborderWidth, 1,
        NULL);

    pd->canvases[0] = x_canvas;
    pd->canvases[1] = y_canvas;
    pd->canvases[2] = z_canvas;

    XtAddEventHandler(x_canvas, ExposureMask, False, plot_expose_handler, x_plot_data);
    XtAddEventHandler(y_canvas, ExposureMask, False, plot_expose_handler, y_plot_data);
    XtAddEventHandler(z_canvas, ExposureMask, False, plot_expose_handler, z_plot_data);

    Widget close_button = XtVaCreateManagedWidget("Close",
        commandWidgetClass, popup_form,
        XtNfromVert, z_canvas,
        NULL);
    XtAddCallback(close_button, XtNcallback, close_line_profile_callback, pd);

    /* Layer toggle button next to Close */
    Widget layer_button = XtVaCreateManagedWidget("layer_btn",
        commandWidgetClass, popup_form,
        XtNlabel, "Layer",
        XtNfromVert, z_canvas,
        XtNfromHoriz, close_button,
        NULL);
    XtAddCallback(layer_button, XtNcallback, toggle_line_layer_callback, pd);

    XtPopup(popup_shell, XtGrabNone);
}

/* Popup data for profile (3 plots) */
typedef struct {
    Widget shell;
    Widget mean_canvas, std_canvas, skewness_canvas;
    PlotData *mean_plot;
    PlotData *std_plot;
    PlotData *skewness_plot;
    double *phys_values;    /* physical coordinate values (e.g., height in m) */
    double *layer_values;   /* 1-indexed layer numbers */
    double phys_min, phys_max;
    double layer_min, layer_max;
    int show_layer;         /* 0 = physical coords (default), 1 = layer index */
    char phys_label[32];    /* axis label for physical mode */
    char layer_label[32];   /* axis label for layer-index mode, e.g. "Z Layer" */
} ProfilePopupData;

/* Callback to destroy profile popup and free data */
void close_profile_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    ProfilePopupData *popup_data = (ProfilePopupData *)client_data;

    if (popup_data) {
        if (popup_data->mean_plot) {
            if (popup_data->mean_plot->data) free(popup_data->mean_plot->data);
            /* x_values is owned by popup_data, not individual plots */
            free(popup_data->mean_plot);
        }
        if (popup_data->std_plot) {
            if (popup_data->std_plot->data) free(popup_data->std_plot->data);
            free(popup_data->std_plot);
        }
        if (popup_data->skewness_plot) {
            if (popup_data->skewness_plot->data) free(popup_data->skewness_plot->data);
            free(popup_data->skewness_plot);
        }
        if (popup_data->phys_values) free(popup_data->phys_values);
        if (popup_data->layer_values) free(popup_data->layer_values);
        XtDestroyWidget(popup_data->shell);
        free(popup_data);
    }
}

/* Toggle between physical-coordinate and layer-index display in profile popup */
void toggle_zlayer_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    ProfilePopupData *popup_data = (ProfilePopupData *)client_data;
    popup_data->show_layer = !popup_data->show_layer;

    double *new_x    = popup_data->show_layer ? popup_data->layer_values : popup_data->phys_values;
    double  new_xmin = popup_data->show_layer ? popup_data->layer_min    : popup_data->phys_min;
    double  new_xmax = popup_data->show_layer ? popup_data->layer_max    : popup_data->phys_max;
    const char *new_label = popup_data->show_layer ? popup_data->layer_label : popup_data->phys_label;

    /* Update x axis in all three plots */
    popup_data->mean_plot->x_values = new_x;
    popup_data->mean_plot->xmin = new_xmin;
    popup_data->mean_plot->xmax = new_xmax;
    snprintf(popup_data->mean_plot->xlabel, sizeof(popup_data->mean_plot->xlabel), "%s", new_label);

    popup_data->std_plot->x_values = new_x;
    popup_data->std_plot->xmin = new_xmin;
    popup_data->std_plot->xmax = new_xmax;
    snprintf(popup_data->std_plot->xlabel, sizeof(popup_data->std_plot->xlabel), "%s", new_label);

    popup_data->skewness_plot->x_values = new_x;
    popup_data->skewness_plot->xmin = new_xmin;
    popup_data->skewness_plot->xmax = new_xmax;
    snprintf(popup_data->skewness_plot->xlabel, sizeof(popup_data->skewness_plot->xlabel), "%s", new_label);

    /* Update button label to show what clicking it will switch to */
    XtVaSetValues(w, XtNlabel, popup_data->show_layer ? popup_data->phys_label : popup_data->layer_label, NULL);

    /* Force redraw of all three canvases */
    Display *dpy = XtDisplay(popup_data->mean_canvas);
    if (XtWindow(popup_data->mean_canvas))
        XClearArea(dpy, XtWindow(popup_data->mean_canvas), 0, 0, 0, 0, True);
    if (XtWindow(popup_data->std_canvas))
        XClearArea(dpy, XtWindow(popup_data->std_canvas), 0, 0, 0, 0, True);
    if (XtWindow(popup_data->skewness_canvas))
        XClearArea(dpy, XtWindow(popup_data->skewness_canvas), 0, 0, 0, 0, True);
}

/* Show slice statistics (mean and std) along current axis */
void show_slice_statistics(PlotfileData *pf) {
    const char *axis_names[] = {"X", "Y", "Z"};
    int axis = pf->slice_axis;
    int n_slices = pf->grid_dims[axis];

    /* Determine slice dimensions */
    int slice_dim1, slice_dim2;
    if (axis == 2) {  /* Z slices: X-Y planes */
        slice_dim1 = pf->grid_dims[0];
        slice_dim2 = pf->grid_dims[1];
    } else if (axis == 1) {  /* Y slices: X-Z planes */
        slice_dim1 = pf->grid_dims[0];
        slice_dim2 = pf->grid_dims[2];
    } else {  /* X slices: Y-Z planes */
        slice_dim1 = pf->grid_dims[1];
        slice_dim2 = pf->grid_dims[2];
    }
    int slice_size = slice_dim1 * slice_dim2;

    /* Allocate arrays for mean, std, and skewness */
    double *means = (double *)malloc(n_slices * sizeof(double));
    double *stds = (double *)malloc(n_slices * sizeof(double));
    double *skewness = (double *)malloc(n_slices * sizeof(double));
    double *layer_indices = (double *)malloc(n_slices * sizeof(double));

    /* Physical coordinate values (e.g., height in m for Z axis) */
    double *phys_values = (double *)malloc(n_slices * sizeof(double));
    double dphys = (pf->prob_hi[axis] - pf->prob_lo[axis]) / n_slices;
    for (int s = 0; s < n_slices; s++)
        phys_values[s] = pf->prob_lo[axis] + (s + 0.5) * dphys;
    double phys_min = phys_values[0];
    double phys_max = phys_values[n_slices - 1];

    /* Calculate mean, std, and skewness for each slice */
    for (int s = 0; s < n_slices; s++) {
        layer_indices[s] = s + 1;  /* 1-indexed for display */

        double sum = 0.0;
        double sum_sq = 0.0;

        /* First pass: calculate mean and variance */
        for (int j = 0; j < slice_dim2; j++) {
            for (int i = 0; i < slice_dim1; i++) {
                int idx;
                if (axis == 2) {  /* Z slice */
                    idx = s * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + i;
                } else if (axis == 1) {  /* Y slice */
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + s * pf->grid_dims[0] + i;
                } else {  /* X slice */
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + i * pf->grid_dims[0] + s;
                }
                double val = pf->data[idx];
                sum += val;
                sum_sq += val * val;
            }
        }

        means[s] = sum / slice_size;
        double variance = (sum_sq / slice_size) - (means[s] * means[s]);
        stds[s] = (variance > 0) ? sqrt(variance) : 0.0;

        /* Second pass: calculate skewness (third moment) */
        double sum_third = 0.0;
        for (int j = 0; j < slice_dim2; j++) {
            for (int i = 0; i < slice_dim1; i++) {
                int idx;
                if (axis == 2) {  /* Z slice */
                    idx = s * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + i;
                } else if (axis == 1) {  /* Y slice */
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + s * pf->grid_dims[0] + i;
                } else {  /* X slice */
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + i * pf->grid_dims[0] + s;
                }
                double val = pf->data[idx];
                double diff = val - means[s];
                sum_third += diff * diff * diff;
            }
        }

        /* Skewness = E[(X - mu)^3] / sigma^3 */
        if (stds[s] > 0) {
            double std3 = stds[s] * stds[s] * stds[s];
            skewness[s] = (sum_third / slice_size) / std3;
        } else {
            skewness[s] = 0.0;
        }
    }

    /* Determine physical axis label: "Height (m)" for Z, otherwise axis name with units */
    const char *axis_labels[] = {"X (m)", "Y (m)", "Height (m)"};
    const char *layer_labels[] = {"X Layer", "Y Layer", "Z Layer"};
    const char *phys_label = axis_labels[axis];
    const char *layer_label = layer_labels[axis];

    /* Create plot data for mean - default x axis = physical coordinates */
    PlotData *mean_plot = (PlotData *)malloc(sizeof(PlotData));
    mean_plot->n_points = n_slices;
    mean_plot->data = means;
    mean_plot->x_values = phys_values;
    mean_plot->vmin = 1e30;
    mean_plot->vmax = -1e30;
    for (int i = 0; i < n_slices; i++) {
        if (means[i] < mean_plot->vmin) mean_plot->vmin = means[i];
        if (means[i] > mean_plot->vmax) mean_plot->vmax = means[i];
    }
    mean_plot->xmin = phys_min;
    mean_plot->xmax = phys_max;
    snprintf(mean_plot->title, sizeof(mean_plot->title), "%s Mean along %s axis",
             pf->variables[pf->current_var], axis_names[axis]);
    snprintf(mean_plot->xlabel, sizeof(mean_plot->xlabel), "%s", phys_label);
    snprintf(mean_plot->vlabel, sizeof(mean_plot->vlabel), "%s Mean", pf->variables[pf->current_var]);

    /* Create plot data for std */
    PlotData *std_plot = (PlotData *)malloc(sizeof(PlotData));
    std_plot->n_points = n_slices;
    std_plot->data = stds;
    std_plot->x_values = phys_values;
    std_plot->vmin = 1e30;
    std_plot->vmax = -1e30;
    for (int i = 0; i < n_slices; i++) {
        if (stds[i] < std_plot->vmin) std_plot->vmin = stds[i];
        if (stds[i] > std_plot->vmax) std_plot->vmax = stds[i];
    }
    std_plot->xmin = phys_min;
    std_plot->xmax = phys_max;
    snprintf(std_plot->title, sizeof(std_plot->title), "%s Std Dev along %s axis",
             pf->variables[pf->current_var], axis_names[axis]);
    snprintf(std_plot->xlabel, sizeof(std_plot->xlabel), "%s", phys_label);
    snprintf(std_plot->vlabel, sizeof(std_plot->vlabel), "%s Std", pf->variables[pf->current_var]);

    /* Create plot data for skewness */
    PlotData *skewness_plot = (PlotData *)malloc(sizeof(PlotData));
    skewness_plot->n_points = n_slices;
    skewness_plot->data = skewness;
    skewness_plot->x_values = phys_values;
    skewness_plot->vmin = 1e30;
    skewness_plot->vmax = -1e30;
    for (int i = 0; i < n_slices; i++) {
        if (skewness[i] < skewness_plot->vmin) skewness_plot->vmin = skewness[i];
        if (skewness[i] > skewness_plot->vmax) skewness_plot->vmax = skewness[i];
    }
    skewness_plot->xmin = phys_min;
    skewness_plot->xmax = phys_max;
    snprintf(skewness_plot->title, sizeof(skewness_plot->title), "%s Skewness along %s axis",
             pf->variables[pf->current_var], axis_names[axis]);
    snprintf(skewness_plot->xlabel, sizeof(skewness_plot->xlabel), "%s", phys_label);
    snprintf(skewness_plot->vlabel, sizeof(skewness_plot->vlabel), "%s Skewness", pf->variables[pf->current_var]);

    /* Create popup data structure */
    ProfilePopupData *popup_data = (ProfilePopupData *)malloc(sizeof(ProfilePopupData));
    popup_data->mean_plot = mean_plot;
    popup_data->std_plot = std_plot;
    popup_data->skewness_plot = skewness_plot;
    popup_data->phys_values = phys_values;
    popup_data->layer_values = layer_indices;
    popup_data->phys_min = phys_min;
    popup_data->phys_max = phys_max;
    popup_data->layer_min = 1;
    popup_data->layer_max = n_slices;
    popup_data->show_layer = 0;
    snprintf(popup_data->phys_label, sizeof(popup_data->phys_label), "%s", phys_label);
    snprintf(popup_data->layer_label, sizeof(popup_data->layer_label), "%s", layer_label);

    /* Create popup shell - wider for 3 side-by-side plots */
    Widget popup_shell = XtVaCreatePopupShell("Slice Statistics",
        transientShellWidgetClass, toplevel,
        XtNwidth, 1200,
        XtNheight, 450,
        NULL);

    popup_data->shell = popup_shell;

    Widget popup_form = XtVaCreateManagedWidget("form",
        formWidgetClass, popup_shell,
        NULL);

    /* Mean plot canvas - left */
    Widget mean_canvas = XtVaCreateManagedWidget("mean_plot",
        simpleWidgetClass, popup_form,
        XtNwidth, 380,
        XtNheight, 350,
        XtNborderWidth, 1,
        NULL);

    /* Std plot canvas - middle (next to mean) */
    Widget std_canvas = XtVaCreateManagedWidget("std_plot",
        simpleWidgetClass, popup_form,
        XtNfromHoriz, mean_canvas,
        XtNwidth, 380,
        XtNheight, 350,
        XtNborderWidth, 1,
        NULL);

    /* Kurtosis plot canvas - right (next to std) */
    Widget skewness_canvas = XtVaCreateManagedWidget("skewness_plot",
        simpleWidgetClass, popup_form,
        XtNfromHoriz, std_canvas,
        XtNwidth, 380,
        XtNheight, 350,
        XtNborderWidth, 1,
        NULL);

    /* Store canvas widgets in popup_data for toggle redraws */
    popup_data->mean_canvas = mean_canvas;
    popup_data->std_canvas = std_canvas;
    popup_data->skewness_canvas = skewness_canvas;

    /* Add expose event handlers - using horizontal plot (layer on Y, value on X) */
    XtAddEventHandler(mean_canvas, ExposureMask, False, horizontal_plot_expose_handler, mean_plot);
    XtAddEventHandler(std_canvas, ExposureMask, False, horizontal_plot_expose_handler, std_plot);
    XtAddEventHandler(skewness_canvas, ExposureMask, False, horizontal_plot_expose_handler, skewness_plot);

    /* Close button - below the plots */
    Widget close_button = XtVaCreateManagedWidget("Close",
        commandWidgetClass, popup_form,
        XtNfromVert, mean_canvas,
        NULL);

    XtAddCallback(close_button, XtNcallback, close_profile_popup_callback, popup_data);

    /* Layer toggle button - next to Close; label is axis-dependent */
    Widget zlayer_button = XtVaCreateManagedWidget("layer_btn",
        commandWidgetClass, popup_form,
        XtNlabel, layer_label,
        XtNfromVert, mean_canvas,
        XtNfromHoriz, close_button,
        NULL);

    XtAddCallback(zlayer_button, XtNcallback, toggle_zlayer_callback, popup_data);

    /* Time-Height Contour button (only useful with multiple timesteps) */
    if (n_timesteps > 1) {
        Widget thc_btn = XtVaCreateManagedWidget("Time-Profile Contour",
            commandWidgetClass, popup_form,
            XtNfromVert,  mean_canvas,
            XtNfromHoriz, zlayer_button,
            NULL);
        XtAddCallback(thc_btn, XtNcallback, time_height_contour_callback, (XtPointer)pf);
    }

    /* Show popup */
    XtPopup(popup_shell, XtGrabNone);
}

/* ============================================================
 * Time-Height Contour popup
 * ============================================================ */

typedef struct {
    Widget shell;
    Widget canvas;
    double *contour_data; /* [ntimes * nz] horizontal-mean values */
    double *times;        /* [ntimes] simulation times */
    double *z_values;     /* [nz] z physical coordinates */
    int ntimes, nz;
    double vmin, vmax;
    char title[128];
    char var_name[64];
} TimeHeightContourData;

static void draw_time_height_contour(Display *dpy, Window win, GC gc2,
                                     TimeHeightContourData *thc,
                                     int width, int height) {
    XSetForeground(dpy, gc2, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, gc2, 0, 0, width, height);
    XSetForeground(dpy, gc2, BlackPixel(dpy, screen));
    if (font) XSetFont(dpy, gc2, font->fid);

    if (!thc || thc->ntimes < 1 || thc->nz < 1) {
        XDrawString(dpy, win, gc2, 10, 20, "No data", 7);
        XFlush(dpy); return;
    }

    int plot_left = 75, plot_right = width - 95;
    int plot_top  = 30, plot_bottom = height - 45;
    int plot_w = plot_right - plot_left;
    int plot_h = plot_bottom - plot_top;
    if (plot_w <= 0 || plot_h <= 0) return;

    double tmin = thc->times[0], tmax = thc->times[thc->ntimes - 1];
    if (tmin == tmax) tmax = tmin + 1.0;
    double zmin = thc->z_values[0], zmax = thc->z_values[thc->nz - 1];
    if (zmin == zmax) zmax = zmin + 1.0;
    double vmin = thc->vmin, vmax = thc->vmax;
    if (vmin == vmax) { vmin -= 0.5; vmax += 0.5; }

    /* Render contour image */
    XImage *img = XCreateImage(dpy, DefaultVisual(dpy, screen),
                               DefaultDepth(dpy, screen), ZPixmap, 0,
                               NULL, plot_w, plot_h, 32, 0);
    if (img) {
        img->data = (char *)malloc(img->bytes_per_line * plot_h);
        if (img->data) {
            for (int py = 0; py < plot_h; py++) {
                /* Map pixel row → z (top=zmax, bottom=zmin) */
                double z_frac = 1.0 - (double)py / (plot_h - 1);
                double z_target = zmin + z_frac * (zmax - zmin);
                /* Find nearest z index */
                int zi = 0;
                double best = fabs(thc->z_values[0] - z_target);
                for (int i = 1; i < thc->nz; i++) {
                    double d = fabs(thc->z_values[i] - z_target);
                    if (d < best) { best = d; zi = i; }
                }
                for (int px = 0; px < plot_w; px++) {
                    /* Map pixel col → time index */
                    double t_frac = (double)px / (plot_w - 1);
                    double t_target = tmin + t_frac * (tmax - tmin);
                    int ti = 0;
                    double bdt = fabs(thc->times[0] - t_target);
                    for (int i = 1; i < thc->ntimes; i++) {
                        double dt = fabs(thc->times[i] - t_target);
                        if (dt < bdt) { bdt = dt; ti = i; }
                    }
                    double v = thc->contour_data[ti * thc->nz + zi];
                    double t = (v - vmin) / (vmax - vmin);
                    if (t < 0.0) t = 0.0; if (t > 1.0) t = 1.0;
                    RGB rgb = viridis_colormap(t);
                    unsigned long px_col = ((unsigned long)(rgb.r)<<16)|
                                           ((unsigned long)(rgb.g)<<8)|(rgb.b);
                    XPutPixel(img, px, py, px_col);
                }
            }
            XPutImage(dpy, win, gc2, img, 0, 0, plot_left, plot_top, plot_w, plot_h);
            free(img->data); img->data = NULL;
        }
        XDestroyImage(img);
    }

    /* Axes */
    XSetForeground(dpy, gc2, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, gc2, plot_left, plot_top, plot_w, plot_h);

    char lbl[64];
    /* X (time) ticks */
    for (int i = 0; i <= 5; i++) {
        double t = tmin + (tmax - tmin) * i / 5;
        int xp = plot_left + (int)((double)i / 5 * plot_w);
        XDrawLine(dpy, win, gc2, xp, plot_bottom, xp, plot_bottom + 3);
        profile_fmt_val(lbl, sizeof(lbl), t);
        int lw = font ? XTextWidth(font, lbl, strlen(lbl)) : 40;
        XDrawString(dpy, win, gc2, xp - lw/2, plot_bottom + 16, lbl, strlen(lbl));
    }
    /* Y (z) ticks */
    for (int i = 0; i <= 6; i++) {
        double z = zmin + (zmax - zmin) * i / 6;
        int yp = plot_bottom - (int)((z - zmin) / (zmax - zmin) * plot_h);
        XDrawLine(dpy, win, gc2, plot_left - 3, yp, plot_left, yp);
        profile_fmt_val(lbl, sizeof(lbl), z);
        int lw = font ? XTextWidth(font, lbl, strlen(lbl)) : 40;
        XDrawString(dpy, win, gc2, plot_left - lw - 5, yp + 4, lbl, strlen(lbl));
    }
    /* Axis labels */
    const char *xl = "time (s)";
    int xlw = font ? XTextWidth(font, xl, strlen(xl)) : 50;
    XDrawString(dpy, win, gc2, plot_left + (plot_w - xlw)/2, plot_bottom + 32, xl, strlen(xl));
    XDrawString(dpy, win, gc2, 2, plot_top + 10, "z (m)", 5);
    /* Title */
    XDrawString(dpy, win, gc2, plot_left + 4, plot_top - 2, thc->title, strlen(thc->title));

    /* Colorbar */
    int cb_x = plot_right + 8, cb_w = 14;
    for (int py = 0; py < plot_h; py++) {
        double t = 1.0 - (double)py / (plot_h - 1);
        RGB rgb = viridis_colormap(t);
        unsigned long px_col = ((unsigned long)(rgb.r)<<16)|((unsigned long)(rgb.g)<<8)|(rgb.b);
        XSetForeground(dpy, gc2, px_col);
        XFillRectangle(dpy, win, gc2, cb_x, plot_top + py, cb_w, 1);
    }
    XSetForeground(dpy, gc2, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, gc2, cb_x, plot_top, cb_w, plot_h);
    int lx = cb_x + cb_w + 4;
    profile_fmt_val(lbl, sizeof(lbl), vmax);
    XDrawLine(dpy, win, gc2, cb_x, plot_top, cb_x + cb_w + 3, plot_top);
    XDrawString(dpy, win, gc2, lx, plot_top + 10, lbl, strlen(lbl));
    profile_fmt_val(lbl, sizeof(lbl), (vmin + vmax) * 0.5);
    int mid_y = plot_top + plot_h / 2;
    XDrawLine(dpy, win, gc2, cb_x, mid_y, cb_x + cb_w + 3, mid_y);
    XDrawString(dpy, win, gc2, lx, mid_y + 4, lbl, strlen(lbl));
    profile_fmt_val(lbl, sizeof(lbl), vmin);
    XDrawLine(dpy, win, gc2, cb_x, plot_bottom, cb_x + cb_w + 3, plot_bottom);
    XDrawString(dpy, win, gc2, lx, plot_bottom + 4, lbl, strlen(lbl));

    XFlush(dpy);
}

static void time_height_expose_handler(Widget w, XtPointer client_data,
                                       XEvent *event, Boolean *cont) {
    (void)event; (void)cont;
    TimeHeightContourData *thc = (TimeHeightContourData *)client_data;
    Window win = XtWindow(w);
    unsigned int wd, ht, bw, depth;
    Window root; int x, y;
    XGetGeometry(display, win, &root, &x, &y, &wd, &ht, &bw, &depth);
    GC gc2 = XCreateGC(display, win, 0, NULL);
    if (font) XSetFont(display, gc2, font->fid);
    draw_time_height_contour(display, win, gc2, thc, (int)wd, (int)ht);
    XFreeGC(display, gc2);
}

static void time_height_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)call_data;
    TimeHeightContourData *thc = (TimeHeightContourData *)client_data;
    if (thc) {
        XtDestroyWidget(thc->shell);
        if (thc->contour_data) free(thc->contour_data);
        if (thc->times)        free(thc->times);
        if (thc->z_values)     free(thc->z_values);
        free(thc);
    }
}

void show_time_height_contour(PlotfileData *pf) {
    if (n_timesteps < 2) {
        /* nothing to show for a single timestep */
        fprintf(stderr, "Time-Height Contour requires multiple timesteps.\n");
        return;
    }

    int var_idx  = pf->current_var;
    int saved_ts = current_timestep;

    /* We use axis 2 (z) as the height axis always */
    int nz = pf->grid_dims[2];
    (void)pf->grid_dims[0]; (void)pf->grid_dims[1]; /* nx/ny read per-timestep below */

    TimeHeightContourData *thc = (TimeHeightContourData *)calloc(1, sizeof(TimeHeightContourData));
    thc->ntimes       = n_timesteps;
    thc->nz           = nz;
    thc->contour_data = (double *)malloc(n_timesteps * nz * sizeof(double));
    thc->times        = (double *)malloc(n_timesteps * sizeof(double));
    thc->z_values     = (double *)malloc(nz * sizeof(double));

    /* Fill z physical coords from current plotfile */
    double dz = (pf->prob_hi[2] - pf->prob_lo[2]) / nz;
    for (int k = 0; k < nz; k++)
        thc->z_values[k] = pf->prob_lo[2] + (k + 0.5) * dz;

    thc->vmin =  1e300;
    thc->vmax = -1e300;

    /* Iterate over all timesteps, compute horizontal mean per z-level */
    for (int ti = 0; ti < n_timesteps; ti++) {
        /* Load this timestep's data without triggering GUI updates */
        PlotfileData tmp_pf;
        memset(&tmp_pf, 0, sizeof(tmp_pf));
        strncpy(tmp_pf.plotfile_dir, timestep_paths[ti], MAX_PATH - 1);
        tmp_pf.current_var  = var_idx;
        tmp_pf.slice_axis   = pf->slice_axis;
        tmp_pf.slice_idx    = pf->slice_idx;
        tmp_pf.colormap     = pf->colormap;
        tmp_pf.current_level = 0;

        if (read_header(&tmp_pf) < 0)        { thc->times[ti] = ti; continue; }
        if (read_cell_h(&tmp_pf) < 0)        { thc->times[ti] = tmp_pf.time; continue; }
        if (read_variable_data(&tmp_pf, var_idx) < 0) { thc->times[ti] = tmp_pf.time; continue; }

        thc->times[ti] = tmp_pf.time;

        int tnx = tmp_pf.grid_dims[0];
        int tny = tmp_pf.grid_dims[1];
        int tnz = tmp_pf.grid_dims[2];
        int use_nz = (tnz < nz) ? tnz : nz;

        for (int k = 0; k < use_nz; k++) {
            double sum = 0.0;
            long count = 0;
            for (int j = 0; j < tny; j++) {
                for (int i = 0; i < tnx; i++) {
                    long idx = (long)k * tnx * tny + (long)j * tnx + i;
                    if (idx < (long)tnx * tny * tnz) {
                        sum += tmp_pf.data[idx];
                        count++;
                    }
                }
            }
            double mean_val = (count > 0) ? sum / count : 0.0;
            thc->contour_data[ti * nz + k] = mean_val;
            if (mean_val < thc->vmin) thc->vmin = mean_val;
            if (mean_val > thc->vmax) thc->vmax = mean_val;
        }
        /* Zero-fill any extra z levels if tnz < nz */
        for (int k = use_nz; k < nz; k++)
            thc->contour_data[ti * nz + k] = 0.0;

        if (tmp_pf.data) { free(tmp_pf.data); tmp_pf.data = NULL; }
    }

    /* Restore original state */
    strncpy(pf->plotfile_dir, timestep_paths[saved_ts], MAX_PATH - 1);

    snprintf(thc->title, sizeof(thc->title),
             "Time-Height Contour: %s (horiz. mean)", pf->variables[var_idx]);
    strncpy(thc->var_name, pf->variables[var_idx], sizeof(thc->var_name) - 1);

    /* Create popup */
    Widget shell = XtVaCreatePopupShell("Time-Height Contour",
        transientShellWidgetClass, toplevel,
        XtNwidth,  860,
        XtNheight, 520,
        NULL);
    thc->shell = shell;

    Widget frm = XtVaCreateManagedWidget("form", formWidgetClass, shell, NULL);

    Widget cnv = XtVaCreateManagedWidget("thContour",
        simpleWidgetClass, frm,
        XtNwidth,  840, XtNheight, 460,
        XtNborderWidth, 1,
        XtNtop,    XawChainTop,
        XtNbottom, XawChainBottom,
        XtNleft,   XawChainLeft,
        XtNright,  XawChainRight,
        NULL);
    thc->canvas = cnv;
    XtAddEventHandler(cnv, ExposureMask, False, time_height_expose_handler, thc);

    Widget close_btn = XtVaCreateManagedWidget("Close",
        commandWidgetClass, frm,
        XtNfromVert, cnv,
        NULL);
    XtAddCallback(close_btn, XtNcallback, time_height_close_callback, thc);

    XtPopup(shell, XtGrabNone);
}

static void time_height_contour_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)call_data;
    PlotfileData *pf = (PlotfileData *)client_data;
    if (pf && pf->data) show_time_height_contour(pf);
}

/* Profile button callback */
void profile_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data) {
        show_slice_statistics(global_pf);
    }
}

/* ============================================================
 * 2D FFT Energy Spectrum
 * ============================================================ */

typedef struct { double r, i; } Fft2DComplex;

/* In-place radix-2 Cooley-Tukey FFT; N must be a power of 2, forward transform */
static void fft1d_ct(Fft2DComplex *x, int N) {
    /* Bit-reversal permutation */
    for (int i = 1, j = 0; i < N; i++) {
        int bit = N >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { Fft2DComplex tmp = x[i]; x[i] = x[j]; x[j] = tmp; }
    }
    /* Butterfly stages */
    for (int len = 2; len <= N; len <<= 1) {
        double ang = -2.0 * M_PI / len;
        double wr0 = cos(ang), wi0 = sin(ang);
        for (int i = 0; i < N; i += len) {
            double wr = 1.0, wi = 0.0;
            for (int k = 0; k < len / 2; k++) {
                Fft2DComplex u = x[i + k];
                Fft2DComplex v = { x[i+k+len/2].r*wr - x[i+k+len/2].i*wi,
                                   x[i+k+len/2].r*wi + x[i+k+len/2].i*wr };
                x[i+k].r = u.r + v.r;         x[i+k].i = u.i + v.i;
                x[i+k+len/2].r = u.r - v.r;   x[i+k+len/2].i = u.i - v.i;
                double tmp_wr = wr*wr0 - wi*wi0;
                wi = wr*wi0 + wi*wr0;  wr = tmp_wr;
            }
        }
    }
}

typedef struct {
    Widget shell;
    Widget canvas;
    Widget method_fft_btn;   /* 2D FFT method button */
    Widget method_wk_btn;    /* Wiener-Khinchin method button */
    PlotfileData *pf;
    double *k_vals;   /* Physical wavenumber for each bin (rad/unit) */
    double *e_vals;   /* Energy per bin: sum of |F|^2 / (Nx*Ny)^2 */
    double *e_vals_y; /* For WK: vertical spectrum (NULL for 2D FFT) */
    int n_bins;       /* Number of valid bins */
    int n_bins_y;     /* For WK: number of bins in vertical spectrum */
    int fft_nx;       /* FFT grid size (power-of-2) used in X */
    int fft_ny;       /* FFT grid size (power-of-2) used in Y */
    int method;       /* 0 = 2D FFT, 1 = Wiener-Khinchin */
} FFTPopupData;

static FFTPopupData *g_fft_popup = NULL;

static void close_fft2d_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    FFTPopupData *popup = (FFTPopupData *)client_data;
    if (popup) {
        if (popup->k_vals) free(popup->k_vals);
        if (popup->e_vals) free(popup->e_vals);
        if (popup->e_vals_y) free(popup->e_vals_y);
        XtDestroyWidget(popup->shell);
        free(popup);
        if (g_fft_popup == popup) g_fft_popup = NULL;
    }
}

static void compute_2dfft_spectrum(FFTPopupData *popup) {
    PlotfileData *pf = popup->pf;
    if (!current_slice_data || slice_width <= 0 || slice_height <= 0) return;

    /* Square region: min of dimensions, then largest power-of-2 (cap at 2048) */
    int S = (slice_width < slice_height) ? slice_width : slice_height;
    int N = 1;
    while (N * 2 <= S && N < 2048) N <<= 1;
    int Nx = N, Ny = N;
    if (N < 4) return;

    /* Center offsets for extracting the square region */
    int ox = (slice_width - N) / 2;
    int oy = (slice_height - N) / 2;

    Fft2DComplex *buf = (Fft2DComplex *)malloc((size_t)Nx * Ny * sizeof(Fft2DComplex));
    if (!buf) return;

    /* Fill buf with Hann-windowed slice data (subtract mean to reduce DC spike) */
    double mean_val = 0.0;
    for (int iy = 0; iy < N; iy++)
        for (int ix = 0; ix < N; ix++)
            mean_val += current_slice_data[(oy + iy) * slice_width + (ox + ix)];
    mean_val /= (double)(N * N);

    for (int iy = 0; iy < N; iy++) {
        double wy = 0.5 * (1.0 - cos(2.0 * M_PI * iy / (N - 1)));
        for (int ix = 0; ix < N; ix++) {
            double wx = 0.5 * (1.0 - cos(2.0 * M_PI * ix / (N - 1)));
            buf[iy * N + ix].r = (current_slice_data[(oy + iy) * slice_width + (ox + ix)] - mean_val) * wx * wy;
            buf[iy * N + ix].i = 0.0;
        }
    }

    /* 2D FFT: row-wise, then column-wise */
    int tmp_len = (Nx > Ny) ? Nx : Ny;
    Fft2DComplex *tmp = (Fft2DComplex *)malloc(tmp_len * sizeof(Fft2DComplex));
    if (!tmp) { free(buf); return; }

    for (int iy = 0; iy < Ny; iy++) {
        memcpy(tmp, &buf[iy * Nx], Nx * sizeof(Fft2DComplex));
        fft1d_ct(tmp, Nx);
        memcpy(&buf[iy * Nx], tmp, Nx * sizeof(Fft2DComplex));
    }
    for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) tmp[iy] = buf[iy * Nx + ix];
        fft1d_ct(tmp, Ny);
        for (int iy = 0; iy < Ny; iy++) buf[iy * Nx + ix] = tmp[iy];
    }
    free(tmp);

    /* Determine physical domain sizes for the slice axes */
    int axis = pf->slice_axis;
    int ax1 = (axis == 0) ? 1 : 0;   /* first  slice axis index in prob_lo/hi */
    int ax2 = (axis == 1) ? 2 : (axis == 0 ? 2 : 1); /* second slice axis */
    double Lx_fft = (pf->prob_hi[ax1] - pf->prob_lo[ax1]) * Nx / pf->grid_dims[ax1];
    double Ly_fft = (pf->prob_hi[ax2] - pf->prob_lo[ax2]) * Ny / pf->grid_dims[ax2];

    /* Minimum dk for radial binning (rad per physical unit) */
    double dk_min = 2.0 * M_PI / ((Lx_fft > Ly_fft) ? Lx_fft : Ly_fft);

    /* Max k we can reliably represent (inscribed circle of the first Brillouin zone) */
    double kx_nyq = M_PI * Nx / Lx_fft;
    double ky_nyq = M_PI * Ny / Ly_fft;
    double k_max_reliable = (kx_nyq < ky_nyq) ? kx_nyq : ky_nyq;
    int n_bins = (int)(k_max_reliable / dk_min) + 2;
    if (n_bins > 4096) n_bins = 4096;

    double *e_sum = (double *)calloc(n_bins, sizeof(double));
    if (!e_sum) { free(buf); return; }

    double norm2 = 1.0 / ((double)Nx * (double)Ny);
    norm2 *= norm2;

    for (int iy = 0; iy < Ny; iy++) {
        int ky_idx = (iy <= Ny / 2) ? iy : iy - Ny;
        double ky_phys = 2.0 * M_PI * ky_idx / Ly_fft;
        for (int ix = 0; ix < Nx; ix++) {
            int kx_idx = (ix <= Nx / 2) ? ix : ix - Nx;
            double kx_phys = 2.0 * M_PI * kx_idx / Lx_fft;
            double k_phys = sqrt(kx_phys * kx_phys + ky_phys * ky_phys);
            int bin = (int)(k_phys / dk_min);
            if (bin > 0 && bin < n_bins) {
                Fft2DComplex f = buf[iy * Nx + ix];
                e_sum[bin] += (f.r * f.r + f.i * f.i) * norm2;
            }
        }
    }
    free(buf);

    /* Record FFT grid size for display */
    popup->fft_nx = Nx;
    popup->fft_ny = Ny;

    /* Build output k and E arrays (skip empty bins) */
    if (popup->k_vals) free(popup->k_vals);
    if (popup->e_vals) free(popup->e_vals);
    popup->k_vals = (double *)malloc(n_bins * sizeof(double));
    popup->e_vals = (double *)malloc(n_bins * sizeof(double));
    popup->n_bins = 0;

    for (int b = 1; b < n_bins; b++) {
        if (e_sum[b] > 0.0) {
            popup->k_vals[popup->n_bins] = b * dk_min;
            popup->e_vals[popup->n_bins] = e_sum[b];
            popup->n_bins++;
        }
    }
    free(e_sum);
}

/* Wiener-Khinchin method: row/column periodograms (Ma et al. 2024 approach) */
static void compute_wk_spectrum(FFTPopupData *popup) {
    PlotfileData *pf = popup->pf;
    if (!current_slice_data || slice_width <= 0 || slice_height <= 0) return;

    int W = slice_width, H = slice_height;
    popup->fft_nx = W;
    popup->fft_ny = H;

    /* Global mean subtraction */
    double mean_val = 0.0;
    for (int i = 0; i < W * H; i++)
        mean_val += current_slice_data[i];
    mean_val /= (double)(W * H);

    /* Determine physical domain size for the primary (horizontal) slice axis */
    int axis = pf->slice_axis;
    int ax1 = (axis == 0) ? 1 : 0;
    double Lx = pf->prob_hi[ax1] - pf->prob_lo[ax1];
    /* Note: For WK, we use the same k scale (based on Lx) for both spectra,
       matching spectracam's approach of a single shared k axis */

    /* E_x: average |FFT_row|² over all rows */
    int kMaxX = W / 2;
    double *Ex = (double *)calloc(kMaxX + 1, sizeof(double));
    if (!Ex) return;

    /* Allocate row buffers for FFT */
    Fft2DComplex *row_buf = (Fft2DComplex *)malloc(W * sizeof(Fft2DComplex));
    if (!row_buf) { free(Ex); return; }

    /* Find next power of 2 >= W for zero-padded FFT */
    int W_padded = 1;
    while (W_padded < W) W_padded <<= 1;
    Fft2DComplex *row_padded = (Fft2DComplex *)malloc(W_padded * sizeof(Fft2DComplex));
    if (!row_padded) { free(Ex); free(row_buf); return; }

    for (int y = 0; y < H; y++) {
        /* Fill buffer with mean-subtracted row data, zero-pad */
        for (int x = 0; x < W; x++) {
            row_padded[x].r = current_slice_data[y * W + x] - mean_val;
            row_padded[x].i = 0.0;
        }
        for (int x = W; x < W_padded; x++) {
            row_padded[x].r = 0.0;
            row_padded[x].i = 0.0;
        }
        fft1d_ct(row_padded, W_padded);
        /* Accumulate power for first kMaxX bins (scaled to original W) */
        for (int k = 1; k <= kMaxX; k++) {
            int idx = (k * W_padded) / W;  /* Scale bin index */
            if (idx < W_padded)
                Ex[k] += row_padded[idx].r * row_padded[idx].r + row_padded[idx].i * row_padded[idx].i;
        }
    }
    for (int k = 1; k <= kMaxX; k++) Ex[k] /= H;
    free(row_padded);
    free(row_buf);

    /* E_y: average |FFT_col|² over all columns */
    int kMaxY = H / 2;
    double *Ey = (double *)calloc(kMaxY + 1, sizeof(double));
    if (!Ey) { free(Ex); return; }

    int H_padded = 1;
    while (H_padded < H) H_padded <<= 1;
    Fft2DComplex *col_padded = (Fft2DComplex *)malloc(H_padded * sizeof(Fft2DComplex));
    if (!col_padded) { free(Ex); free(Ey); return; }

    for (int x = 0; x < W; x++) {
        for (int y = 0; y < H; y++) {
            col_padded[y].r = current_slice_data[y * W + x] - mean_val;
            col_padded[y].i = 0.0;
        }
        for (int y = H; y < H_padded; y++) {
            col_padded[y].r = 0.0;
            col_padded[y].i = 0.0;
        }
        fft1d_ct(col_padded, H_padded);
        for (int k = 1; k <= kMaxY; k++) {
            int idx = (k * H_padded) / H;
            if (idx < H_padded)
                Ey[k] += col_padded[idx].r * col_padded[idx].r + col_padded[idx].i * col_padded[idx].i;
        }
    }
    for (int k = 1; k <= kMaxY; k++) Ey[k] /= W;
    free(col_padded);

    /* Store results with physical wavenumber */
    if (popup->k_vals) free(popup->k_vals);
    if (popup->e_vals) free(popup->e_vals);
    if (popup->e_vals_y) free(popup->e_vals_y);

    popup->k_vals = (double *)malloc((kMaxX + 1) * sizeof(double));
    popup->e_vals = (double *)malloc((kMaxX + 1) * sizeof(double));
    popup->e_vals_y = (double *)malloc((kMaxY + 1) * sizeof(double));
    popup->n_bins = 0;
    popup->n_bins_y = 0;

    /* Physical wavenumber: k = 2π * index / L */
    for (int k = 1; k <= kMaxX; k++) {
        if (Ex[k] > 0.0) {
            popup->k_vals[popup->n_bins] = 2.0 * M_PI * k / Lx;
            popup->e_vals[popup->n_bins] = Ex[k];
            popup->n_bins++;
        }
    }
    for (int k = 1; k <= kMaxY; k++) {
        if (Ey[k] > 0.0) {
            popup->e_vals_y[popup->n_bins_y] = Ey[k];
            popup->n_bins_y++;
        }
    }

    free(Ex);
    free(Ey);
}

/* Recompute spectrum based on current method */
static void recompute_fft_spectrum(FFTPopupData *popup) {
    if (popup->method == 0) {
        if (popup->e_vals_y) { free(popup->e_vals_y); popup->e_vals_y = NULL; }
        popup->n_bins_y = 0;
        compute_2dfft_spectrum(popup);
    } else {
        compute_wk_spectrum(popup);
    }
}

/* Forward declaration for button callback */
static void draw_fft_spectrum(FFTPopupData *popup);

/* Method button callbacks */
static void fft_method_fft2d_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    FFTPopupData *popup = (FFTPopupData *)client_data;
    if (popup->method == 0) return;  /* Already selected */
    popup->method = 0;
    /* Update button labels to show selection */
    XtVaSetValues(popup->method_fft_btn, XtNlabel, "[2D FFT]", NULL);
    XtVaSetValues(popup->method_wk_btn, XtNlabel, "Wiener-Khinchin", NULL);
    recompute_fft_spectrum(popup);
    draw_fft_spectrum(popup);
}

static void fft_method_wk_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    FFTPopupData *popup = (FFTPopupData *)client_data;
    if (popup->method == 1) return;  /* Already selected */
    popup->method = 1;
    XtVaSetValues(popup->method_fft_btn, XtNlabel, "2D FFT", NULL);
    XtVaSetValues(popup->method_wk_btn, XtNlabel, "[Wiener-Khinchin]", NULL);
    recompute_fft_spectrum(popup);
    draw_fft_spectrum(popup);
}

static void draw_fft_spectrum(FFTPopupData *popup) {
    if (!popup || !popup->canvas || !XtIsRealized(popup->canvas)) return;
    if (popup->n_bins < 2) return;

    Window win = XtWindow(popup->canvas);
    Dimension cw, ch;
    XtVaGetValues(popup->canvas, XtNwidth, &cw, XtNheight, &ch, NULL);
    int width = (int)cw, height = (int)ch;

    GC gc2 = XCreateGC(display, win, 0, NULL);
    if (font) XSetFont(display, gc2, font->fid);

    /* Background and border */
    XSetForeground(display, gc2, WhitePixel(display, screen));
    XFillRectangle(display, win, gc2, 0, 0, width, height);
    XSetForeground(display, gc2, BlackPixel(display, screen));
    XDrawRectangle(display, win, gc2, 0, 0, width - 1, height - 1);

    /* Title */
    PlotfileData *pf = popup->pf;
    char title[256];
    snprintf(title, sizeof(title), "Energy Spectrum: %s  (log-log)",
             (pf && pf->current_var >= 0) ? pf->variables[pf->current_var] : "");
    XDrawString(display, win, gc2, 10, 18, title, strlen(title));

    /* Plot area margins (leave 130 px below pb for axis labels + method notes) */
    int pl = 80, pr = width - 20;
    int pt = 30, pb = height - 130;
    int pw = pr - pl, ph = pb - pt;
    if (pw <= 10 || ph <= 10) { XFreeGC(display, gc2); return; }

    /* Compute log10 range of k and E */
    double lkmin = 1e30, lkmax = -1e30, lemin = 1e30, lemax = -1e30;
    for (int b = 0; b < popup->n_bins; b++) {
        if (popup->k_vals[b] > 0.0 && popup->e_vals[b] > 0.0) {
            double lk = log10(popup->k_vals[b]);
            double le = log10(popup->e_vals[b]);
            if (lk < lkmin) lkmin = lk;
            if (lk > lkmax) lkmax = lk;
            if (le < lemin) lemin = le;
            if (le > lemax) lemax = le;
        }
    }
    /* Include WK vertical spectrum in range calculation */
    if (popup->method == 1 && popup->e_vals_y && popup->n_bins_y > 0) {
        int n_y = (popup->n_bins_y < popup->n_bins) ? popup->n_bins_y : popup->n_bins;
        for (int b = 0; b < n_y; b++) {
            if (popup->e_vals_y[b] > 0.0) {
                double le = log10(popup->e_vals_y[b]);
                if (le < lemin) lemin = le;
                if (le > lemax) lemax = le;
            }
        }
    }
    if (lkmax <= lkmin || lemax <= lemin) { XFreeGC(display, gc2); return; }

    double kpad = (lkmax - lkmin) * 0.05;
    double epad = (lemax - lemin) * 0.05;
    lkmin -= kpad;  lkmax += kpad;
    lemin -= epad;  lemax += epad;

    /* Axes */
    XDrawLine(display, win, gc2, pl, pb, pr, pb);
    XDrawLine(display, win, gc2, pl, pt, pl, pb);

#define FFT_K2X(lk) (pl + (int)(((lk) - lkmin) / (lkmax - lkmin) * pw))
#define FFT_E2Y(le) (pb - (int)(((le) - lemin) / (lemax - lemin) * ph))

    /* X-axis ticks */
    char lbl[64];
    for (int i = 0; i <= 5; i++) {
        double lk = lkmin + (lkmax - lkmin) * i / 5.0;
        int xp = FFT_K2X(lk);
        XDrawLine(display, win, gc2, xp, pb, xp, pb + 4);
        snprintf(lbl, sizeof(lbl), "%.1e", pow(10.0, lk));
        if (font) {
            int lw = XTextWidth(font, lbl, strlen(lbl));
            XDrawString(display, win, gc2, xp - lw / 2, pb + 16, lbl, strlen(lbl));
        }
    }
    /* Y-axis ticks */
    for (int i = 0; i <= 5; i++) {
        double le = lemin + (lemax - lemin) * i / 5.0;
        int yp = FFT_E2Y(le);
        XDrawLine(display, win, gc2, pl - 4, yp, pl, yp);
        snprintf(lbl, sizeof(lbl), "%.1e", pow(10.0, le));
        if (font) {
            int lw = XTextWidth(font, lbl, strlen(lbl));
            XDrawString(display, win, gc2, pl - lw - 6, yp + 4, lbl, strlen(lbl));
        }
    }
    /* Axis labels */
    if (font) {
        const char *xl = "k  (rad / length unit)";
        int lw = XTextWidth(font, xl, strlen(xl));
        XDrawString(display, win, gc2, pl + pw / 2 - lw / 2, pb + 34, xl, strlen(xl));
        const char *yl = "E(k)";
        XDrawString(display, win, gc2, 2, pt + ph / 2, yl, strlen(yl));
    }

    /* Plot primary spectrum (cyan for both methods) */
    XColor cyan_c, dummy_c;
    unsigned long cyan_pix = BlackPixel(display, screen);
    if (XAllocNamedColor(display, DefaultColormap(display, screen), "#00c8e0", &cyan_c, &dummy_c))
        cyan_pix = cyan_c.pixel;
    XSetForeground(display, gc2, cyan_pix);
    XSetLineAttributes(display, gc2, 2, LineSolid, CapButt, JoinMiter);

    int px = -1, py = -1;
    for (int b = 0; b < popup->n_bins; b++) {
        if (popup->k_vals[b] > 0.0 && popup->e_vals[b] > 0.0) {
            int xp = FFT_K2X(log10(popup->k_vals[b]));
            int yp = FFT_E2Y(log10(popup->e_vals[b]));
            if (xp >= pl && xp <= pr && yp >= pt && yp <= pb) {
                if (px >= 0) XDrawLine(display, win, gc2, px, py, xp, yp);
                px = xp;  py = yp;
            } else { px = -1;  py = -1; }
        }
    }

    /* For Wiener-Khinchin, also plot vertical spectrum (green) */
    if (popup->method == 1 && popup->e_vals_y && popup->n_bins_y > 0) {
        XColor green_c;
        unsigned long green_pix = BlackPixel(display, screen);
        if (XAllocNamedColor(display, DefaultColormap(display, screen), "#4ade80", &green_c, &dummy_c))
            green_pix = green_c.pixel;
        XSetForeground(display, gc2, green_pix);

        px = -1; py = -1;
        /* For WK vertical, use same k_vals array (they share wavenumber scale) */
        int n_y = (popup->n_bins_y < popup->n_bins) ? popup->n_bins_y : popup->n_bins;
        for (int b = 0; b < n_y; b++) {
            if (popup->k_vals[b] > 0.0 && popup->e_vals_y[b] > 0.0) {
                int xp = FFT_K2X(log10(popup->k_vals[b]));
                int yp = FFT_E2Y(log10(popup->e_vals_y[b]));
                if (xp >= pl && xp <= pr && yp >= pt && yp <= pb) {
                    if (px >= 0) XDrawLine(display, win, gc2, px, py, xp, yp);
                    px = xp;  py = yp;
                } else { px = -1;  py = -1; }
            }
        }
    }

    /* Reference slope lines anchored at ~1/5 of the spectrum (left-middle region) */
    int anchor_idx = popup->n_bins / 5;  /* 1/5 of the way through */
    if (anchor_idx < 1) anchor_idx = 1;
    if (anchor_idx >= popup->n_bins) anchor_idx = popup->n_bins - 1;
    double lk_anchor = log10(popup->k_vals[anchor_idx]);
    double le_anchor = log10(popup->e_vals[anchor_idx]);
    /* Ensure valid anchor values */
    if (popup->k_vals[anchor_idx] <= 0.0 || popup->e_vals[anchor_idx] <= 0.0) {
        /* Fallback: search for first valid point near 1/3 */
        lk_anchor = (lkmin + lkmax) * 0.5;
        le_anchor = (lemin + lemax) * 0.5;
        for (int b = anchor_idx; b < popup->n_bins && b < anchor_idx + 10; b++) {
            if (popup->k_vals[b] > 0.0 && popup->e_vals[b] > 0.0) {
                lk_anchor = log10(popup->k_vals[b]);
                le_anchor = log10(popup->e_vals[b]);
                break;
            }
        }
    }

    /* Helper to draw a reference line */
    #define DRAW_REF_LINE(slope, r, g, b, dash1, dash2, label_str, label_offset) do { \
        double C = le_anchor - (slope) * lk_anchor; \
        double lk1 = lkmin, lk2 = lkmax; \
        double le1 = C + (slope) * lk1; \
        double le2 = C + (slope) * lk2; \
        if (le1 > lemax) { lk1 = (lemax - C) / (slope); le1 = lemax; } \
        else if (le1 < lemin) { lk1 = (lemin - C) / (slope); le1 = lemin; } \
        if (le2 > lemax) { lk2 = (lemax - C) / (slope); le2 = lemax; } \
        else if (le2 < lemin) { lk2 = (lemin - C) / (slope); le2 = lemin; } \
        if (lk1 < lk2) { \
            XColor rc; unsigned long rpix = BlackPixel(display, screen); \
            rc.red = (r) << 8; rc.green = (g) << 8; rc.blue = (b) << 8; rc.flags = DoRed | DoGreen | DoBlue; \
            if (XAllocColor(display, DefaultColormap(display, screen), &rc)) rpix = rc.pixel; \
            XSetForeground(display, gc2, rpix); \
            char dash_list[] = {dash1, dash2}; \
            XSetLineAttributes(display, gc2, 2, LineOnOffDash, CapButt, JoinMiter); \
            XSetDashes(display, gc2, 0, dash_list, 2); \
            XDrawLine(display, win, gc2, FFT_K2X(lk1), FFT_E2Y(le1), FFT_K2X(lk2), FFT_E2Y(le2)); \
            if (font) { \
                double lk_lbl = lk2 - (lk2 - lk1) * 0.15; /* 85% along the line */ \
                double le_lbl = C + (slope) * lk_lbl; \
                int kx = FFT_K2X(lk_lbl) + 4; /* Right of this point */ \
                int ky = FFT_E2Y(le_lbl) + 2 + (label_offset); /* Above the line */ \
                if (kx > pl + 20 && kx < pr - 5 && ky > pt + 5 && ky < pb - 5) \
                    XDrawString(display, win, gc2, kx, ky, label_str, strlen(label_str)); \
            } \
        } \
    } while(0)

    /* k^(-1): Purple #a78bfa (167, 139, 250) - Batchelor */
    DRAW_REF_LINE(-1.0, 167, 139, 250, 4, 4, "k^(-1)", -10);

    /* k^(-5/3): Orange #f0a500 (240, 165, 0) - Kolmogorov */
    DRAW_REF_LINE(-5.0/3.0, 240, 165, 0, 8, 5, "k^(-5/3)", -10);

    /* k^(-3): Hot pink #ff69b4 (255, 105, 180) - Kraichnan */
    DRAW_REF_LINE(-3.0, 255, 105, 180, 6, 4, "k^(-3)", -10);

    #undef DRAW_REF_LINE

#undef FFT_K2X
#undef FFT_E2Y

    /* Method notes drawn below the x-axis label */
    if (font) {
        XSetForeground(display, gc2, BlackPixel(display, screen));
        XSetLineAttributes(display, gc2, 1, LineSolid, CapButt, JoinMiter);

        /* Separator line */
        int sep_y = pb + 50;
        XDrawLine(display, win, gc2, pl, sep_y, pr, sep_y);

        char note[512];
        int ny = sep_y + 14;
        const int line_h = 14;
        int max_text_width = pr - pl - 10;

        /* Helper macro to draw wrapped text */
        #define DRAW_WRAPPED(text) do { \
            const char *src = (text); \
            int src_len = strlen(src); \
            int start = 0; \
            while (start < src_len) { \
                int end = src_len; \
                while (end > start && XTextWidth(font, src + start, end - start) > max_text_width) { \
                    /* Find last space before max width */ \
                    int last_space = end; \
                    for (int i = end - 1; i > start; i--) { \
                        if (src[i] == ' ') { last_space = i; break; } \
                    } \
                    end = (last_space > start) ? last_space : end - 1; \
                } \
                XDrawString(display, win, gc2, pl, ny, src + start, end - start); \
                ny += line_h; \
                start = end; \
                while (start < src_len && src[start] == ' ') start++; \
            } \
        } while(0)

        if (popup->method == 0) {
            snprintf(note, sizeof(note),
                     "Method: 2D Cooley-Tukey FFT on '%s', "
                     "centered square crop to %d x %d (power-of-2), "
                     "Hann windowed, mean subtracted.",
                     (pf && pf->current_var >= 0) ? pf->variables[pf->current_var] : "field",
                     popup->fft_nx, popup->fft_ny);
            DRAW_WRAPPED(note);

            snprintf(note, sizeof(note),
                     "E(k) = sum|F(kx,ky)|^2 / (Nx*Ny)^2 "
                     "summed over annular shells of physical wavenumber k.");
            DRAW_WRAPPED(note);
        } else {
            snprintf(note, sizeof(note),
                     "Method: Wiener-Khinchin (Ma et al. 2024) on '%s', "
                     "full %d x %d slice, mean subtracted.",
                     (pf && pf->current_var >= 0) ? pf->variables[pf->current_var] : "field",
                     popup->fft_nx, popup->fft_ny);
            DRAW_WRAPPED(note);

            snprintf(note, sizeof(note),
                     "Cyan: E_x(k) horizontal avg, Green: E_y(k) vertical avg. "
                     "Overlap indicates isotropy.");
            DRAW_WRAPPED(note);
        }

        snprintf(note, sizeof(note),
                 "Ref. slopes: purple k^(-1) Batchelor, "
                 "orange k^(-5/3) Kolmogorov, pink k^(-3) Kraichnan.");
        DRAW_WRAPPED(note);

        #undef DRAW_WRAPPED
    }

    XFreeGC(display, gc2);
    XFlush(display);
}

static void fft2d_canvas_expose_handler(Widget w, XtPointer client_data,
                                         XEvent *event, Boolean *cont) {
    if (event->type != Expose) return;
    draw_fft_spectrum((FFTPopupData *)client_data);
}

void show_2dfft(PlotfileData *pf) {
    FFTPopupData *popup = (FFTPopupData *)calloc(1, sizeof(FFTPopupData));
    popup->pf = pf;
    popup->method = 0;  /* Default to 2D FFT */
    compute_2dfft_spectrum(popup);

    /* More square aspect ratio (spectracam style) */
    Widget sh = XtVaCreatePopupShell("Energy Spectrum",
        transientShellWidgetClass, toplevel,
        XtNwidth, 760, XtNheight, 720, NULL);
    popup->shell = sh;

    Widget frm = XtVaCreateManagedWidget("form", formWidgetClass, sh, NULL);

    /* Method selection buttons */
    Widget method_label = XtVaCreateManagedWidget("Method:",
        labelWidgetClass, frm,
        XtNborderWidth, 0, NULL);

    Widget fft_btn = XtVaCreateManagedWidget("[2D FFT]",
        commandWidgetClass, frm,
        XtNfromHoriz, method_label, NULL);
    XtAddCallback(fft_btn, XtNcallback, fft_method_fft2d_callback, popup);
    popup->method_fft_btn = fft_btn;

    Widget wk_btn = XtVaCreateManagedWidget("Wiener-Khinchin",
        commandWidgetClass, frm,
        XtNfromHoriz, fft_btn, NULL);
    XtAddCallback(wk_btn, XtNcallback, fft_method_wk_callback, popup);
    popup->method_wk_btn = wk_btn;

    Widget cv = XtVaCreateManagedWidget("fftCanvas",
        simpleWidgetClass, frm,
        XtNwidth, 740, XtNheight, 620, XtNborderWidth, 1,
        XtNfromVert, method_label, NULL);
    popup->canvas = cv;
    XtAddEventHandler(cv, ExposureMask, False, fft2d_canvas_expose_handler, popup);

    Widget close_btn = XtVaCreateManagedWidget("Close",
        commandWidgetClass, frm,
        XtNfromVert, cv, NULL);
    XtAddCallback(close_btn, XtNcallback, close_fft2d_popup_callback, popup);

    g_fft_popup = popup;
    XtPopup(sh, XtGrabNone);
}

void fft2d_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data)
        show_2dfft(global_pf);
}

/* Popup data for distribution histogram */
typedef struct {
    Widget shell;
    Widget canvas;
    Widget layer_button;
    Widget domain_button;
    double *bin_counts;
    double *bin_centers;
    int n_bins;
    double count_max;
    double bin_min, bin_max;
    double mean, std, skewness;
    char title[256];
    char xlabel[128];
    int mode;  /* 0=Layer, 1=Domain */
    PlotfileData *pf;  /* Reference to plotfile data */
} DistributionPopupData;

/* Global pointer to current distribution popup */
static DistributionPopupData *g_distrib_popup = NULL;

/* Callback to destroy distribution popup and free data */
void close_distribution_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    DistributionPopupData *popup_data = (DistributionPopupData *)client_data;

    if (popup_data) {
        if (popup_data->bin_counts) free(popup_data->bin_counts);
        if (popup_data->bin_centers) free(popup_data->bin_centers);
        XtDestroyWidget(popup_data->shell);
        free(popup_data);
        if (g_distrib_popup == popup_data) {
            g_distrib_popup = NULL;
        }
    }
}

/* Compute distribution histogram for Layer or Domain mode */
void compute_distribution_data(DistributionPopupData *popup, int mode) {
    PlotfileData *pf = popup->pf;
    if (!pf || !pf->data) return;

    const char *axis_names[] = {"X", "Y", "Z"};
    int axis = pf->slice_axis;
    int slice_idx = pf->slice_idx;

    /* Free old data if exists */
    if (popup->bin_counts) { free(popup->bin_counts); popup->bin_counts = NULL; }
    if (popup->bin_centers) { free(popup->bin_centers); popup->bin_centers = NULL; }

    int data_size;
    double *data_array;

    if (mode == 0) {
        /* Layer mode: extract slice */
        int slice_dim1, slice_dim2;
        if (axis == 2) {
            slice_dim1 = pf->grid_dims[0];
            slice_dim2 = pf->grid_dims[1];
        } else if (axis == 1) {
            slice_dim1 = pf->grid_dims[0];
            slice_dim2 = pf->grid_dims[2];
        } else {
            slice_dim1 = pf->grid_dims[1];
            slice_dim2 = pf->grid_dims[2];
        }
        data_size = slice_dim1 * slice_dim2;
        data_array = (double *)malloc(data_size * sizeof(double));

        int k = 0;
        for (int j = 0; j < slice_dim2; j++) {
            for (int i = 0; i < slice_dim1; i++) {
                int idx;
                if (axis == 2) {
                    idx = slice_idx * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + i;
                } else if (axis == 1) {
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + slice_idx * pf->grid_dims[0] + i;
                } else {
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + i * pf->grid_dims[0] + slice_idx;
                }
                data_array[k++] = pf->data[idx];
            }
        }

        snprintf(popup->title, sizeof(popup->title), "%s Distribution - %s Layer %d (Level %d)",
                 pf->variables[pf->current_var], axis_names[axis], slice_idx + 1, pf->current_level);
    } else {
        /* Domain mode: use entire domain */
        data_size = pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2];
        data_array = (double *)malloc(data_size * sizeof(double));
        for (int i = 0; i < data_size; i++) {
            data_array[i] = pf->data[i];
        }

        snprintf(popup->title, sizeof(popup->title), "%s Distribution - Entire Domain (Level %d)",
                 pf->variables[pf->current_var], pf->current_level);
    }

    /* Calculate statistics */
    double sum = 0.0, sum_sq = 0.0;
    double data_min = 1e30, data_max = -1e30;

    for (int i = 0; i < data_size; i++) {
        double val = data_array[i];
        sum += val;
        sum_sq += val * val;
        if (val < data_min) data_min = val;
        if (val > data_max) data_max = val;
    }

    double mean = sum / data_size;
    double variance = (sum_sq / data_size) - (mean * mean);
    double std = (variance > 0) ? sqrt(variance) : 0.0;

    /* Calculate skewness */
    double sum_third = 0.0;
    for (int i = 0; i < data_size; i++) {
        double diff = data_array[i] - mean;
        sum_third += diff * diff * diff;
    }
    double skewness = 0.0;
    if (std > 0) {
        double std3 = std * std * std;
        skewness = (sum_third / data_size) / std3;
    }

    /* Determine number of bins using Sturges' rule */
    int n_bins = (int)(1 + 3.322 * log10((double)data_size));
    if (n_bins < 10) n_bins = 10;
    if (n_bins > 100) n_bins = 100;

    /* Create histogram */
    double *bin_counts = (double *)calloc(n_bins, sizeof(double));
    double *bin_centers = (double *)malloc(n_bins * sizeof(double));
    double bin_width = (data_max - data_min) / n_bins;
    if (bin_width == 0) bin_width = 1.0;

    for (int i = 0; i < n_bins; i++) {
        bin_centers[i] = data_min + (i + 0.5) * bin_width;
    }

    /* Count values in each bin */
    for (int i = 0; i < data_size; i++) {
        int bin = (int)((data_array[i] - data_min) / bin_width);
        if (bin < 0) bin = 0;
        if (bin >= n_bins) bin = n_bins - 1;
        bin_counts[bin]++;
    }

    /* Find max count for scaling */
    double count_max = 0;
    for (int i = 0; i < n_bins; i++) {
        if (bin_counts[i] > count_max) count_max = bin_counts[i];
    }
    if (count_max == 0) count_max = 1;

    free(data_array);

    /* Store results in popup data */
    popup->bin_counts = bin_counts;
    popup->bin_centers = bin_centers;
    popup->n_bins = n_bins;
    popup->count_max = count_max;
    popup->bin_min = data_min;
    popup->bin_max = data_max;
    popup->mean = mean;
    popup->std = std;
    popup->skewness = skewness;
    popup->mode = mode;
    snprintf(popup->xlabel, sizeof(popup->xlabel), "%s", pf->variables[pf->current_var]);
}

/* Redraw distribution histogram */
void redraw_distribution_histogram(DistributionPopupData *popup) {
    if (!popup || !popup->canvas) return;
    if (!XtIsRealized(popup->canvas)) return;

    Window win = XtWindow(popup->canvas);
    if (!win) return;

    Dimension width, height;
    XtVaGetValues(popup->canvas, XtNwidth, &width, XtNheight, &height, NULL);

    GC plot_gc = XCreateGC(display, win, 0, NULL);
    if (font) XSetFont(display, plot_gc, font->fid);

    draw_histogram(display, win, plot_gc, popup->bin_counts, popup->bin_centers,
                   popup->n_bins, width, height, popup->count_max,
                   popup->bin_min, popup->bin_max, popup->title, popup->xlabel,
                   popup->mean, popup->std, popup->skewness);

    XFreeGC(display, plot_gc);
    XFlush(display);  /* Force immediate display update */
}

/* Distribution mode button callback */
void distrib_mode_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int mode = (int)(long)client_data;
    if (!g_distrib_popup) return;

    compute_distribution_data(g_distrib_popup, mode);
    redraw_distribution_histogram(g_distrib_popup);
}

/* Update distribution popup if it's open (called when layer changes) */
void update_distribution_histogram(int mode) {
    if (!g_distrib_popup) return;

    /* Use current mode if mode is -1, otherwise use specified mode */
    int update_mode = (mode < 0) ? g_distrib_popup->mode : mode;
    compute_distribution_data(g_distrib_popup, update_mode);
    redraw_distribution_histogram(g_distrib_popup);
}

/* Distribution canvas expose handler */
void distrib_canvas_expose_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != Expose) return;
    DistributionPopupData *popup = (DistributionPopupData *)client_data;
    if (popup) {
        redraw_distribution_histogram(popup);
    }
}

/* Draw histogram on a window */
void draw_histogram(Display *dpy, Window win, GC plot_gc, double *bin_counts, double *bin_centers,
                    int n_bins, int width, int height, double count_max,
                    double bin_min, double bin_max, const char *title, const char *xlabel,
                    double mean, double std, double skewness) {
    /* Clear background */
    XSetForeground(dpy, plot_gc, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, plot_gc, 0, 0, width, height);

    /* Draw border */
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, plot_gc, 0, 0, width - 1, height - 1);

    /* Draw title */
    if (font) {
        XSetFont(dpy, plot_gc, font->fid);
        XDrawString(dpy, win, plot_gc, 10, 20, title, strlen(title));
    }

    /* Plot area */
    int plot_left = 70;
    int plot_right = width - 20;
    int plot_top = 40;
    int plot_bottom = height - 80;
    int plot_width = plot_right - plot_left;
    int plot_height = plot_bottom - plot_top;

    if (plot_width <= 0 || plot_height <= 0 || n_bins < 1) return;

    /* Draw axes */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_bottom, plot_right, plot_bottom);  /* x-axis */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_top, plot_left, plot_bottom);      /* y-axis */

    /* Draw y-axis (count) ticks and labels */
    char label[64];
    int num_y_ticks = 4;
    for (int i = 0; i <= num_y_ticks; i++) {
        double y_val = count_max * i / num_y_ticks;
        int y_pos = plot_bottom - (int)(plot_height * i / num_y_ticks);

        XDrawLine(dpy, win, plot_gc, plot_left - 3, y_pos, plot_left, y_pos);
        snprintf(label, sizeof(label), "%.0f", y_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(dpy, win, plot_gc, plot_left - label_width - 5, y_pos + 4, label, strlen(label));
    }

    /* Draw x-axis (value) ticks and labels */
    int num_x_ticks = 5;
    for (int i = 0; i <= num_x_ticks; i++) {
        double x_val = bin_min + (bin_max - bin_min) * i / num_x_ticks;
        int x_pos = plot_left + (int)(plot_width * i / num_x_ticks);

        XDrawLine(dpy, win, plot_gc, x_pos, plot_bottom, x_pos, plot_bottom + 3);
        snprintf(label, sizeof(label), "%.2e", x_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(dpy, win, plot_gc, x_pos - label_width / 2, plot_bottom + 14, label, strlen(label));
    }

    /* Draw x-axis label */
    if (xlabel && xlabel[0]) {
        int xlabel_width = XTextWidth(font, xlabel, strlen(xlabel));
        XDrawString(dpy, win, plot_gc, plot_left + (plot_width - xlabel_width) / 2,
                    plot_bottom + 30, xlabel, strlen(xlabel));
    }

    /* Draw histogram bars */
    XSetForeground(dpy, plot_gc, 0x4444FF);  /* Blue */
    double bin_width = (bin_max - bin_min) / n_bins;
    int bar_width = plot_width / n_bins;
    if (bar_width < 1) bar_width = 1;

    for (int i = 0; i < n_bins; i++) {
        int x = plot_left + (int)((bin_centers[i] - bin_min - bin_width/2) / (bin_max - bin_min) * plot_width);
        int bar_height = (int)(bin_counts[i] / count_max * plot_height);
        if (bar_height < 0) bar_height = 0;
        int y = plot_bottom - bar_height;

        XFillRectangle(dpy, win, plot_gc, x, y, bar_width - 1, bar_height);
    }

    /* Draw statistics text */
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    char stats[256];
    snprintf(stats, sizeof(stats), "Mean: %.4e   Std: %.4e   Skewness: %.4f", mean, std, skewness);
    XDrawString(dpy, win, plot_gc, plot_left, plot_bottom + 55, stats, strlen(stats));

    XFlush(dpy);
}

/* Draw SDM histogram with um units, log scale support, kurtosis, cutoff info */
void draw_sdm_histogram(Display *dpy, Window win, GC plot_gc, HistogramData *hist,
                         int width, int height, int log_x, int log_y, const char *ylabel) {
    /* Clear background */
    XSetForeground(dpy, plot_gc, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, plot_gc, 0, 0, width, height);

    /* Draw border */
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, plot_gc, 0, 0, width - 1, height - 1);

    if (!hist || hist->n_bins < 1) {
        /* No data — show message */
        if (font) XSetFont(dpy, plot_gc, font->fid);
        const char *msg = "No particles in this timestep";
        int mw = font ? XTextWidth(font, msg, strlen(msg)) : 0;
        XDrawString(dpy, win, plot_gc, (width - mw) / 2, height / 2, msg, strlen(msg));
        XFlush(dpy);
        return;
    }

    /* Draw title */
    if (font) {
        XSetFont(dpy, plot_gc, font->fid);
        XDrawString(dpy, win, plot_gc, 10, 20, hist->title, strlen(hist->title));
    }

    /* Plot area with wider left margin for y-axis label */
    int plot_left = 100;
    int plot_right = width - 20;
    int plot_top = 40;
    int plot_bottom = height - 100;
    int plot_width = plot_right - plot_left;
    int plot_height = plot_bottom - plot_top;

    if (plot_width <= 0 || plot_height <= 0) return;

    /* Draw axes */
    XDrawLine(dpy, win, plot_gc, plot_left, plot_bottom, plot_right, plot_bottom);
    XDrawLine(dpy, win, plot_gc, plot_left, plot_top, plot_left, plot_bottom);

    /* Y-axis label (drawn horizontally above y-axis) */
    if (ylabel && ylabel[0]) {
        XDrawString(dpy, win, plot_gc, plot_left, plot_top - 8, ylabel, strlen(ylabel));
    }

    /* Determine Y range */
    double y_max = hist->count_max;
    double y_min_display = 0;
    if (log_y) {
        /* For log Y, find min positive value */
        double min_pos = y_max;
        for (int i = 0; i < hist->n_bins; i++) {
            if (hist->bin_counts[i] > 0 && hist->bin_counts[i] < min_pos)
                min_pos = hist->bin_counts[i];
        }
        y_min_display = pow(10.0, floor(log10(min_pos > 0 ? min_pos : 1)));
        y_max = pow(10.0, ceil(log10(y_max > 0 ? y_max : 1)));
        if (y_min_display >= y_max) y_min_display = y_max / 10.0;
    }

    /* Draw Y-axis ticks */
    char label[64];
    if (log_y) {
        double log_ymin = log10(y_min_display);
        double log_ymax = log10(y_max);
        int imin = (int)floor(log_ymin);
        int imax = (int)ceil(log_ymax);
        for (int i = imin; i <= imax; i++) {
            double y_val = pow(10.0, i);
            if (y_val < y_min_display || y_val > y_max) continue;
            double frac = (log10(y_val) - log_ymin) / (log_ymax - log_ymin);
            int y_pos = plot_bottom - (int)(plot_height * frac);
            XDrawLine(dpy, win, plot_gc, plot_left - 3, y_pos, plot_left, y_pos);
            snprintf(label, sizeof(label), "1e%d", i);
            int lw = XTextWidth(font, label, strlen(label));
            XDrawString(dpy, win, plot_gc, plot_left - lw - 5, y_pos + 4, label, strlen(label));
        }
    } else {
        int num_y_ticks = 4;
        for (int i = 0; i <= num_y_ticks; i++) {
            double y_val = y_max * i / num_y_ticks;
            int y_pos = plot_bottom - (int)(plot_height * i / num_y_ticks);
            XDrawLine(dpy, win, plot_gc, plot_left - 3, y_pos, plot_left, y_pos);
            /* Smart formatting based on value magnitude */
            if (y_val == 0) {
                snprintf(label, sizeof(label), "0");
            } else if (y_val >= 1e6 || y_val < 0.001) {
                snprintf(label, sizeof(label), "%.1e", y_val);
            } else if (y_val < 0.1) {
                snprintf(label, sizeof(label), "%.3f", y_val);
            } else if (y_val < 1) {
                snprintf(label, sizeof(label), "%.2f", y_val);
            } else if (y_val < 100) {
                snprintf(label, sizeof(label), "%.1f", y_val);
            } else {
                snprintf(label, sizeof(label), "%.0f", y_val);
            }
            int lw = XTextWidth(font, label, strlen(label));
            XDrawString(dpy, win, plot_gc, plot_left - lw - 5, y_pos + 4, label, strlen(label));
        }
    }

    /* Determine X range */
    double x_min = hist->bin_min;
    double x_max = hist->bin_max;
    double log_xmin = 0, log_xmax = 1;
    if (log_x) {
        log_xmin = (x_min > 0) ? log10(x_min) : log10(x_max) - 3;
        log_xmax = (x_max > 0) ? log10(x_max) : 0;
        if (log_xmin >= log_xmax) log_xmin = log_xmax - 1;
    } else {
        /* In linear mode, start x-axis from 0 unless user specified a min */
        /* hist->bin_min will be non-zero if user set xlim_min */
        if (x_min <= 0) x_min = 0;
    }

    /* Draw X-axis ticks */
    if (log_x) {
        int imin = (int)floor(log_xmin);
        int imax = (int)ceil(log_xmax);
        for (int i = imin; i <= imax; i++) {
            double x_val = pow(10.0, i);
            double frac = (log10(x_val) - log_xmin) / (log_xmax - log_xmin);
            if (frac < 0 || frac > 1) continue;
            int x_pos = plot_left + (int)(plot_width * frac);
            XDrawLine(dpy, win, plot_gc, x_pos, plot_bottom, x_pos, plot_bottom + 3);
            snprintf(label, sizeof(label), "1e%d", i);
            int lw = XTextWidth(font, label, strlen(label));
            XDrawString(dpy, win, plot_gc, x_pos - lw / 2, plot_bottom + 14, label, strlen(label));
        }
    } else {
        int num_x_ticks = 5;
        for (int i = 0; i <= num_x_ticks; i++) {
            double x_val = x_min + (x_max - x_min) * i / num_x_ticks;
            int x_pos = plot_left + (int)(plot_width * i / num_x_ticks);
            XDrawLine(dpy, win, plot_gc, x_pos, plot_bottom, x_pos, plot_bottom + 3);
            snprintf(label, sizeof(label), "%.2f", x_val);
            int lw = XTextWidth(font, label, strlen(label));
            XDrawString(dpy, win, plot_gc, x_pos - lw / 2, plot_bottom + 14, label, strlen(label));
        }
    }

    /* Draw x-axis label */
    {
        const char *xlab = "radius (um)";
        int xlab_w = XTextWidth(font, xlab, strlen(xlab));
        XDrawString(dpy, win, plot_gc, plot_left + (plot_width - xlab_w) / 2,
                    plot_bottom + 30, xlab, strlen(xlab));
    }

    /* Draw histogram bars */
    XSetForeground(dpy, plot_gc, 0x4444FF);
    double bin_width = (hist->n_bins > 1) ? (x_max - x_min) / hist->n_bins : 1.0;

    double log_ymin_d = log_y ? log10(y_min_display) : 0;
    double log_ymax_d = log_y ? log10(y_max) : y_max;

    for (int i = 0; i < hist->n_bins; i++) {
        double bcount = hist->bin_counts[i];
        if (bcount <= 0 && log_y) continue;

        int bar_x, bar_w, bar_h, bar_y;
        double left_edge, right_edge;

        /* Use explicit bin_edges if available, otherwise compute from center and uniform width */
        if (hist->bin_edges) {
            left_edge = hist->bin_edges[i];
            right_edge = hist->bin_edges[i + 1];
        } else {
            double bc = hist->bin_centers[i];
            left_edge = bc - bin_width / 2;
            right_edge = bc + bin_width / 2;
        }

        if (log_x) {
            if (left_edge <= 0) left_edge = x_min > 0 ? x_min : right_edge / 10;
            if (right_edge <= 0) continue;
            double frac_l = (log10(left_edge) - log_xmin) / (log_xmax - log_xmin);
            double frac_r = (log10(right_edge) - log_xmin) / (log_xmax - log_xmin);
            if (frac_l < 0) frac_l = 0;
            if (frac_r > 1) frac_r = 1;
            bar_x = plot_left + (int)(plot_width * frac_l);
            bar_w = (int)(plot_width * (frac_r - frac_l));
            if (bar_w < 1) bar_w = 1;
        } else {
            double frac_l = (left_edge - x_min) / (x_max - x_min);
            double frac_r = (right_edge - x_min) / (x_max - x_min);
            if (frac_l < 0) frac_l = 0;
            if (frac_r > 1) frac_r = 1;
            bar_x = plot_left + (int)(plot_width * frac_l);
            bar_w = (int)(plot_width * (frac_r - frac_l));
            if (bar_w < 1) bar_w = 1;
        }

        if (log_y) {
            double log_val = log10(bcount);
            double frac = (log_val - log_ymin_d) / (log_ymax_d - log_ymin_d);
            if (frac < 0) frac = 0;
            bar_h = (int)(frac * plot_height);
        } else {
            bar_h = (int)(bcount / y_max * plot_height);
        }
        if (bar_h < 0) bar_h = 0;
        bar_y = plot_bottom - bar_h;

        XFillRectangle(dpy, win, plot_gc, bar_x, bar_y, bar_w, bar_h);
    }

    /* Draw statistics text (two lines) */
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    char stats[512];
    snprintf(stats, sizeof(stats), "Mean: %.4f um   Std: %.4f um   Skew: %.4f   Kurt: %.4f",
             hist->mean, hist->std, hist->skewness, hist->kurtosis);
    XDrawString(dpy, win, plot_gc, plot_left, plot_bottom + 50, stats, strlen(stats));

    /* Second line: cutoff info if applicable */
    if (hist->xlabel[0] && strcmp(hist->xlabel, "radius (um)") != 0) {
        /* xlabel used to carry cutoff info */
        XDrawString(dpy, win, plot_gc, plot_left, plot_bottom + 68, hist->xlabel, strlen(hist->xlabel));
    }

    XFlush(dpy);
}

/* Expose event handler for histogram canvas */
void histogram_expose_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != Expose) return;

    HistogramData *hist_data = (HistogramData *)client_data;
    if (!hist_data || !hist_data->bin_counts) return;

    Window win = XtWindow(w);
    if (!win) return;

    Dimension width, height;
    XtVaGetValues(w, XtNwidth, &width, XtNheight, &height, NULL);

    GC plot_gc = XCreateGC(display, win, 0, NULL);
    draw_histogram(display, win, plot_gc, hist_data->bin_counts, hist_data->bin_centers,
                   hist_data->n_bins, width, height, hist_data->count_max,
                   hist_data->bin_min, hist_data->bin_max, hist_data->title, hist_data->xlabel,
                   hist_data->mean, hist_data->std, hist_data->skewness);
    XFreeGC(display, plot_gc);
}

/* Close callback that also frees histogram data */
void close_histogram_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    HistogramData *hist_data = (HistogramData *)client_data;

    if (hist_data) {
        if (hist_data->bin_counts) free(hist_data->bin_counts);
        if (hist_data->bin_centers) free(hist_data->bin_centers);

        /* Get the shell widget and destroy it */
        Widget shell = XtParent(XtParent(w));
        XtDestroyWidget(shell);

        free(hist_data);
    }
}

/* Show distribution histogram for current slice or domain */
void show_distribution(PlotfileData *pf) {
    /* Create popup data structure */
    DistributionPopupData *popup = (DistributionPopupData *)calloc(1, sizeof(DistributionPopupData));
    popup->pf = pf;
    popup->mode = 0;  /* Start in Layer mode */

    /* Compute initial histogram for Layer mode */
    compute_distribution_data(popup, 0);

    /* Create popup shell */
    Widget popup_shell = XtVaCreatePopupShell("Distribution",
        transientShellWidgetClass, toplevel,
        XtNwidth, 620,
        XtNheight, 420,
        NULL);
    popup->shell = popup_shell;

    Widget popup_form = XtVaCreateManagedWidget("form",
        formWidgetClass, popup_shell,
        NULL);

    /* Mode selection box */
    Widget mode_box = XtVaCreateManagedWidget("modeBox",
        boxWidgetClass, popup_form,
        XtNorientation, XtorientHorizontal,
        XtNborderWidth, 0,
        NULL);

    /* Mode label */
    XtVaCreateManagedWidget("Mode:",
        labelWidgetClass, mode_box,
        XtNborderWidth, 0,
        NULL);

    /* Layer button */
    Widget layer_btn = XtVaCreateManagedWidget("Layer",
        commandWidgetClass, mode_box,
        NULL);
    XtAddCallback(layer_btn, XtNcallback, distrib_mode_callback, (XtPointer)0);
    popup->layer_button = layer_btn;

    /* Domain button */
    Widget domain_btn = XtVaCreateManagedWidget("Domain",
        commandWidgetClass, mode_box,
        NULL);
    XtAddCallback(domain_btn, XtNcallback, distrib_mode_callback, (XtPointer)1);
    popup->domain_button = domain_btn;

    /* Histogram canvas */
    Widget hist_canvas = XtVaCreateManagedWidget("histogram",
        simpleWidgetClass, popup_form,
        XtNfromVert, mode_box,
        XtNwidth, 600,
        XtNheight, 320,
        XtNborderWidth, 1,
        NULL);
    popup->canvas = hist_canvas;

    /* Add expose event handler */
    XtAddEventHandler(hist_canvas, ExposureMask, False, distrib_canvas_expose_handler, popup);

    /* Close button */
    Widget close_button = XtVaCreateManagedWidget("Close",
        commandWidgetClass, popup_form,
        XtNfromVert, hist_canvas,
        NULL);

    XtAddCallback(close_button, XtNcallback, close_distribution_popup_callback, popup);

    /* Set global pointer for mode callbacks */
    g_distrib_popup = popup;

    /* Show popup */
    XtPopup(popup_shell, XtGrabNone);
}

/* Distribution button callback */
void distribution_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data) {
        show_distribution(global_pf);
    }
}

/* Helper function to find variable index by name */
int find_variable_index(PlotfileData *pf, const char *name) {
    for (int i = 0; i < pf->n_vars; i++) {
        if (strcmp(pf->variables[i], name) == 0) {
            return i;
        }
    }
    return -1;
}

/* Helper function to find a velocity component with fallback patterns */
int find_velocity_component(PlotfileData *pf, const char *primary, char fallback_char) {
    int idx;
    char patterns[4][64];
    
    /* Try primary name (e.g., "x_velocity") */
    if ((idx = find_variable_index(pf, primary)) >= 0) {
        return idx;
    }
    
    /* Try simple single letter (e.g., "u", "v", "w") */
    snprintf(patterns[0], sizeof(patterns[0]), "%c", fallback_char);
    if ((idx = find_variable_index(pf, patterns[0])) >= 0) {
        return idx;
    }
    
    /* Try patterns with underscore prefix (e.g., "u_gas", "v_gas", "w_gas") */
    for (int i = 0; i < pf->n_vars; i++) {
        if (strlen(pf->variables[i]) >= 2 && 
            pf->variables[i][0] == fallback_char && 
            pf->variables[i][1] == '_') {
            return i;
        }
    }
    
    return -1;
}

/* Get default component names based on current slice axis */
void get_default_quiver_components(PlotfileData *pf, char *x_comp, char *y_comp) {
    int x_idx, y_idx;
    const char *primary_x, *primary_y;
    char fallback_x, fallback_y;
    
    /* Determine primary names and fallback characters based on slice axis */
    switch (pf->slice_axis) {
        case 0:  /* X plane - show Y and Z velocity */
            primary_x = "y_velocity"; fallback_x = 'v';
            primary_y = "z_velocity"; fallback_y = 'w';
            break;
        case 1:  /* Y plane - show X and Z velocity */
            primary_x = "x_velocity"; fallback_x = 'u';
            primary_y = "z_velocity"; fallback_y = 'w';
            break;
        case 2:  /* Z plane - show X and Y velocity */
        default:
            primary_x = "x_velocity"; fallback_x = 'u';
            primary_y = "y_velocity"; fallback_y = 'v';
            break;
    }
    
    /* Find components with fallback logic */
    x_idx = find_velocity_component(pf, primary_x, fallback_x);
    y_idx = find_velocity_component(pf, primary_y, fallback_y);
    
    /* Set component names based on what was found */
    if (x_idx >= 0) {
        strcpy(x_comp, pf->variables[x_idx]);
    } else {
        strcpy(x_comp, primary_x);  /* Fall back to primary name */
    }
    
    if (y_idx >= 0) {
        strcpy(y_comp, pf->variables[y_idx]);
    } else {
        strcpy(y_comp, primary_y);  /* Fall back to primary name */
    }
}

/* Quiver button callback */
void quiver_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data) {
        /* Get default component names and enable quiver immediately */
        char default_x[64], default_y[64];
        get_default_quiver_components(global_pf, default_x, default_y);
        
        /* Find default indices */
        quiver_data.x_comp_index = find_variable_index(global_pf, default_x);
        quiver_data.y_comp_index = find_variable_index(global_pf, default_y);
        
        if (quiver_data.x_comp_index >= 0 && quiver_data.y_comp_index >= 0) {
            quiver_data.enabled = 1;
            /* Trigger immediate render with default settings */
            render_slice(global_pf);
        } else {
            fprintf(stderr, "Warning: Could not find default velocity components\n");
        }
        
        /* Show dialog for adjustments */
        show_quiver_dialog(global_pf);
    }
}

/* Show quiver component selection dialog */
void show_quiver_dialog(PlotfileData *pf) {
    Arg args[20];
    int n;
    Widget form, label, close_button;
    Widget density_minus, density_plus, scale_minus, scale_plus;
    Widget color_black, color_white, color_red, color_blue;
    char density_text[32], scale_text[32];
    
    /* Don't create multiple dialogs */
    if (quiver_data.shell) {
        XtPopup(quiver_data.shell, XtGrabNone);
        return;
    }
    
    /* Get default component names */
    char default_x[64], default_y[64];
    get_default_quiver_components(pf, default_x, default_y);
    
    /* Find default indices */
    quiver_data.x_comp_index = find_variable_index(pf, default_x);
    quiver_data.y_comp_index = find_variable_index(pf, default_y);
    
    /* Create popup shell */
    quiver_data.shell = XtVaCreatePopupShell("Quiver Options",
        transientShellWidgetClass, toplevel,
        XtNwidth, 400,
        XtNheight, 350,
        NULL);
    
    /* Main form */
    n = 0;
    form = XtCreateManagedWidget("form", formWidgetClass, quiver_data.shell, args, n);
    
    /* Title label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Quiver Options:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    label = XtCreateManagedWidget("titleLabel", labelWidgetClass, form, args, n);
    
    /* X component label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "X Component:"); n++;
    XtSetArg(args[n], XtNfromVert, label); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget x_label = XtCreateManagedWidget("xLabel", labelWidgetClass, form, args, n);
    
    /* X component button */
    n = 0;
    XtSetArg(args[n], XtNlabel, default_x); n++;
    XtSetArg(args[n], XtNfromVert, x_label); n++;
    XtSetArg(args[n], XtNwidth, 150); n++;
    quiver_data.x_comp_text = XtCreateManagedWidget("xCompButton", commandWidgetClass, form, args, n);
    XtAddCallback(quiver_data.x_comp_text, XtNcallback, show_variable_selector, (XtPointer)1);
    
    /* Y component label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Y Component:"); n++;
    XtSetArg(args[n], XtNfromVert, quiver_data.x_comp_text); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget y_label = XtCreateManagedWidget("yLabel", labelWidgetClass, form, args, n);
    
    /* Y component button */
    n = 0;
    XtSetArg(args[n], XtNlabel, default_y); n++;
    XtSetArg(args[n], XtNfromVert, y_label); n++;
    XtSetArg(args[n], XtNwidth, 150); n++;
    quiver_data.y_comp_text = XtCreateManagedWidget("yCompButton", commandWidgetClass, form, args, n);
    XtAddCallback(quiver_data.y_comp_text, XtNcallback, show_variable_selector, (XtPointer)0);
    
    /* Density control */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Density:"); n++;
    XtSetArg(args[n], XtNfromVert, quiver_data.y_comp_text); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget density_title = XtCreateManagedWidget("densityTitle", labelWidgetClass, form, args, n);
    
    /* Density minus button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "-"); n++;
    XtSetArg(args[n], XtNfromVert, density_title); n++;
    XtSetArg(args[n], XtNwidth, 30); n++;
    density_minus = XtCreateManagedWidget("densityMinus", commandWidgetClass, form, args, n);
    XtAddCallback(density_minus, XtNcallback, quiver_density_callback, (XtPointer)-1);
    
    /* Density display */
    snprintf(density_text, sizeof(density_text), "Density: %d", quiver_data.density);
    n = 0;
    XtSetArg(args[n], XtNlabel, density_text); n++;
    XtSetArg(args[n], XtNfromVert, density_title); n++;
    XtSetArg(args[n], XtNfromHoriz, density_minus); n++;
    XtSetArg(args[n], XtNwidth, 100); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    quiver_data.density_label = XtCreateManagedWidget("densityLabel", labelWidgetClass, form, args, n);
    
    /* Density plus button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "+"); n++;
    XtSetArg(args[n], XtNfromVert, density_title); n++;
    XtSetArg(args[n], XtNfromHoriz, quiver_data.density_label); n++;
    XtSetArg(args[n], XtNwidth, 30); n++;
    density_plus = XtCreateManagedWidget("densityPlus", commandWidgetClass, form, args, n);
    XtAddCallback(density_plus, XtNcallback, quiver_density_callback, (XtPointer)1);
    
    /* Scale control */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Scale:"); n++;
    XtSetArg(args[n], XtNfromVert, density_minus); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget scale_title = XtCreateManagedWidget("scaleTitle", labelWidgetClass, form, args, n);
    
    /* Scale minus button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "-"); n++;
    XtSetArg(args[n], XtNfromVert, scale_title); n++;
    XtSetArg(args[n], XtNwidth, 30); n++;
    scale_minus = XtCreateManagedWidget("scaleMinus", commandWidgetClass, form, args, n);
    XtAddCallback(scale_minus, XtNcallback, quiver_scale_callback, (XtPointer)-1);
    
    /* Scale display */
    snprintf(scale_text, sizeof(scale_text), "Scale: %.1f", quiver_data.scale);
    n = 0;
    XtSetArg(args[n], XtNlabel, scale_text); n++;
    XtSetArg(args[n], XtNfromVert, scale_title); n++;
    XtSetArg(args[n], XtNfromHoriz, scale_minus); n++;
    XtSetArg(args[n], XtNwidth, 100); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    quiver_data.scale_label = XtCreateManagedWidget("scaleLabel", labelWidgetClass, form, args, n);
    
    /* Scale plus button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "+"); n++;
    XtSetArg(args[n], XtNfromVert, scale_title); n++;
    XtSetArg(args[n], XtNfromHoriz, quiver_data.scale_label); n++;
    XtSetArg(args[n], XtNwidth, 30); n++;
    scale_plus = XtCreateManagedWidget("scalePlus", commandWidgetClass, form, args, n);
    XtAddCallback(scale_plus, XtNcallback, quiver_scale_callback, (XtPointer)1);
    
    /* Color control */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Color:"); n++;
    XtSetArg(args[n], XtNfromVert, scale_minus); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget color_title = XtCreateManagedWidget("colorTitle", labelWidgetClass, form, args, n);
    
    /* Color buttons */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Black"); n++;
    XtSetArg(args[n], XtNfromVert, color_title); n++;
    XtSetArg(args[n], XtNwidth, 60); n++;
    color_black = XtCreateManagedWidget("colorBlack", commandWidgetClass, form, args, n);
    XtAddCallback(color_black, XtNcallback, quiver_color_callback, (XtPointer)0);
    
    n = 0;
    XtSetArg(args[n], XtNlabel, "White"); n++;
    XtSetArg(args[n], XtNfromVert, color_title); n++;
    XtSetArg(args[n], XtNfromHoriz, color_black); n++;
    XtSetArg(args[n], XtNwidth, 60); n++;
    color_white = XtCreateManagedWidget("colorWhite", commandWidgetClass, form, args, n);
    XtAddCallback(color_white, XtNcallback, quiver_color_callback, (XtPointer)1);
    
    n = 0;
    XtSetArg(args[n], XtNlabel, "Red"); n++;
    XtSetArg(args[n], XtNfromVert, color_title); n++;
    XtSetArg(args[n], XtNfromHoriz, color_white); n++;
    XtSetArg(args[n], XtNwidth, 60); n++;
    color_red = XtCreateManagedWidget("colorRed", commandWidgetClass, form, args, n);
    XtAddCallback(color_red, XtNcallback, quiver_color_callback, (XtPointer)2);
    
    n = 0;
    XtSetArg(args[n], XtNlabel, "Blue"); n++;
    XtSetArg(args[n], XtNfromVert, color_title); n++;
    XtSetArg(args[n], XtNfromHoriz, color_red); n++;
    XtSetArg(args[n], XtNwidth, 60); n++;
    color_blue = XtCreateManagedWidget("colorBlue", commandWidgetClass, form, args, n);
    XtAddCallback(color_blue, XtNcallback, quiver_color_callback, (XtPointer)3);
    
    /* Remove button - removes quiver and closes dialog */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Remove"); n++;
    XtSetArg(args[n], XtNfromVert, color_black); n++;
    close_button = XtCreateManagedWidget("removeButton", commandWidgetClass, form, args, n);
    XtAddCallback(close_button, XtNcallback, quiver_remove_callback, NULL);
    
    /* Show dialog */
    XtPopup(quiver_data.shell, XtGrabNone);
}

/* Apply quiver settings callback */
void quiver_apply_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    /* Component indices are already set by variable selector */
    if (quiver_data.x_comp_index >= 0 && quiver_data.y_comp_index >= 0) {
        quiver_data.enabled = 1;
        /* Trigger redraw to show quiver overlay */
        render_slice(global_pf);
    } else {
        quiver_data.enabled = 0;
        fprintf(stderr, "Warning: Invalid variable selection\n");
    }
}

/* Close quiver dialog callback */
void quiver_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (quiver_data.shell) {
        XtDestroyWidget(quiver_data.shell);
        quiver_data.shell = NULL;
        quiver_data.x_comp_text = NULL;
        quiver_data.y_comp_text = NULL;
        quiver_data.density_label = NULL;
        quiver_data.scale_label = NULL;
    }
}

/* Remove quiver callback - disables quiver and closes dialog */
void quiver_remove_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    /* Disable quiver */
    quiver_data.enabled = 0;
    
    /* Trigger redraw to remove quiver overlay */
    if (global_pf) {
        render_slice(global_pf);
    }
    
    /* Close dialog */
    quiver_close_callback(w, client_data, call_data);
}

/* Density adjustment callback */
void quiver_density_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int direction = (int)(long)client_data;
    char density_text[32];
    
    quiver_data.density += direction;
    if (quiver_data.density < 1) quiver_data.density = 1;
    if (quiver_data.density > 5) quiver_data.density = 5;
    
    snprintf(density_text, sizeof(density_text), "Density: %d", quiver_data.density);
    XtVaSetValues(quiver_data.density_label, XtNlabel, density_text, NULL);
    
    /* Trigger redraw if quiver is enabled */
    if (quiver_data.enabled && global_pf) {
        render_slice(global_pf);
    }
}

/* Scale adjustment callback */
void quiver_scale_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int direction = (int)(long)client_data;
    char scale_text[32];
    
    quiver_data.scale += direction * 0.2;
    if (quiver_data.scale < 0.2) quiver_data.scale = 0.2;
    if (quiver_data.scale > 3.0) quiver_data.scale = 3.0;
    
    snprintf(scale_text, sizeof(scale_text), "Scale: %.1f", quiver_data.scale);
    XtVaSetValues(quiver_data.scale_label, XtNlabel, scale_text, NULL);
    
    /* Trigger redraw if quiver is enabled */
    if (quiver_data.enabled && global_pf) {
        render_slice(global_pf);
    }
}

/* Color selection callback */
void quiver_color_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    quiver_data.color = (int)(long)client_data;
    /* Trigger redraw if quiver is enabled */
    if (quiver_data.enabled && global_pf) {
        render_slice(global_pf);
    }
}

/* Show variable selection popup */
void show_variable_selector(Widget w, XtPointer client_data, XtPointer call_data) {
    int for_x_component = (int)(long)client_data;
    Arg args[20];
    int n;
    Widget form, label, close_button;
    char title[64];

    /* Don't create multiple selectors */
    if (var_select_data.shell) {
        XtDestroyWidget(var_select_data.shell);
        if (var_select_data.var_buttons) {
            free(var_select_data.var_buttons);
            var_select_data.var_buttons = NULL;
        }
    }

    var_select_data.selecting_for_x = for_x_component;
    var_select_data.n_vars = global_pf->n_vars;

    /* Create popup shell */
    snprintf(title, sizeof(title), "Select %s Component", for_x_component ? "X" : "Y");
    var_select_data.shell = XtVaCreatePopupShell(title,
        transientShellWidgetClass, toplevel,
        XtNwidth, 250,
        XtNheight, 300,
        NULL);
    
    /* Main form */
    n = 0;
    form = XtCreateManagedWidget("form", formWidgetClass, var_select_data.shell, args, n);
    
    /* Title label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Available Variables:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    label = XtCreateManagedWidget("titleLabel", labelWidgetClass, form, args, n);
    
    /* Create box widget for variable buttons (simple scrollable list) */
    n = 0;
    XtSetArg(args[n], XtNfromVert, label); n++;
    XtSetArg(args[n], XtNorientation, XtorientVertical); n++;
    XtSetArg(args[n], XtNwidth, 230); n++;
    XtSetArg(args[n], XtNheight, 200); n++;
    Widget var_box = XtCreateManagedWidget("varBox", boxWidgetClass, form, args, n);
    
    /* Create buttons for each variable */
    var_select_data.var_buttons = (Widget *)malloc(global_pf->n_vars * sizeof(Widget));
    
    for (int i = 0; i < global_pf->n_vars; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, global_pf->variables[i]); n++;
        XtSetArg(args[n], XtNwidth, 200); n++;
        var_select_data.var_buttons[i] = XtCreateManagedWidget(global_pf->variables[i], 
                                                              commandWidgetClass, var_box, args, n);
        XtAddCallback(var_select_data.var_buttons[i], XtNcallback, 
                     variable_select_callback, (XtPointer)(long)i);
    }
    
    /* Close button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Close"); n++;
    XtSetArg(args[n], XtNfromVert, var_box); n++;
    close_button = XtCreateManagedWidget("closeButton", commandWidgetClass, form, args, n);
    XtAddCallback(close_button, XtNcallback, variable_selector_close_callback, NULL);
    
    /* Show popup */
    XtPopup(var_select_data.shell, XtGrabNone);
}

/* Variable selection callback */
void variable_select_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int var_index = (int)(long)client_data;
    
    /* Update the appropriate component */
    if (var_select_data.selecting_for_x) {
        quiver_data.x_comp_index = var_index;
        /* Update button label */
        XtVaSetValues(quiver_data.x_comp_text, XtNlabel, global_pf->variables[var_index], NULL);
    } else {
        quiver_data.y_comp_index = var_index;
        /* Update button label */
        XtVaSetValues(quiver_data.y_comp_text, XtNlabel, global_pf->variables[var_index], NULL);
    }
    
    /* Trigger immediate redraw if quiver is enabled */
    if (quiver_data.enabled && global_pf) {
        render_slice(global_pf);
    }
    
    /* Close the selector */
    variable_selector_close_callback(w, NULL, NULL);
}

/* Close variable selector callback */
void variable_selector_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (var_select_data.shell) {
        XtDestroyWidget(var_select_data.shell);
        var_select_data.shell = NULL;
        if (var_select_data.var_buttons) {
            free(var_select_data.var_buttons);
            var_select_data.var_buttons = NULL;
        }
    }
}

/* Draw an arrow from (x1,y1) to (x2,y2) */
void draw_arrow(Display *dpy, Drawable win, GC graphics_gc, int x1, int y1, int x2, int y2) {
    /* Draw main line */
    XDrawLine(dpy, win, graphics_gc, x1, y1, x2, y2);
    
    /* Calculate arrow head */
    double angle = atan2(y2 - y1, x2 - x1);
    double head_len = 4.0;  /* Arrow head length */
    double head_angle = 0.5;  /* Arrow head angle */
    
    int head_x1 = x2 - (int)(head_len * cos(angle - head_angle));
    int head_y1 = y2 - (int)(head_len * sin(angle - head_angle));
    int head_x2 = x2 - (int)(head_len * cos(angle + head_angle));
    int head_y2 = y2 - (int)(head_len * sin(angle + head_angle));
    
    /* Draw arrow head */
    XDrawLine(dpy, win, graphics_gc, x2, y2, head_x1, head_y1);
    XDrawLine(dpy, win, graphics_gc, x2, y2, head_x2, head_y2);
}

/* Render quiver overlay */
void render_quiver_overlay(PlotfileData *pf) {
    if (!quiver_data.enabled || quiver_data.x_comp_index < 0 || quiver_data.y_comp_index < 0) {
        return;
    }
    
    /* Get current slice dimensions */
    int width, height;
    if (pf->slice_axis == 2) {
        width = pf->grid_dims[0];
        height = pf->grid_dims[1];
    } else if (pf->slice_axis == 1) {
        width = pf->grid_dims[0];
        height = pf->grid_dims[2];
    } else {
        width = pf->grid_dims[1];
        height = pf->grid_dims[2];
    }
    
    /* Read component data */
    double *x_comp_data = (double *)malloc(pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2] * sizeof(double));
    double *y_comp_data = (double *)malloc(pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2] * sizeof(double));
    
    /* Save current variable and read component data */
    int saved_var = pf->current_var;
    
    pf->current_var = quiver_data.x_comp_index;
    read_variable_data(pf, quiver_data.x_comp_index);
    memcpy(x_comp_data, pf->data, pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2] * sizeof(double));
    
    pf->current_var = quiver_data.y_comp_index;
    read_variable_data(pf, quiver_data.y_comp_index);
    memcpy(y_comp_data, pf->data, pf->grid_dims[0] * pf->grid_dims[1] * pf->grid_dims[2] * sizeof(double));
    
    /* Restore original variable */
    pf->current_var = saved_var;
    read_variable_data(pf, saved_var);
    
    /* Extract slices for both components */
    double *x_slice = (double *)malloc(width * height * sizeof(double));
    double *y_slice = (double *)malloc(width * height * sizeof(double));
    
    extract_slice_from_data(x_comp_data, pf, x_slice, pf->slice_axis, pf->slice_idx);
    extract_slice_from_data(y_comp_data, pf, y_slice, pf->slice_axis, pf->slice_idx);

    /* Map coordinates when map mode is enabled */
    int use_map_coords = 0;
    double *x_coord_slice = NULL;
    double *y_coord_slice = NULL;
    if (pf->map_mode && map_has_bounds) {
        int lon_idx = find_variable_index(pf, "lon_m");
        int lat_idx = find_variable_index(pf, "lat_m");
        if (lon_idx >= 0 && lat_idx >= 0) {
            x_coord_slice = (double *)malloc(width * height * sizeof(double));
            y_coord_slice = (double *)malloc(width * height * sizeof(double));

            if (pf->slice_axis == 2) {
                /* Z-slice: lon/lat */
                read_variable_data(pf, lon_idx);
                extract_slice_from_data(pf->data, pf, x_coord_slice, pf->slice_axis, pf->slice_idx);
                read_variable_data(pf, lat_idx);
                extract_slice_from_data(pf->data, pf, y_coord_slice, pf->slice_axis, pf->slice_idx);
            } else if (pf->slice_axis == 1) {
                /* Y-slice: lon vs Z */
                read_variable_data(pf, lon_idx);
                extract_slice_from_data(pf->data, pf, x_coord_slice, pf->slice_axis, pf->slice_idx);
                if (pf->use_z_phys && current_z_phys_slice) {
                    memcpy(y_coord_slice, current_z_phys_slice,
                           (size_t)width * height * sizeof(double));
                } else {
                    for (int jj = 0; jj < height; jj++) {
                        for (int ii = 0; ii < width; ii++) {
                            int idx = jj * width + ii;
                            double z_coord = pf->prob_lo[2] + (jj + 0.5) *
                                (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                            y_coord_slice[idx] = z_coord;
                        }
                    }
                }
            } else {
                /* X-slice: lat vs Z */
                read_variable_data(pf, lat_idx);
                extract_slice_from_data(pf->data, pf, x_coord_slice, pf->slice_axis, pf->slice_idx);
                if (pf->use_z_phys && current_z_phys_slice) {
                    memcpy(y_coord_slice, current_z_phys_slice,
                           (size_t)width * height * sizeof(double));
                } else {
                    for (int jj = 0; jj < height; jj++) {
                        for (int ii = 0; ii < width; ii++) {
                            int idx = jj * width + ii;
                            double z_coord = pf->prob_lo[2] + (jj + 0.5) *
                                (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                            y_coord_slice[idx] = z_coord;
                        }
                    }
                }
            }

            /* Restore current variable */
            read_variable_data(pf, saved_var);
            use_map_coords = 1;
        }
    }
    
    /* Find max magnitude for scaling */
    double max_mag = 0.0;
    for (int i = 0; i < width * height; i++) {
        double mag = sqrt(x_slice[i] * x_slice[i] + y_slice[i] * y_slice[i]);
        if (mag > max_mag) max_mag = mag;
    }
    
    if (max_mag == 0.0) {
        free(x_comp_data);
        free(y_comp_data);
        free(x_slice);
        free(y_slice);
        return;
    }
    
    /* Set up drawing parameters */
    unsigned long arrow_color;
    switch (quiver_data.color) {
        case 1: arrow_color = WhitePixel(display, screen); break;  /* White */
        case 2: arrow_color = 0xFF0000; break;  /* Red */
        case 3: arrow_color = 0x0000FF; break;  /* Blue */
        default: arrow_color = BlackPixel(display, screen); break;  /* Black */
    }
    XSetForeground(display, gc, arrow_color);
    XSetLineAttributes(display, gc, 1, LineSolid, CapRound, JoinRound);
    
    /* Draw arrows with user-controlled density */
    /* Map density 1-5 to skip values with much wider range */
    int skip;
    switch (quiver_data.density) {
        case 1: skip = (width > 100 || height > 100) ? 20 : 16; break;  /* Very sparse */
        case 2: skip = (width > 100 || height > 100) ? 12 : 10; break;  /* Sparse */
        case 3: skip = (width > 100 || height > 100) ? 8 : 6; break;    /* Medium */
        case 4: skip = (width > 100 || height > 100) ? 5 : 4; break;    /* Dense */
        case 5: skip = (width > 100 || height > 100) ? 3 : 2; break;    /* Very dense */
        default: skip = 8; break;
    }
    
    double scale = 15.0 * quiver_data.scale;  /* User-controlled arrow scale */
    
    for (int j = skip/2; j < height; j += skip) {
        for (int i = skip/2; i < width; i += skip) {
            int idx = j * width + i;
            double u = x_slice[idx] / max_mag;
            double v = y_slice[idx] / max_mag;
            
            if (fabs(u) < 1e-10 && fabs(v) < 1e-10) continue;
            
            int screen_x, screen_y;
            int arrow_dx, arrow_dy;
            if (use_map_coords) {
                double x_coord = x_coord_slice[idx];
                double y_coord = y_coord_slice[idx];
                if (x_coord < map_last_lon_min || x_coord > map_last_lon_max ||
                    y_coord < map_last_lat_min || y_coord > map_last_lat_max) {
                    continue;
                }
                screen_x = render_offset_x + (int)((x_coord - map_last_lon_min) /
                                                  (map_last_lon_max - map_last_lon_min) * render_width);
                screen_y = render_offset_y + (int)((map_last_lat_max - y_coord) /
                                                  (map_last_lat_max - map_last_lat_min) * render_height);

                int i_prev = (i > 0) ? i - 1 : i;
                int i_next = (i + 1 < width) ? i + 1 : i;
                int j_prev = (j > 0) ? j - 1 : j;
                int j_next = (j + 1 < height) ? j + 1 : j;
                if (i_prev == i_next || j_prev == j_next) continue;

                int idx_i_prev = j * width + i_prev;
                int idx_i_next = j * width + i_next;
                int idx_j_prev = j_prev * width + i;
                int idx_j_next = j_next * width + i;

                double xi_prev = x_coord_slice[idx_i_prev];
                double yi_prev = y_coord_slice[idx_i_prev];
                double xi_next = x_coord_slice[idx_i_next];
                double yi_next = y_coord_slice[idx_i_next];
                double xj_prev = x_coord_slice[idx_j_prev];
                double yj_prev = y_coord_slice[idx_j_prev];
                double xj_next = x_coord_slice[idx_j_next];
                double yj_next = y_coord_slice[idx_j_next];

                int sx_i_prev = render_offset_x + (int)((xi_prev - map_last_lon_min) /
                                                       (map_last_lon_max - map_last_lon_min) * render_width);
                int sy_i_prev = render_offset_y + (int)((map_last_lat_max - yi_prev) /
                                                       (map_last_lat_max - map_last_lat_min) * render_height);
                int sx_i_next = render_offset_x + (int)((xi_next - map_last_lon_min) /
                                                       (map_last_lon_max - map_last_lon_min) * render_width);
                int sy_i_next = render_offset_y + (int)((map_last_lat_max - yi_next) /
                                                       (map_last_lat_max - map_last_lat_min) * render_height);

                int sx_j_prev = render_offset_x + (int)((xj_prev - map_last_lon_min) /
                                                       (map_last_lon_max - map_last_lon_min) * render_width);
                int sy_j_prev = render_offset_y + (int)((map_last_lat_max - yj_prev) /
                                                       (map_last_lat_max - map_last_lat_min) * render_height);
                int sx_j_next = render_offset_x + (int)((xj_next - map_last_lon_min) /
                                                       (map_last_lon_max - map_last_lon_min) * render_width);
                int sy_j_next = render_offset_y + (int)((map_last_lat_max - yj_next) /
                                                       (map_last_lat_max - map_last_lat_min) * render_height);

                double basis_ix = 0.5 * (sx_i_next - sx_i_prev);
                double basis_iy = 0.5 * (sy_i_next - sy_i_prev);
                double basis_jx = 0.5 * (sx_j_next - sx_j_prev);
                double basis_jy = 0.5 * (sy_j_next - sy_j_prev);

                double mag_i = sqrt(basis_ix * basis_ix + basis_iy * basis_iy);
                double mag_j = sqrt(basis_jx * basis_jx + basis_jy * basis_jy);
                if (mag_i < 1e-6 || mag_j < 1e-6) continue;

                basis_ix /= mag_i;
                basis_iy /= mag_i;
                basis_jx /= mag_j;
                basis_jy /= mag_j;

                arrow_dx = (int)(scale * (u * basis_ix + v * basis_jx));
                arrow_dy = (int)(scale * (u * basis_iy + v * basis_jy));
            } else {
                /* Convert data coordinates to screen coordinates */
                screen_x = render_offset_x + (int)((double)i * render_width / width);
                if (render_uses_z_phys && current_z_phys_slice) {
                    double z = current_z_phys_slice[idx];
                    screen_y = render_offset_y + (int)((render_phys_ymax - z) /
                        (render_phys_ymax - render_phys_ymin) * render_height);
                } else {
                    /* Flip Y to match regular-grid image rendering. */
                    int flipped_j = height - 1 - j;
                    screen_y = render_offset_y + (int)((double)flipped_j * render_height / height);
                }

                arrow_dx = (int)(u * scale);
                arrow_dy = (int)(-v * scale);  /* Flip Y to match screen coordinates */
            }
            
            draw_arrow(display, canvas, gc, screen_x, screen_y, 
                      screen_x + arrow_dx, screen_y + arrow_dy);
        }
    }
    
    /* Cleanup */
    free(x_comp_data);
    free(y_comp_data);
    free(x_slice);
    free(y_slice);
    if (x_coord_slice) free(x_coord_slice);
    if (y_coord_slice) free(y_coord_slice);
}

/* Render map overlay with US coastline */
static int draw_geojson_coastline(const char *path, double lon_min, double lon_max, double lat_min, double lat_max,
                                 int offset_x, int offset_y, int render_w, int render_h, GC coastline_gc) {
    FILE *fp = fopen(path, "r");
    if (!fp) return 0;

    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    if (fsize <= 0) {
        fclose(fp);
        return 0;
    }
    fseek(fp, 0, SEEK_SET);

    char *buf = (char *)malloc((size_t)fsize + 1);
    if (!buf) {
        fclose(fp);
        return 0;
    }
    size_t nread = fread(buf, 1, (size_t)fsize, fp);
    buf[nread] = '\0';
    fclose(fp);

    int use_360 = (lon_min >= 0.0 && lon_max > 180.0);

    int depth = 0;
    int in_coords = 0;
    int coords_pending = 0;
    int coords_depth = -1;
    int line_depth = -1;

    double prev_lon = 0.0, prev_lat = 0.0;
    int have_prev = 0;

    double point_vals[2];
    int nums_in_point = 0;

    for (char *p = buf; *p; p++) {
        if (!in_coords) {
            if (*p == 'c' && strncmp(p, "coordinates", 11) == 0) {
                coords_pending = 1;
                p += 10;
                continue;
            }
        }

        if (*p == '[') {
            depth++;
            if (coords_pending && !in_coords) {
                in_coords = 1;
                coords_pending = 0;
                coords_depth = depth;
                line_depth = -1;
                have_prev = 0;
                nums_in_point = 0;
            }
            continue;
        }

        if (*p == ']') {
            depth--;
            if (in_coords) {
                if (line_depth >= 0 && depth < line_depth) {
                    have_prev = 0;
                    line_depth = -1;
                    nums_in_point = 0;
                }
                if (coords_depth >= 0 && depth < coords_depth) {
                    in_coords = 0;
                    coords_depth = -1;
                    line_depth = -1;
                    have_prev = 0;
                    nums_in_point = 0;
                }
            }
            continue;
        }

        if (in_coords && (*p == '-' || (*p >= '0' && *p <= '9'))) {
            char *endptr = NULL;
            double val = strtod(p, &endptr);
            if (endptr && endptr != p) {
                if (line_depth < 0) line_depth = depth - 1;

                point_vals[nums_in_point++] = val;
                if (nums_in_point == 2) {
                    double lon = point_vals[0];
                    double lat = point_vals[1];

                    if (use_360 && lon < 0.0) lon += 360.0;

                    if (have_prev) {
                        double dlon = fabs(lon - prev_lon);
                        double dlat = fabs(lat - prev_lat);
                        if (dlon > 30.0 || dlat > 30.0) {
                            have_prev = 0;
                        } else {
                            if ((lon >= lon_min && lon <= lon_max && lat >= lat_min && lat <= lat_max) ||
                                (prev_lon >= lon_min && prev_lon <= lon_max && prev_lat >= lat_min && prev_lat <= lat_max)) {
                                int x1 = offset_x + (int)((prev_lon - lon_min) / (lon_max - lon_min) * render_w);
                                int y1 = offset_y + (int)((lat_max - prev_lat) / (lat_max - lat_min) * render_h);
                                int x2 = offset_x + (int)((lon - lon_min) / (lon_max - lon_min) * render_w);
                                int y2 = offset_y + (int)((lat_max - lat) / (lat_max - lat_min) * render_h);
                                XDrawLine(display, canvas, coastline_gc, x1, y1, x2, y2);
                            }
                        }
                    }

                    prev_lon = lon;
                    prev_lat = lat;
                    have_prev = 1;
                    nums_in_point = 0;
                }

                p = endptr - 1;
            }
        }
    }

    free(buf);
    return 1;
}

void render_map_overlay(PlotfileData *pf, double lon_min, double lon_max, double lat_min, double lat_max) {
    /* Use the same rendering area as the data */
    extern int render_offset_x, render_offset_y, render_width, render_height;
    
    printf("Map overlay: bounds [%.2f,%.2f] x [%.2f,%.2f], render area %dx%d at (%d,%d)\n", 
           lon_min, lon_max, lat_min, lat_max, render_width, render_height, render_offset_x, render_offset_y);
    
    /* Create GC for coastline drawing */
    GC coastline_gc = XCreateGC(display, canvas, 0, NULL);
    if (map_color_pixel == 0) update_map_color_pixel();
    XSetForeground(display, coastline_gc, map_color_pixel);
    XSetLineAttributes(display, coastline_gc, 3, LineSolid, CapButt, JoinMiter);  /* Thick line for visibility */

    /* Apply zoom clip to coastline GC */
    if (zoom_level > 1.0) {
        XRectangle clip = {vis_area_x, vis_area_y, vis_area_w, vis_area_h};
        XSetClipRectangles(display, coastline_gc, 0, 0, &clip, 1, Unsorted);
    }
    
    if (!map_coastlines_enabled) {
        XFreeGC(display, coastline_gc);
        return;
    }

    /* Prefer high-resolution GeoJSON coastlines when available */
    scan_coastline_files();
    int drew_any = 0;
    for (int i = 0; i < n_coastlines; i++) {
        CoastlineEntry *ce = &coastlines[i];
        if (!ce->enabled) continue;
        if (draw_geojson_coastline(ce->filename, lon_min, lon_max, lat_min, lat_max,
                                   render_offset_x, render_offset_y, render_width, render_height,
                                   coastline_gc)) {
            drew_any = 1;
        }
    }

    XFreeGC(display, coastline_gc);
}

/* Helper function to extract slice from arbitrary data array */
void extract_slice_from_data(double *data, PlotfileData *pf, double *slice, int axis, int idx) {
    int width, height;
    if (axis == 2) {
        width = pf->grid_dims[0];
        height = pf->grid_dims[1];
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                slice[j * width + i] = data[idx * width * height + j * width + i];
            }
        }
    } else if (axis == 1) {
        width = pf->grid_dims[0];
        height = pf->grid_dims[2];
        for (int k = 0; k < height; k++) {
            for (int i = 0; i < width; i++) {
                slice[k * width + i] = data[k * pf->grid_dims[0] * pf->grid_dims[1] + idx * pf->grid_dims[0] + i];
            }
        }
    } else {
        width = pf->grid_dims[1];
        height = pf->grid_dims[2];
        for (int k = 0; k < height; k++) {
            for (int j = 0; j < width; j++) {
                slice[k * width + j] = data[k * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + idx];
            }
        }
    }
}

/* Popup data for time series (3 plots) */
typedef struct {
    Widget shell;
    PlotData *mean_plot;
    PlotData *std_plot;
    PlotData *skewness_plot;
} TimeSeriesPopupData;

/* Callback to destroy time series popup and free data */
void close_time_series_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    TimeSeriesPopupData *popup_data = (TimeSeriesPopupData *)client_data;

    if (popup_data) {
        if (popup_data->mean_plot) {
            if (popup_data->mean_plot->data) free(popup_data->mean_plot->data);
            if (popup_data->mean_plot->x_values) free(popup_data->mean_plot->x_values);
            free(popup_data->mean_plot);
        }
        if (popup_data->std_plot) {
            if (popup_data->std_plot->data) free(popup_data->std_plot->data);
            /* Note: std_plot->x_values is shared with mean_plot, already freed */
            free(popup_data->std_plot);
        }
        if (popup_data->skewness_plot) {
            if (popup_data->skewness_plot->data) free(popup_data->skewness_plot->data);
            /* Note: skewness_plot->x_values is shared with mean_plot, already freed */
            free(popup_data->skewness_plot);
        }
        XtDestroyWidget(popup_data->shell);
        free(popup_data);
    }
}

/* Show time series statistics (mean, std, skewness) for the current slice across all timesteps */
void show_time_series(PlotfileData *pf) {
    if (n_timesteps <= 1) return;

    const char *axis_names[] = {"X", "Y", "Z"};
    int axis = pf->slice_axis;
    int slice_idx = pf->slice_idx;
    int current_var = pf->current_var;

    /* Save current state */
    int original_timestep = current_timestep;
    char original_dir[MAX_PATH];
    strncpy(original_dir, pf->plotfile_dir, MAX_PATH - 1);

    /* Determine slice dimensions */
    int slice_dim1, slice_dim2;
    if (axis == 2) {
        slice_dim1 = pf->grid_dims[0];
        slice_dim2 = pf->grid_dims[1];
    } else if (axis == 1) {
        slice_dim1 = pf->grid_dims[0];
        slice_dim2 = pf->grid_dims[2];
    } else {
        slice_dim1 = pf->grid_dims[1];
        slice_dim2 = pf->grid_dims[2];
    }
    int slice_size = slice_dim1 * slice_dim2;

    /* Allocate arrays for time series statistics */
    double *means = (double *)malloc(n_timesteps * sizeof(double));
    double *stds = (double *)malloc(n_timesteps * sizeof(double));
    double *skewness = (double *)malloc(n_timesteps * sizeof(double));
    double *time_indices = (double *)malloc(n_timesteps * sizeof(double));

    printf("Computing time series statistics for %d timesteps...\n", n_timesteps);

    /* Loop through all timesteps */
    for (int t = 0; t < n_timesteps; t++) {
        time_indices[t] = t + 1;  /* 1-indexed for display */

        /* Load this timestep's data */
        strncpy(pf->plotfile_dir, timestep_paths[t], MAX_PATH - 1);
        read_header(pf);
        pf->n_boxes = 0;
        read_cell_h(pf);
        read_variable_data(pf, current_var);

        /* Calculate statistics for the slice */
        double sum = 0.0, sum_sq = 0.0;

        for (int j = 0; j < slice_dim2; j++) {
            for (int i = 0; i < slice_dim1; i++) {
                int idx;
                if (axis == 2) {
                    idx = slice_idx * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + i;
                } else if (axis == 1) {
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + slice_idx * pf->grid_dims[0] + i;
                } else {
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + i * pf->grid_dims[0] + slice_idx;
                }
                double val = pf->data[idx];
                sum += val;
                sum_sq += val * val;
            }
        }

        means[t] = sum / slice_size;
        double variance = (sum_sq / slice_size) - (means[t] * means[t]);
        stds[t] = (variance > 0) ? sqrt(variance) : 0.0;

        /* Second pass: calculate skewness (third moment) */
        double sum_third = 0.0;
        for (int j = 0; j < slice_dim2; j++) {
            for (int i = 0; i < slice_dim1; i++) {
                int idx;
                if (axis == 2) {
                    idx = slice_idx * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + i;
                } else if (axis == 1) {
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + slice_idx * pf->grid_dims[0] + i;
                } else {
                    idx = j * pf->grid_dims[0] * pf->grid_dims[1] + i * pf->grid_dims[0] + slice_idx;
                }
                double val = pf->data[idx];
                double diff = val - means[t];
                sum_third += diff * diff * diff;
            }
        }

        /* Skewness = E[(X - mu)^3] / sigma^3 */
        if (stds[t] > 0) {
            double std3 = stds[t] * stds[t] * stds[t];
            skewness[t] = (sum_third / slice_size) / std3;
        } else {
            skewness[t] = 0.0;
        }

        if ((t + 1) % 10 == 0 || t == n_timesteps - 1) {
            printf("  Processed %d/%d timesteps\n", t + 1, n_timesteps);
        }
    }

    /* Restore original state */
    strncpy(pf->plotfile_dir, original_dir, MAX_PATH - 1);
    current_timestep = original_timestep;
    read_header(pf);
    pf->n_boxes = 0;
    read_cell_h(pf);
    read_variable_data(pf, current_var);

    /* Create plot data for mean */
    PlotData *mean_plot = (PlotData *)malloc(sizeof(PlotData));
    mean_plot->n_points = n_timesteps;
    mean_plot->data = means;
    mean_plot->x_values = (double *)malloc(n_timesteps * sizeof(double));
    memcpy(mean_plot->x_values, time_indices, n_timesteps * sizeof(double));
    mean_plot->vmin = 1e30;
    mean_plot->vmax = -1e30;
    for (int i = 0; i < n_timesteps; i++) {
        if (means[i] < mean_plot->vmin) mean_plot->vmin = means[i];
        if (means[i] > mean_plot->vmax) mean_plot->vmax = means[i];
    }
    mean_plot->xmin = 1;
    mean_plot->xmax = n_timesteps;
    snprintf(mean_plot->title, sizeof(mean_plot->title), "%s Mean (%s Layer %d)",
             pf->variables[pf->current_var], axis_names[axis], slice_idx + 1);
    snprintf(mean_plot->xlabel, sizeof(mean_plot->xlabel), "Timestep");
    snprintf(mean_plot->vlabel, sizeof(mean_plot->vlabel), "Mean");

    /* Create plot data for std */
    PlotData *std_plot = (PlotData *)malloc(sizeof(PlotData));
    std_plot->n_points = n_timesteps;
    std_plot->data = stds;
    std_plot->x_values = time_indices;  /* Share with mean, will be freed once */
    std_plot->vmin = 1e30;
    std_plot->vmax = -1e30;
    for (int i = 0; i < n_timesteps; i++) {
        if (stds[i] < std_plot->vmin) std_plot->vmin = stds[i];
        if (stds[i] > std_plot->vmax) std_plot->vmax = stds[i];
    }
    std_plot->xmin = 1;
    std_plot->xmax = n_timesteps;
    snprintf(std_plot->title, sizeof(std_plot->title), "%s Std Dev (%s Layer %d)",
             pf->variables[pf->current_var], axis_names[axis], slice_idx + 1);
    snprintf(std_plot->xlabel, sizeof(std_plot->xlabel), "Timestep");
    snprintf(std_plot->vlabel, sizeof(std_plot->vlabel), "Std Dev");

    /* Create plot data for skewness */
    PlotData *skewness_plot = (PlotData *)malloc(sizeof(PlotData));
    skewness_plot->n_points = n_timesteps;
    skewness_plot->data = skewness;
    skewness_plot->x_values = time_indices;  /* Share with mean, will be freed once */
    skewness_plot->vmin = 1e30;
    skewness_plot->vmax = -1e30;
    for (int i = 0; i < n_timesteps; i++) {
        if (skewness[i] < skewness_plot->vmin) skewness_plot->vmin = skewness[i];
        if (skewness[i] > skewness_plot->vmax) skewness_plot->vmax = skewness[i];
    }
    skewness_plot->xmin = 1;
    skewness_plot->xmax = n_timesteps;
    snprintf(skewness_plot->title, sizeof(skewness_plot->title), "%s Skewness (%s Layer %d)",
             pf->variables[pf->current_var], axis_names[axis], slice_idx + 1);
    snprintf(skewness_plot->xlabel, sizeof(skewness_plot->xlabel), "Timestep");
    snprintf(skewness_plot->vlabel, sizeof(skewness_plot->vlabel), "Skewness");

    /* Create popup data structure */
    TimeSeriesPopupData *popup_data = (TimeSeriesPopupData *)malloc(sizeof(TimeSeriesPopupData));
    popup_data->mean_plot = mean_plot;
    popup_data->std_plot = std_plot;
    popup_data->skewness_plot = skewness_plot;

    /* Create popup shell - wider for 3 side-by-side plots */
    Widget popup_shell = XtVaCreatePopupShell("Time Series Statistics",
        transientShellWidgetClass, toplevel,
        XtNwidth, 1200,
        XtNheight, 450,
        NULL);

    popup_data->shell = popup_shell;

    Widget popup_form = XtVaCreateManagedWidget("form",
        formWidgetClass, popup_shell,
        NULL);

    /* Mean plot canvas - left */
    Widget mean_canvas = XtVaCreateManagedWidget("mean_plot",
        simpleWidgetClass, popup_form,
        XtNwidth, 380,
        XtNheight, 350,
        XtNborderWidth, 1,
        NULL);

    /* Std plot canvas - middle */
    Widget std_canvas = XtVaCreateManagedWidget("std_plot",
        simpleWidgetClass, popup_form,
        XtNfromHoriz, mean_canvas,
        XtNwidth, 380,
        XtNheight, 350,
        XtNborderWidth, 1,
        NULL);

    /* Skewness plot canvas - right */
    Widget skewness_canvas = XtVaCreateManagedWidget("skewness_plot",
        simpleWidgetClass, popup_form,
        XtNfromHoriz, std_canvas,
        XtNwidth, 380,
        XtNheight, 350,
        XtNborderWidth, 1,
        NULL);

    /* Add expose event handlers - using standard plot (timestep on X, value on Y) */
    XtAddEventHandler(mean_canvas, ExposureMask, False, plot_expose_handler, mean_plot);
    XtAddEventHandler(std_canvas, ExposureMask, False, plot_expose_handler, std_plot);
    XtAddEventHandler(skewness_canvas, ExposureMask, False, plot_expose_handler, skewness_plot);

    /* Close button */
    Widget close_button = XtVaCreateManagedWidget("Close",
        commandWidgetClass, popup_form,
        XtNfromVert, mean_canvas,
        NULL);

    XtAddCallback(close_button, XtNcallback, close_time_series_popup_callback, popup_data);

    /* Show popup */
    XtPopup(popup_shell, XtGrabNone);

    printf("Time series statistics displayed.\n");
}

/* --------------------------- Video playback --------------------------- */
static double video_value_for_scale(double v, int mode) {
    if (mode == 1) return log10(v);
    if (mode == 2) return log10(-v);
    return v;
}

static int video_value_valid(double v, int mode) {
    return isfinite(v) && (mode == 0 || (mode == 1 ? v > 0.0 : v < 0.0));
}

static void video_status(VideoState *vs, const char *text) {
    if (vs && vs->status_label) XtVaSetValues(vs->status_label, XtNlabel, text, NULL);
}

static int video_prepare_plot(VideoState *vs, int ti, PlotfileData *tmp,
                              int *var_idx, int *slice_idx) {
    memset(tmp, 0, sizeof(*tmp));
    strncpy(tmp->plotfile_dir, timestep_paths[ti], MAX_PATH - 1);
    if (read_header(tmp) < 0) return -1;
    *var_idx = find_variable_index(tmp, vs->variable_name);
    if (*var_idx < 0) return -1;
    tmp->current_level = vs->requested_level;
    if (tmp->current_level >= tmp->n_levels) {
        fprintf(stderr, "Video: timestep %d has fewer levels; using level %d\n",
                ti + 1, tmp->n_levels - 1);
        tmp->current_level = tmp->n_levels - 1;
    }
    if (tmp->current_level < 0) tmp->current_level = 0;
    if (read_cell_h(tmp) < 0 || read_variable_data(tmp, *var_idx) < 0) return -1;
    tmp->slice_axis = vs->slice_axis;
    tmp->overlay_mode = vs->overlay_mode;
    tmp->map_mode = vs->map_mode;
    tmp->colormap = vs->colormap;
    double span = tmp->prob_hi[vs->slice_axis] - tmp->prob_lo[vs->slice_axis];
    double cells = (double)(tmp->level_lo[vs->slice_axis] + tmp->grid_dims[vs->slice_axis]);
    double dx = (cells > 0.0) ? span / cells : span;
    *slice_idx = (int)floor((vs->slice_position - tmp->prob_lo[vs->slice_axis]) / dx)
                 - tmp->level_lo[vs->slice_axis];
    if (*slice_idx < 0) *slice_idx = 0;
    if (*slice_idx >= tmp->grid_dims[vs->slice_axis])
        *slice_idx = tmp->grid_dims[vs->slice_axis] - 1;
    tmp->slice_idx = *slice_idx;
    if (vs->overlay_mode && tmp->n_levels > 1) load_all_levels(tmp, *var_idx);
    return 0;
}

static void video_free_plot(PlotfileData *tmp) {
    if (tmp->data) free(tmp->data);
    tmp->data = NULL;
    free_all_levels(tmp);
}

static void video_scan_slice(VideoState *vs, PlotfileData *tmp, int idx,
                             double *lo, double *hi) {
    int w = (vs->slice_axis == 2) ? tmp->grid_dims[0] : tmp->grid_dims[1 == vs->slice_axis ? 2 : 1];
    int h = (vs->slice_axis == 2) ? tmp->grid_dims[1] : tmp->grid_dims[2];
    if (vs->slice_axis == 1) { w = tmp->grid_dims[0]; h = tmp->grid_dims[2]; }
    double *s = (double *)malloc((size_t)w * h * sizeof(double));
    if (!s) return;
    extract_slice(tmp, s, vs->slice_axis, idx);
    unsigned char *mask = (unsigned char *)calloc((size_t)w * h, 1);
    if (!mask) { free(s); return; }
    for (int b = 0; b < tmp->n_boxes; b++) {
        Box *box = &tmp->boxes[b];
        int coord = idx + tmp->level_lo[vs->slice_axis];
        if (coord < box->lo[vs->slice_axis] || coord > box->hi[vs->slice_axis]) continue;
        int ax = (vs->slice_axis == 2) ? 0 : (vs->slice_axis == 1 ? 0 : 1);
        int ay = (vs->slice_axis == 2) ? 1 : 2;
        int x0 = box->lo[ax] - tmp->level_lo[ax], x1 = box->hi[ax] - tmp->level_lo[ax];
        int y0 = box->lo[ay] - tmp->level_lo[ay], y1 = box->hi[ay] - tmp->level_lo[ay];
        if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
        if (x1 >= w) x1 = w - 1; if (y1 >= h) y1 = h - 1;
        for (int y = y0; y <= y1; y++) for (int x = x0; x <= x1; x++) mask[y*w+x] = 1;
    }
    for (int p = 0; p < w*h; p++) if (mask[p] && video_value_valid(s[p], vs->scale_mode)) {
        double v = s[p];
        if (v < *lo) *lo = v; if (v > *hi) *hi = v;
    }
    free(mask); free(s);
}

static unsigned long video_pack_pixel(Visual *visual, RGB c) {
    unsigned long p = 0;
    unsigned long masks[3] = {visual->red_mask, visual->green_mask, visual->blue_mask};
    unsigned char vals[3] = {c.r, c.g, c.b};
    for (int n = 0; n < 3; n++) {
        if (!masks[n]) continue;
        int shift = 0; unsigned long m = masks[n];
        while (!(m & 1UL)) { shift++; m >>= 1; }
        p |= (((unsigned long)vals[n] * m + 127) / 255) << shift;
    }
    return p;
}

static unsigned char video_unpack(unsigned long p, unsigned long mask) {
    if (!mask) return 0;
    int shift = 0; unsigned long m = mask;
    while (!(m & 1UL)) { shift++; m >>= 1; }
    return (unsigned char)((((p & mask) >> shift) * 255UL + m/2) / m);
}

static int video_render_frame(VideoState *vs, PlotfileData *tmp, int idx, VideoFrame *out, int frame_no) {
    int sw = (vs->slice_axis == 2) ? tmp->grid_dims[0] : tmp->grid_dims[1];
    int sh = (vs->slice_axis == 2) ? tmp->grid_dims[1] : tmp->grid_dims[2];
    if (vs->slice_axis == 1) { sw = tmp->grid_dims[0]; sh = tmp->grid_dims[2]; }
    double *slice = (double *)malloc((size_t)sw * sh * sizeof(double));
    if (!slice) return -1;
    extract_slice(tmp, slice, vs->slice_axis, idx);
    Display *d = display; Visual *vis = DefaultVisual(d, screen);
    int depth = DefaultDepth(d, screen), bpp = (depth <= 16 ? 2 : 4);
    char *raw = (char *)calloc((size_t)VIDEO_FRAME_WIDTH * VIDEO_FRAME_HEIGHT, bpp);
    if (!raw) { free(slice); return -1; }
    XImage *im = XCreateImage(d, vis, depth, ZPixmap, 0, raw, VIDEO_FRAME_WIDTH,
                              VIDEO_FRAME_HEIGHT, 32, 0);
    if (!im) { free(raw); free(slice); return -1; }
    unsigned long white = video_pack_pixel(vis, (RGB){255,255,255});
    for (int y = 0; y < VIDEO_FRAME_HEIGHT; y++) for (int x = 0; x < VIDEO_FRAME_WIDTH; x++) XPutPixel(im,x,y,white);
    int ox=80, oy=65, ow=790, oh=570;
    for (int y=0; y<oh; y++) for (int x=0; x<ow; x++) {
        int sx = (int)((long long)x * sw / ow), sy = (int)((long long)(oh-1-y) * sh / oh);
        double v = slice[sy*sw+sx]; double z = 0.0;
        if (video_value_valid(v, vs->scale_mode) && vs->global_vmax > vs->global_vmin)
            z = (video_value_for_scale(v,vs->scale_mode)-video_value_for_scale(vs->global_vmin,vs->scale_mode)) /
                (video_value_for_scale(vs->global_vmax,vs->scale_mode)-video_value_for_scale(vs->global_vmin,vs->scale_mode));
        RGB c = video_value_valid(v,vs->scale_mode) ? get_colormap_rgb(z,vs->colormap) : (RGB){220,220,220};
        XPutPixel(im, ox+x, oy+y, video_pack_pixel(vis,c));
    }
    XPutImage(d, vs->frame_pixmap, gc, im, 0,0,0,0,VIDEO_FRAME_WIDTH,VIDEO_FRAME_HEIGHT);
    char title[256];
    snprintf(title,sizeof(title),"%s   t = %.8g   frame %d/%d",vs->variable_name,tmp->time,frame_no+1,vs->n_frames);
    XSetForeground(d,text_gc,BlackPixel(d,screen)); XDrawString(d,vs->frame_pixmap,text_gc,80,32,title,strlen(title));
    const char *unit = get_variable_unit(vs->variable_name); char cb[160];
    snprintf(cb,sizeof(cb),"%s %s",vs->variable_name,unit); XDrawString(d,vs->frame_pixmap,text_gc,885,55,cb,strlen(cb));
    for (int y=0;y<oh;y++) { RGB c=get_colormap_rgb(1.0-(double)y/oh,vs->colormap); XSetForeground(d,gc,video_pack_pixel(vis,c)); XFillRectangle(d,vs->frame_pixmap,gc,885,oy+y,35,1); }
    XDrawRectangle(d,vs->frame_pixmap,text_gc,ox,oy,ow,oh);
    const char *axes[] = {"X (m)", "Y (m)", "Z (m)"};
    int xa = vs->slice_axis == 2 ? 0 : (vs->slice_axis == 1 ? 0 : 1);
    int ya = 2;
    XDrawString(d,vs->frame_pixmap,text_gc,ox+ow/2-20,oy+oh+28,axes[xa],strlen(axes[xa]));
    XDrawString(d,vs->frame_pixmap,text_gc,18,oy+oh/2,axes[ya],strlen(axes[ya]));
    char range[160];
    snprintf(range,sizeof(range),"range: %.6g to %.6g (%s)",vs->global_vmin,vs->global_vmax,
             vs->scale_mode==1?"Log(+)":(vs->scale_mode==2?"Log(-)":"Linear"));
    XDrawString(d,vs->frame_pixmap,text_gc,80,VIDEO_FRAME_HEIGHT-18,range,strlen(range));
    XDestroyImage(im);
    XImage *got = XGetImage(d,vs->frame_pixmap,0,0,VIDEO_FRAME_WIDTH,VIDEO_FRAME_HEIGHT,AllPlanes,ZPixmap);
    if (!got) { free(slice); return -1; }
    out->rgb=(unsigned char*)malloc((size_t)VIDEO_FRAME_WIDTH*VIDEO_FRAME_HEIGHT*3);
    if (!out->rgb) { XDestroyImage(got); free(slice); return -1; }
    for (int y=0;y<VIDEO_FRAME_HEIGHT;y++) for(int x=0;x<VIDEO_FRAME_WIDTH;x++) { unsigned long p=XGetPixel(got,x,y); size_t q=((size_t)y*VIDEO_FRAME_WIDTH+x)*3; out->rgb[q]=video_unpack(p,vis->red_mask); out->rgb[q+1]=video_unpack(p,vis->green_mask); out->rgb[q+2]=video_unpack(p,vis->blue_mask); }
    out->time=tmp->time; out->timestep_number=timestep_numbers[vs->load_index];
    XDestroyImage(got); free(slice); return 0;
}

static void video_show_frame(VideoState *vs) {
    if (!vs || !vs->frames || !vs->frames[vs->current_frame].rgb) return;
    size_t n=(size_t)vs->frame_width*vs->frame_height*3;
    memcpy(vs->display_rgb,vs->frames[vs->current_frame].rgb,n);
    for (int y=0;y<vs->frame_height;y++) for(int x=0;x<vs->frame_width;x++) {
        size_t q=((size_t)y*vs->frame_width+x)*3;
        XPutPixel(vs->display_image,x,y,video_pack_pixel(DefaultVisual(display,screen),(RGB){vs->display_rgb[q],vs->display_rgb[q+1],vs->display_rgb[q+2]}));
    }
    XPutImage(display,vs->canvas,vs->gc,vs->display_image,0,0,0,0,vs->frame_width,vs->frame_height);
    char text[64]; snprintf(text,sizeof(text),"Frame %d/%d",vs->current_frame+1,vs->n_frames);
    XtVaSetValues(vs->frame_label,XtNlabel,text,NULL);
    float top=vs->n_frames>1?(float)vs->current_frame/(vs->n_frames-1):0.0f;
    XawScrollbarSetThumb(vs->scrubber,top,0.08f);
}

static void video_set_timer(VideoState *vs) {
    if (vs->timer_id) XtRemoveTimeOut(vs->timer_id);
    vs->timer_id=0;
    if (vs->playing) vs->timer_id=XtAppAddTimeOut(XtWidgetToApplicationContext(vs->shell),(unsigned long)(1000.0/vs->fps),video_timer_callback,vs);
}

static Boolean video_load_workproc(XtPointer client_data) {
    VideoState *vs=(VideoState*)client_data;
    if (!vs || vs->closing) return True;
    PlotfileData tmp; int var,idx;
    if (vs->load_index >= n_timesteps) {
        if (vs->load_phase==0) {
            if (!isfinite(vs->global_vmin) || !isfinite(vs->global_vmax)) { video_status(vs,"Error: no finite values for selected scale"); vs->loading=0; return True; }
            if (vs->global_vmin==vs->global_vmax) { double e=fmax(fabs(vs->global_vmin)*1e-6,1e-12); vs->global_vmin-=e; vs->global_vmax+=e; }
            size_t bytes=(size_t)n_timesteps*VIDEO_FRAME_WIDTH*VIDEO_FRAME_HEIGHT*3;
            fprintf(stderr,"Video frame cache: %zu bytes\n",bytes);
            if (bytes>VIDEO_CACHE_LIMIT) { video_status(vs,"Error: video cache exceeds 2 GiB safety limit"); vs->loading=0; return True; }
            vs->frames=(VideoFrame*)calloc(n_timesteps,sizeof(VideoFrame));
            if (!vs->frames) { video_status(vs,"Error: unable to allocate video frame cache"); vs->loading=0; return True; }
            vs->load_phase=1; vs->load_index=0;
        } else {
            vs->loading=0; XtSetSensitive(vs->play_button,True); XtSetSensitive(vs->stop_button,True); XtSetSensitive(vs->save_button,True);
            video_status(vs,"Ready"); video_show_frame(vs); return True;
        }
    }
    char status[128]; snprintf(status,sizeof(status),"Loading video: %s %d/%d",vs->load_phase==0?"scanning range":"rendering frames",vs->load_index+1,n_timesteps); video_status(vs,status);
    if (video_prepare_plot(vs,vs->load_index,&tmp,&var,&idx)<0) { video_free_plot(&tmp); video_status(vs,"Error: could not load a video timestep"); vs->loading=0; return True; }
    if (vs->load_phase==0) video_scan_slice(vs,&tmp,idx,&vs->global_vmin,&vs->global_vmax);
    else if (video_render_frame(vs,&tmp,idx,&vs->frames[vs->load_index],vs->load_index)<0) { video_free_plot(&tmp); video_status(vs,"Error: could not render video frame"); vs->loading=0; return True; }
    video_free_plot(&tmp); vs->load_index++;
    return False;
}

static void video_expose_callback(Widget w, XtPointer cd, XtPointer call_data) { (void)w;(void)call_data; video_show_frame((VideoState*)cd); }
static void video_play_callback(Widget w, XtPointer cd, XtPointer call_data) { (void)w;(void)call_data; VideoState *vs=cd; if(!vs||vs->loading)return; if(vs->current_frame>=vs->n_frames-1)vs->current_frame=0; vs->playing=1; video_show_frame(vs); video_set_timer(vs); }
static void video_stop_callback(Widget w, XtPointer cd, XtPointer call_data) { (void)w;(void)call_data; VideoState *vs=cd; if(vs){vs->playing=0;video_set_timer(vs);} }
static void video_scrub_callback(Widget w, XtPointer cd, XtPointer call_data) { (void)w; VideoState *vs=cd; float *p=(float*)call_data; if(vs&&p&&vs->n_frames>1){vs->current_frame=(int)floor(*p*(vs->n_frames-1)+0.5);if(vs->current_frame<0)vs->current_frame=0;if(vs->current_frame>=vs->n_frames)vs->current_frame=vs->n_frames-1;video_show_frame(vs);} }
static void video_timer_callback(XtPointer cd, XtIntervalId *id) { (void)id; VideoState *vs=cd; if(!vs||!vs->playing)return; if(vs->current_frame>=vs->n_frames-1){vs->playing=0;vs->timer_id=0;video_show_frame(vs);return;} vs->current_frame++; video_show_frame(vs); video_set_timer(vs); }

static void video_close_callback(Widget w, XtPointer cd, XtPointer call_data) {
    (void)w;(void)call_data; VideoState *vs=cd; if(!vs)return; vs->closing=1;
    if(vs->work_id)XtRemoveWorkProc(vs->work_id); if(vs->timer_id)XtRemoveTimeOut(vs->timer_id);
    if(vs->display_image)XDestroyImage(vs->display_image); if(vs->gc)XFreeGC(display,vs->gc); if(vs->frame_pixmap)XFreePixmap(display,vs->frame_pixmap); free(vs->display_rgb);
    if(vs->frames){for(int i=0;i<vs->n_frames;i++)free(vs->frames[i].rgb);free(vs->frames);} if(active_video==vs)active_video=NULL; XtDestroyWidget(vs->shell); free(vs);
}

static void video_export(VideoState *vs, const char *path, int mp4) {
    if(!vs||!path||!*path)return; char out[MAX_PATH]; snprintf(out,sizeof(out),"%s",path);
    const char *ext=mp4?".mp4":".gif"; if(!strstr(out,ext))strncat(out,ext,sizeof(out)-strlen(out)-1);
    if(access(out,F_OK)==0){video_status(vs,"Refusing to overwrite existing file");return;}
    char td[]="/tmp/pltview-video-XXXXXX"; if(!mkdtemp(td)){video_status(vs,"Error: cannot create temporary export directory");return;}
    char file[MAX_PATH]; int ok=1; for(int i=0;i<vs->n_frames;i++){snprintf(file,sizeof(file),"%s/frame_%06d.ppm",td,i+1);FILE*f=fopen(file,"wb");if(!f){ok=0;break;}fprintf(f,"P6\n%d %d\n255\n",vs->frame_width,vs->frame_height);fwrite(vs->frames[i].rgb,1,(size_t)vs->frame_width*vs->frame_height*3,f);fclose(f);}
    if(ok){char rate[32];snprintf(rate,sizeof(rate),"%.6g",vs->fps);char pattern[MAX_PATH];snprintf(pattern,sizeof(pattern),"%s/frame_%%06d.ppm",td);char *av[12];int n=0;av[n++]="ffmpeg";av[n++]="-y";av[n++]="-framerate";av[n++]=rate;av[n++]="-i";av[n++]=pattern;if(mp4){av[n++]="-pix_fmt";av[n++]="yuv420p";}else{av[n++]="-vf";av[n++]="format=rgb24";}av[n++]=out;av[n]=NULL;pid_t p=fork();if(p==0){execvp(av[0],av);_exit(127);}if(p>0&&waitpid(p,NULL,0)==p&&access(out,F_OK)==0)video_status(vs,"Video saved");else video_status(vs,"Error: ffmpeg export failed");}
    for(int i=0;i<vs->n_frames;i++){snprintf(file,sizeof(file),"%s/frame_%06d.ppm",td,i+1);unlink(file);}rmdir(td);
}
static void video_save_gif_callback(Widget w,XtPointer cd,XtPointer c){(void)w;(void)c;VideoState*vs=cd;video_export(vs,XawDialogGetValueString(vs->save_text),0);}
static void video_save_mp4_callback(Widget w,XtPointer cd,XtPointer c){(void)w;(void)c;VideoState*vs=cd;video_export(vs,XawDialogGetValueString(vs->save_text),1);}
static void video_save_callback(Widget w,XtPointer cd,XtPointer c){(void)w;(void)c;VideoState*vs=cd;if(vs&&vs->n_frames>0){vs->save_shell=XtVaCreatePopupShell("Save Video",transientShellWidgetClass,vs->shell,NULL);Widget d=XtVaCreateManagedWidget("Output path",dialogWidgetClass,vs->save_shell,XtNvalue,"animation",NULL);vs->save_text=d;Widget g=XtVaCreateManagedWidget("Save GIF",commandWidgetClass,d,XtNfromVert,d,NULL);XtAddCallback(g,XtNcallback,video_save_gif_callback,vs);Widget m=XtVaCreateManagedWidget("Save MP4",commandWidgetClass,d,XtNfromHoriz,g,NULL);XtAddCallback(m,XtNcallback,video_save_mp4_callback,vs);XtPopup(vs->save_shell,XtGrabNone);}}

static void video_fps_apply(Widget w,XtPointer cd,XtPointer c){(void)w;(void)c;VideoState*vs=cd;char *s=XawDialogGetValueString(vs->fps_text);char*e;double f=strtod(s,&e);if(e==s||*e||f<0.1||f>60.0){video_status(vs,"Invalid FPS (use 0.1 to 60.0)");return;}vs->fps=f;video_status(vs,"FPS updated");if(vs->playing)video_set_timer(vs);}

void video_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w;(void)client_data;(void)call_data; if(!global_pf||n_timesteps<2)return;
    if(active_video){XRaiseWindow(display,XtWindow(active_video->shell));return;}
    VideoState*vs=calloc(1,sizeof(*vs));if(!vs)return;active_video=vs;vs->n_frames=n_timesteps;vs->frame_width=VIDEO_FRAME_WIDTH;vs->frame_height=VIDEO_FRAME_HEIGHT;vs->fps=5.0;vs->slice_axis=global_pf->slice_axis;vs->requested_level=global_pf->current_level;vs->overlay_mode=global_pf->overlay_mode;vs->map_mode=global_pf->map_mode;vs->colormap=global_pf->colormap;vs->scale_mode=scale_mode;strncpy(vs->variable_name,global_pf->variables[global_pf->current_var],63);
    double span=global_pf->prob_hi[vs->slice_axis]-global_pf->prob_lo[vs->slice_axis];double cells=global_pf->level_lo[vs->slice_axis]+global_pf->grid_dims[vs->slice_axis];double dx=span/(cells>0?cells:1);vs->slice_position=global_pf->prob_lo[vs->slice_axis]+(global_pf->level_lo[vs->slice_axis]+global_pf->slice_idx+0.5)*dx;vs->global_vmin=INFINITY;vs->global_vmax=-INFINITY;
    vs->shell=XtVaCreatePopupShell("Video",topLevelShellWidgetClass,toplevel,XtNtitle,"PLTView Video",NULL);Widget f=XtVaCreateManagedWidget("videoForm",formWidgetClass,vs->shell,NULL);vs->status_label=XtVaCreateManagedWidget("Loading video",labelWidgetClass,f,XtNwidth,VIDEO_FRAME_WIDTH,NULL);vs->canvas_widget=XtVaCreateManagedWidget("videoCanvas",simpleWidgetClass,f,XtNfromVert,vs->status_label,XtNwidth,VIDEO_FRAME_WIDTH,XtNheight,VIDEO_FRAME_HEIGHT,NULL);vs->frame_label=XtVaCreateManagedWidget("Frame 0/0",labelWidgetClass,f,XtNfromVert,vs->canvas_widget,NULL);vs->scrubber=XtVaCreateManagedWidget("scrubber",scrollbarWidgetClass,f,XtNfromVert,vs->frame_label,XtNorientation,XtorientHorizontal,XtNwidth,700,NULL);XtAddCallback(vs->scrubber,XtNjumpProc,video_scrub_callback,vs);
    vs->play_button=XtVaCreateManagedWidget("Play",commandWidgetClass,f,XtNfromVert,vs->scrubber,NULL);vs->stop_button=XtVaCreateManagedWidget("Stop",commandWidgetClass,f,XtNfromVert,vs->scrubber,XtNfromHoriz,vs->play_button,NULL);vs->save_button=XtVaCreateManagedWidget("Save",commandWidgetClass,f,XtNfromVert,vs->scrubber,XtNfromHoriz,vs->stop_button,NULL);Widget close=XtVaCreateManagedWidget("Close",commandWidgetClass,f,XtNfromVert,vs->scrubber,XtNfromHoriz,vs->save_button,NULL);vs->fps_text=XtVaCreateManagedWidget("FPS",dialogWidgetClass,f,XtNfromVert,vs->play_button,XtNvalue,"5.0",NULL);Widget apply=XtVaCreateManagedWidget("Apply",commandWidgetClass,f,XtNfromVert,vs->play_button,XtNfromHoriz,vs->fps_text,NULL);XtAddCallback(apply,XtNcallback,video_fps_apply,vs);XtAddCallback(vs->play_button,XtNcallback,video_play_callback,vs);XtAddCallback(vs->stop_button,XtNcallback,video_stop_callback,vs);XtAddCallback(vs->save_button,XtNcallback,video_save_callback,vs);XtAddCallback(close,XtNcallback,video_close_callback,vs);XtAddCallback(vs->canvas_widget,XtNcallback,video_expose_callback,vs);XtRealizeWidget(vs->shell);vs->canvas=XtWindow(vs->canvas_widget);vs->gc=XCreateGC(display,vs->canvas,0,NULL);vs->display_rgb=calloc((size_t)VIDEO_FRAME_WIDTH*VIDEO_FRAME_HEIGHT*3,1);vs->display_image=XCreateImage(display,DefaultVisual(display,screen),DefaultDepth(display,screen),ZPixmap,0,(char*)calloc((size_t)VIDEO_FRAME_WIDTH*VIDEO_FRAME_HEIGHT,(DefaultDepth(display,screen)<=16?2:4)),VIDEO_FRAME_WIDTH,VIDEO_FRAME_HEIGHT,32,0);vs->frame_pixmap=XCreatePixmap(display,vs->canvas,VIDEO_FRAME_WIDTH,VIDEO_FRAME_HEIGHT,DefaultDepth(display,screen));XtPopup(vs->shell,XtGrabNone);XtSetSensitive(vs->play_button,False);XtSetSensitive(vs->stop_button,False);XtSetSensitive(vs->save_button,False);vs->loading=1;vs->work_id=XtAppAddWorkProc(XtWidgetToApplicationContext(vs->shell),video_load_workproc,vs);
}

void cleanup(PlotfileData *pf) {
    if (pf->data) free(pf->data);
    free_z_phys_cache(pf);
    free_all_levels(pf);
    if (pixel_data) free(pixel_data);
    pixel_data_size = 0;
    if (current_slice_data) free(current_slice_data);
    if (current_z_phys_slice) free(current_z_phys_slice);
    free(map_hover_lookup);
    map_hover_lookup = NULL;
    map_hover_lookup_size = 0;
    map_hover_lookup_active = 0;
}

/* ========== SDM Mode GUI and Rendering ========== */

static const char *sdm_metric_labels[SDM_N_METRICS] = {
    "Count", "SD Count", "Concentration", "Mass", "Mean Mult"
};

static const char *sdm_metric_ylabels[SDM_N_METRICS] = {
    "Particle count",
    "Super droplet count",
    "Concentration (#/m3)",
    "Mass (kg)",
    "Mean multiplicity"
};

static const char *sdm_metric_titles[SDM_N_METRICS] = {
    "Droplet Size Distribution - Particle Count",
    "Droplet Size Distribution - Super Droplet Count",
    "Droplet Size Distribution - Number Concentration",
    "Droplet Size Distribution - Mass",
    "Droplet Size Distribution - Mean Multiplicity"
};

/* Compute SDM histogram into the provided HistogramData struct */
void compute_sdm_histogram(ParticleData *pd, HistogramData *hist) {
    if (!pd || pd->n_particles <= 0 || !pd->radius) return;

    /* Convert radius to um and apply cutoff filter */
    int n_used = 0;
    double *radius_um = (double *)malloc(pd->n_particles * sizeof(double));
    double *mult_used = (double *)malloc(pd->n_particles * sizeof(double));
    double *mass_used = (double *)malloc(pd->n_particles * sizeof(double));

    for (int i = 0; i < pd->n_particles; i++) {
        double r_um = pd->radius[i] * 1e6;  /* Convert to um */
        if (pd->cutoff_radius > 0 && r_um <= pd->cutoff_radius) continue;
        radius_um[n_used] = r_um;
        mult_used[n_used] = pd->multiplicity[i];
        mass_used[n_used] = pd->mass[i];
        n_used++;
    }

    if (n_used == 0) {
        free(radius_um); free(mult_used); free(mass_used);
        /* Clear histogram */
        if (hist->bin_counts) { free(hist->bin_counts); hist->bin_counts = NULL; }
        if (hist->bin_centers) { free(hist->bin_centers); hist->bin_centers = NULL; }
        if (hist->bin_edges) { free(hist->bin_edges); hist->bin_edges = NULL; }
        hist->n_bins = 0;
        hist->count_max = 1;
        hist->mean = hist->std = hist->skewness = hist->kurtosis = 0;
        snprintf(hist->title, sizeof(hist->title), "%s", sdm_metric_titles[pd->current_metric]);
        snprintf(hist->xlabel, sizeof(hist->xlabel), "No particles after cutoff");
        return;
    }

    /* Find radius range in um */
    double rmin = radius_um[0], rmax = radius_um[0];
    for (int i = 1; i < n_used; i++) {
        if (radius_um[i] < rmin) rmin = radius_um[i];
        if (radius_um[i] > rmax) rmax = radius_um[i];
    }

    /* Determine number of bins */
    int n_bins;
    double bin_width;
    if (pd->custom_bin_width > 0) {
        bin_width = pd->custom_bin_width;
        n_bins = (int)ceil((rmax - rmin) / bin_width);
        if (n_bins < 1) n_bins = 1;
        if (n_bins > 500) n_bins = 500;
        /* Adjust rmax to fit whole bins */
        rmax = rmin + n_bins * bin_width;
    } else {
        n_bins = (int)(1 + 3.322 * log10((double)n_used));
        if (n_bins < 10) n_bins = 10;
        if (n_bins > 100) n_bins = 100;
        bin_width = (rmax - rmin) / n_bins;
        if (bin_width == 0) bin_width = 1.0;
    }

    /* Log-spaced bins if log_bin enabled and rmin > 0 */
    int use_log_bins = pd->log_bin && rmin > 0;
    double log_rmin = 0, log_rmax = 0, log_bin_width = 0;
    if (use_log_bins) {
        log_rmin = log10(rmin);
        log_rmax = log10(rmax);
        log_bin_width = (log_rmax - log_rmin) / n_bins;
        if (log_bin_width <= 0) log_bin_width = 1.0 / n_bins;
    }

    /* Allocate bin arrays */
    double *bin_counts = (double *)calloc(n_bins, sizeof(double));
    double *bin_centers = (double *)malloc(n_bins * sizeof(double));
    double *bin_edges = (double *)malloc((n_bins + 1) * sizeof(double));
    double *bin_sd_counts = (double *)calloc(n_bins, sizeof(double));
    double *bin_mass = (double *)calloc(n_bins, sizeof(double));

    if (use_log_bins) {
        for (int i = 0; i <= n_bins; i++) {
            bin_edges[i] = pow(10.0, log_rmin + i * log_bin_width);
        }
        for (int i = 0; i < n_bins; i++) {
            double log_center = log_rmin + (i + 0.5) * log_bin_width;
            bin_centers[i] = pow(10.0, log_center);
        }
    } else {
        for (int i = 0; i <= n_bins; i++) {
            bin_edges[i] = rmin + i * bin_width;
        }
        for (int i = 0; i < n_bins; i++) {
            bin_centers[i] = rmin + (i + 0.5) * bin_width;
        }
    }

    /* Accumulate per-bin values */
    for (int i = 0; i < n_used; i++) {
        int bin;
        if (use_log_bins) {
            double log_r = log10(radius_um[i]);
            bin = (int)((log_r - log_rmin) / log_bin_width);
        } else {
            bin = (int)((radius_um[i] - rmin) / bin_width);
        }
        if (bin < 0) bin = 0;
        if (bin >= n_bins) bin = n_bins - 1;

        bin_counts[bin] += mult_used[i];
        bin_sd_counts[bin] += 1.0;
        bin_mass[bin] += mass_used[i] * mult_used[i];
    }

    /* Select which metric to use as the displayed values */
    double *display_values = (double *)malloc(n_bins * sizeof(double));
    for (int i = 0; i < n_bins; i++) {
        switch (pd->current_metric) {
            case SDM_METRIC_PARTICLE_COUNT:
                display_values[i] = bin_counts[i];
                break;
            case SDM_METRIC_SD_COUNT:
                display_values[i] = bin_sd_counts[i];
                break;
            case SDM_METRIC_CONCENTRATION:
                display_values[i] = (pd->domain_volume > 0) ?
                    bin_counts[i] / pd->domain_volume : bin_counts[i];
                break;
            case SDM_METRIC_MASS:
                display_values[i] = bin_mass[i];
                break;
            case SDM_METRIC_MEAN_MULT:
                display_values[i] = (bin_sd_counts[i] > 0) ?
                    bin_counts[i] / bin_sd_counts[i] : 0.0;
                break;
            default:
                display_values[i] = bin_counts[i];
                break;
        }
    }

    /* Find max and total for scaling */
    double count_max = 0;
    double total_sum = 0;
    for (int i = 0; i < n_bins; i++) {
        if (display_values[i] > count_max) count_max = display_values[i];
        total_sum += display_values[i];
    }
    if (count_max == 0) count_max = 1;

    /* Compute statistics on radius in um (weighted by multiplicity) */
    double total_mult = 0, sum_r = 0, sum_r2 = 0;
    for (int i = 0; i < n_used; i++) {
        double w = mult_used[i];
        total_mult += w;
        sum_r += radius_um[i] * w;
        sum_r2 += radius_um[i] * radius_um[i] * w;
    }
    double mean = (total_mult > 0) ? sum_r / total_mult : 0;
    double variance = (total_mult > 0) ? (sum_r2 / total_mult) - (mean * mean) : 0;
    double std = (variance > 0) ? sqrt(variance) : 0;

    double sum_third = 0, sum_fourth = 0;
    for (int i = 0; i < n_used; i++) {
        double diff = radius_um[i] - mean;
        double d2 = diff * diff;
        sum_third += d2 * diff * mult_used[i];
        sum_fourth += d2 * d2 * mult_used[i];
    }
    double skewness = 0, kurtosis = 0;
    if (std > 0 && total_mult > 0) {
        double std3 = std * std * std;
        skewness = (sum_third / total_mult) / std3;
        kurtosis = (sum_fourth / total_mult) / (std * std * std * std) - 3.0;
    }

    /* Fill HistogramData */
    if (hist->bin_counts) free(hist->bin_counts);
    if (hist->bin_centers) free(hist->bin_centers);
    if (hist->bin_edges) { free(hist->bin_edges); hist->bin_edges = NULL; }

    hist->bin_counts = display_values;
    hist->bin_centers = bin_centers;
    hist->bin_edges = bin_edges;
    hist->n_bins = n_bins;
    hist->count_max = count_max;
    hist->total = total_sum;
    hist->bin_min = rmin;
    hist->bin_max = rmax;

    /* Apply xlim overrides if set */
    if (pd->xlim_min > 0) hist->bin_min = pd->xlim_min;
    if (pd->xlim_max > 0) hist->bin_max = pd->xlim_max;

    /* Apply PDF normalization if enabled: PDF = value / (total_sum * bin_width) */
    /* This makes the integral over radius equal to 1 */
    if (pd->pdf_mode) {
        double total_sum = 0;
        for (int i = 0; i < n_bins; i++) {
            total_sum += display_values[i];
        }
        if (total_sum > 0) {
            count_max = 0;
            for (int i = 0; i < n_bins; i++) {
                double bw = bin_edges[i + 1] - bin_edges[i];
                if (bw > 0) {
                    hist->bin_counts[i] = display_values[i] / (total_sum * bw);
                } else {
                    hist->bin_counts[i] = 0;
                }
                if (hist->bin_counts[i] > count_max) count_max = hist->bin_counts[i];
            }
            hist->count_max = (count_max > 0) ? count_max : 1.0;
        }
    }

    hist->mean = mean;
    hist->std = std;
    hist->skewness = skewness;
    hist->kurtosis = kurtosis;

    /* Set title */
    if (pd->pdf_mode) {
        snprintf(hist->title, sizeof(hist->title), "%s (PDF)", sdm_metric_titles[pd->current_metric]);
    } else {
        snprintf(hist->title, sizeof(hist->title), "%s", sdm_metric_titles[pd->current_metric]);
    }

    /* Use xlabel to carry cutoff info for second stats line */
    if (pd->cutoff_radius > 0) {
        snprintf(hist->xlabel, sizeof(hist->xlabel), "Cutoff: %.2f um, %d particles used",
                 pd->cutoff_radius, n_used);
    } else {
        hist->xlabel[0] = '\0';
    }

    free(bin_counts);
    free(bin_sd_counts);
    free(bin_mass);
    free(radius_um);
    free(mult_used);
    free(mass_used);
}

/* Render SDM histogram directly to the SDM canvas */
void render_sdm_histogram(ParticleData *pd) {
    if (!pd || !sdm_canvas || !display) return;

    /* Ensure histogram data exists */
    if (!sdm_hist_data) {
        sdm_hist_data = (HistogramData *)calloc(1, sizeof(HistogramData));
    }

    compute_sdm_histogram(pd, sdm_hist_data);

    /* Get canvas dimensions */
    Dimension width, height;
    XtVaGetValues(sdm_canvas_widget, XtNwidth, &width, XtNheight, &height, NULL);

    /* Use PDF label when in PDF mode */
    const char *ylabel;
    if (pd->pdf_mode) {
        ylabel = "Probability density (1/um)";
    } else {
        ylabel = sdm_metric_ylabels[pd->current_metric];
    }

    GC plot_gc = XCreateGC(display, sdm_canvas, 0, NULL);
    if (font) XSetFont(display, plot_gc, font->fid);
    draw_sdm_histogram(display, sdm_canvas, plot_gc, sdm_hist_data,
                       width, height, pd->log_x, pd->log_y, ylabel);
    XFreeGC(display, plot_gc);
}

/* Update SDM info label */
void update_sdm_info_label(ParticleData *pd, const char *plotfile_dir) {
    if (!sdm_info_label || !pd) return;
    char text[512];
    char total_str[64];
    const char *basename = strrchr(plotfile_dir, '/');
    basename = basename ? basename + 1 : plotfile_dir;

    const char *pdf_str = pd->pdf_mode ? " (PDF)" : "";

    /* Format total based on magnitude */
    double total = sdm_hist_data ? sdm_hist_data->total : 0;
    if (total >= 1e9) {
        snprintf(total_str, sizeof(total_str), "%.3e", total);
    } else if (total >= 1e6) {
        snprintf(total_str, sizeof(total_str), "%.2fM", total / 1e6);
    } else if (total >= 1e3) {
        snprintf(total_str, sizeof(total_str), "%.2fK", total / 1e3);
    } else if (total >= 1.0) {
        snprintf(total_str, sizeof(total_str), "%.2f", total);
    } else {
        snprintf(total_str, sizeof(total_str), "%.3e", total);
    }

    if (n_timesteps > 1) {
        snprintf(text, sizeof(text), "SDM: %s  |  Metric: %s%s  |  Total: %s  |  Step %d/%d",
                 basename, sdm_metric_labels[pd->current_metric], pdf_str, total_str,
                 current_timestep + 1, n_timesteps);
    } else {
        snprintf(text, sizeof(text), "SDM: %s  |  Metric: %s%s  |  Total: %s",
                 basename, sdm_metric_labels[pd->current_metric], pdf_str, total_str);
    }

    Arg args[1];
    XtSetArg(args[0], XtNlabel, text);
    XtSetValues(sdm_info_label, args, 1);
}

/* SDM metric button callback */
void sdm_metric_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int metric = (int)(long)client_data;
    if (global_pd && metric >= 0 && metric < SDM_N_METRICS) {
        global_pd->current_metric = metric;
        render_sdm_histogram(global_pd);
        update_sdm_info_label(global_pd, timestep_paths[current_timestep]);
    }
}

/* SDM timestep switch */
void sdm_switch_timestep(ParticleData *pd, int new_timestep) {
    if (new_timestep < 0 || new_timestep >= n_timesteps) return;
    current_timestep = new_timestep;

    /* Re-read particle data from new timestep */
    read_sdm_header(pd, timestep_paths[current_timestep]);
    pd->domain_volume = compute_domain_volume(timestep_paths[current_timestep]);
    read_sdm_data(pd, timestep_paths[current_timestep]);

    render_sdm_histogram(pd);
    update_sdm_info_label(pd, timestep_paths[current_timestep]);
    update_time_label();
}

/* SDM time navigation button callback */
void sdm_time_nav_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int dir = (int)(long)client_data;
    if (!global_pd || n_timesteps <= 1) return;

    int new_ts;
    if (dir == 0) {  /* prev */
        new_ts = current_timestep - 1;
        if (new_ts < 0) new_ts = n_timesteps - 1;
    } else {  /* next */
        new_ts = current_timestep + 1;
        if (new_ts >= n_timesteps) new_ts = 0;
    }
    sdm_switch_timestep(global_pd, new_ts);
}

/* SDM canvas expose handler */
void sdm_canvas_expose_callback(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != Expose) return;
    if (global_pd) {
        render_sdm_histogram(global_pd);
    }
}

/* SDM LogX toggle callback */
void sdm_logx_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pd) return;
    global_pd->log_x = !global_pd->log_x;
    render_sdm_histogram(global_pd);
}

/* SDM LogY toggle callback */
void sdm_logy_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pd) return;
    global_pd->log_y = !global_pd->log_y;
    render_sdm_histogram(global_pd);
}

/* SDM LogBin toggle callback - log-spaced bins (constant width in log space) */
void sdm_logbin_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pd) return;
    global_pd->log_bin = !global_pd->log_bin;
    render_sdm_histogram(global_pd);
}

/* SDM PDF toggle callback */
void sdm_pdf_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pd) return;
    global_pd->pdf_mode = !global_pd->pdf_mode;
    render_sdm_histogram(global_pd);
    update_sdm_info_label(global_pd, timestep_paths[current_timestep]);
}

/* SDM Settings apply callback */
void sdm_settings_apply_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pd) return;

    Arg args[1];
    String cutoff_str, binwidth_str, xlim_min_str, xlim_max_str;

    if (sdm_settings_text_cutoff) {
        XtSetArg(args[0], XtNstring, &cutoff_str);
        XtGetValues(sdm_settings_text_cutoff, args, 1);
        if (cutoff_str && strlen(cutoff_str) > 0)
            global_pd->cutoff_radius = atof(cutoff_str);
        else
            global_pd->cutoff_radius = 0;
    }

    if (sdm_settings_text_binwidth) {
        XtSetArg(args[0], XtNstring, &binwidth_str);
        XtGetValues(sdm_settings_text_binwidth, args, 1);
        if (binwidth_str && strlen(binwidth_str) > 0)
            global_pd->custom_bin_width = atof(binwidth_str);
        else
            global_pd->custom_bin_width = 0;
    }

    if (sdm_settings_text_xlim_min) {
        XtSetArg(args[0], XtNstring, &xlim_min_str);
        XtGetValues(sdm_settings_text_xlim_min, args, 1);
        if (xlim_min_str && strlen(xlim_min_str) > 0)
            global_pd->xlim_min = atof(xlim_min_str);
        else
            global_pd->xlim_min = 0;
    }

    if (sdm_settings_text_xlim_max) {
        XtSetArg(args[0], XtNstring, &xlim_max_str);
        XtGetValues(sdm_settings_text_xlim_max, args, 1);
        if (xlim_max_str && strlen(xlim_max_str) > 0)
            global_pd->xlim_max = atof(xlim_max_str);
        else
            global_pd->xlim_max = 0;
    }

    /* Close dialog */
    if (sdm_dialog_shell) {
        XtPopdown(sdm_dialog_shell);
        XtDestroyWidget(sdm_dialog_shell);
        sdm_dialog_shell = NULL;
    }
    sdm_dialog_active = 0;
    sdm_active_text_widget = NULL;
    sdm_settings_text_cutoff = NULL;
    sdm_settings_text_binwidth = NULL;
    sdm_settings_text_xlim_min = NULL;
    sdm_settings_text_xlim_max = NULL;

    render_sdm_histogram(global_pd);
    update_sdm_info_label(global_pd, timestep_paths[current_timestep]);
}

/* SDM Settings close callback */
void sdm_settings_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (sdm_dialog_shell) {
        XtPopdown(sdm_dialog_shell);
        XtDestroyWidget(sdm_dialog_shell);
        sdm_dialog_shell = NULL;
    }
    sdm_dialog_active = 0;
    sdm_active_text_widget = NULL;
    sdm_settings_text_cutoff = NULL;
    sdm_settings_text_binwidth = NULL;
    sdm_settings_text_xlim_min = NULL;
    sdm_settings_text_xlim_max = NULL;
}

/* SDM Settings cutoff focus callback */
void sdm_cutoff_focus_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    sdm_active_text_widget = sdm_settings_text_cutoff;
    sdm_active_field = 0;
}

/* SDM Settings binwidth focus callback */
void sdm_binwidth_focus_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    sdm_active_text_widget = sdm_settings_text_binwidth;
    sdm_active_field = 1;
}

/* SDM Settings text widget click-to-focus event handler */
void sdm_text_click_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != ButtonPress || !sdm_dialog_shell) {
        *continue_dispatch = True;
        return;
    }

    /* Set keyboard focus to clicked widget using the global dialog shell */
    XtSetKeyboardFocus(sdm_dialog_shell, w);
    Time time_val = CurrentTime;
    XtCallAcceptFocus(w, &time_val);

    /* Also set X input focus directly */
    XSetInputFocus(display, XtWindow(w), RevertToParent, CurrentTime);

    /* Update active widget tracking */
    sdm_active_text_widget = w;
    if (w == sdm_settings_text_cutoff) {
        sdm_active_field = 0;
    } else if (w == sdm_settings_text_binwidth) {
        sdm_active_field = 1;
    } else if (w == sdm_settings_text_xlim_min) {
        sdm_active_field = 2;
    } else if (w == sdm_settings_text_xlim_max) {
        sdm_active_field = 3;
    }

    *continue_dispatch = True;  /* Allow normal text widget processing */
}

/* SDM Settings button callback - open popup dialog */
void sdm_settings_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Arg args[10];
    int n;
    Widget form_w, label, button_w;
    char cutoff_str[64], binwidth_str[64], xlim_min_str[64], xlim_max_str[64];

    /* Format current values */
    if (global_pd && global_pd->cutoff_radius > 0)
        snprintf(cutoff_str, sizeof(cutoff_str), "%.4f", global_pd->cutoff_radius);
    else
        cutoff_str[0] = '\0';

    if (global_pd && global_pd->custom_bin_width > 0)
        snprintf(binwidth_str, sizeof(binwidth_str), "%.4f", global_pd->custom_bin_width);
    else
        binwidth_str[0] = '\0';

    if (global_pd && global_pd->xlim_min > 0)
        snprintf(xlim_min_str, sizeof(xlim_min_str), "%.4f", global_pd->xlim_min);
    else
        xlim_min_str[0] = '\0';

    if (global_pd && global_pd->xlim_max > 0)
        snprintf(xlim_max_str, sizeof(xlim_max_str), "%.4f", global_pd->xlim_max);
    else
        xlim_max_str[0] = '\0';

    n = 0;
    XtSetArg(args[n], XtNtitle, "SDM Settings"); n++;
    sdm_dialog_shell = XtCreatePopupShell("sdmSettings", transientShellWidgetClass, toplevel, args, n);

    n = 0;
    form_w = XtCreateManagedWidget("form", formWidgetClass, sdm_dialog_shell, args, n);

    /* Title label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Histogram settings:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    label = XtCreateManagedWidget("title", labelWidgetClass, form_w, args, n);

    /* Cutoff label */
    n = 0;
    XtSetArg(args[n], XtNfromVert, label); n++;
    XtSetArg(args[n], XtNlabel, "Cutoff (um):"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget cutoff_label = XtCreateManagedWidget("cutoffLabel", labelWidgetClass, form_w, args, n);

    /* Cutoff text input */
    n = 0;
    XtSetArg(args[n], XtNfromVert, label); n++;
    XtSetArg(args[n], XtNfromHoriz, cutoff_label); n++;
    XtSetArg(args[n], XtNwidth, 120); n++;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetArg(args[n], XtNstring, cutoff_str); n++;
    sdm_settings_text_cutoff = XtCreateManagedWidget("cutoffInput", asciiTextWidgetClass, form_w, args, n);
    XtAddEventHandler(sdm_settings_text_cutoff, ButtonPressMask, False, sdm_text_click_handler, NULL);

    /* Bin width label */
    n = 0;
    XtSetArg(args[n], XtNfromVert, cutoff_label); n++;
    XtSetArg(args[n], XtNlabel, "Bin width (um):"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget bw_label = XtCreateManagedWidget("bwLabel", labelWidgetClass, form_w, args, n);

    /* Bin width text input */
    n = 0;
    XtSetArg(args[n], XtNfromVert, cutoff_label); n++;
    XtSetArg(args[n], XtNfromHoriz, bw_label); n++;
    XtSetArg(args[n], XtNwidth, 120); n++;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetArg(args[n], XtNstring, binwidth_str); n++;
    sdm_settings_text_binwidth = XtCreateManagedWidget("bwInput", asciiTextWidgetClass, form_w, args, n);
    XtAddEventHandler(sdm_settings_text_binwidth, ButtonPressMask, False, sdm_text_click_handler, NULL);

    /* X-lim min label */
    n = 0;
    XtSetArg(args[n], XtNfromVert, bw_label); n++;
    XtSetArg(args[n], XtNlabel, "X-min (um):"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget xlim_min_label = XtCreateManagedWidget("xlimMinLabel", labelWidgetClass, form_w, args, n);

    /* X-lim min text input */
    n = 0;
    XtSetArg(args[n], XtNfromVert, bw_label); n++;
    XtSetArg(args[n], XtNfromHoriz, xlim_min_label); n++;
    XtSetArg(args[n], XtNwidth, 120); n++;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetArg(args[n], XtNstring, xlim_min_str); n++;
    sdm_settings_text_xlim_min = XtCreateManagedWidget("xlimMinInput", asciiTextWidgetClass, form_w, args, n);
    XtAddEventHandler(sdm_settings_text_xlim_min, ButtonPressMask, False, sdm_text_click_handler, NULL);

    /* X-lim max label */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_min_label); n++;
    XtSetArg(args[n], XtNlabel, "X-max (um):"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget xlim_max_label = XtCreateManagedWidget("xlimMaxLabel", labelWidgetClass, form_w, args, n);

    /* X-lim max text input */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_min_label); n++;
    XtSetArg(args[n], XtNfromHoriz, xlim_max_label); n++;
    XtSetArg(args[n], XtNwidth, 120); n++;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetArg(args[n], XtNstring, xlim_max_str); n++;
    sdm_settings_text_xlim_max = XtCreateManagedWidget("xlimMaxInput", asciiTextWidgetClass, form_w, args, n);
    XtAddEventHandler(sdm_settings_text_xlim_max, ButtonPressMask, False, sdm_text_click_handler, NULL);

    /* Apply button */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_max_label); n++;
    XtSetArg(args[n], XtNlabel, "Apply"); n++;
    button_w = XtCreateManagedWidget("apply", commandWidgetClass, form_w, args, n);
    XtAddCallback(button_w, XtNcallback, sdm_settings_apply_callback, NULL);

    /* Close button */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_max_label); n++;
    XtSetArg(args[n], XtNfromHoriz, button_w); n++;
    XtSetArg(args[n], XtNlabel, "Close"); n++;
    button_w = XtCreateManagedWidget("close", commandWidgetClass, form_w, args, n);
    XtAddCallback(button_w, XtNcallback, sdm_settings_close_callback, NULL);

    XtRealizeWidget(sdm_dialog_shell);
    XtPopup(sdm_dialog_shell, XtGrabNonexclusive);

    /* Set keyboard focus to cutoff text input */
    XtSetKeyboardFocus(sdm_dialog_shell, sdm_settings_text_cutoff);
    XSync(display, False);
    Time time_val = CurrentTime;
    XtCallAcceptFocus(sdm_settings_text_cutoff, &time_val);

    /* Also set X input focus to the text widget's window */
    XSetInputFocus(display, XtWindow(sdm_settings_text_cutoff), RevertToParent, CurrentTime);

    sdm_dialog_active = 1;
    sdm_active_text_widget = sdm_settings_text_cutoff;
    sdm_active_field = 0;
}

/* Initialize SDM GUI */
void init_sdm_gui(ParticleData *pd, const char *plotfile_dir, int argc, char **argv) {
    Arg args[20];
    int n;
    Widget button;

    global_pd = pd;

    toplevel = XtAppInitialize(NULL, "PLTView-SDM", NULL, 0, &argc, argv, NULL, NULL, 0);
    display = XtDisplay(toplevel);
    screen = DefaultScreen(display);

    /* Load font */
    font = XLoadQueryFont(display, "fixed");
    if (!font) font = XLoadQueryFont(display, "*");

    /* Main form */
    n = 0;
    XtSetArg(args[n], XtNwidth, 750); n++;
    XtSetArg(args[n], XtNheight, 600); n++;
    form = XtCreateManagedWidget("form", formWidgetClass, toplevel, args, n);

    /* Info label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "SDM - Loading..."); n++;
    XtSetArg(args[n], XtNwidth, 730); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    sdm_info_label = XtCreateManagedWidget("info", labelWidgetClass, form, args, n);

    /* Histogram canvas */
    n = 0;
    XtSetArg(args[n], XtNfromVert, sdm_info_label); n++;
    XtSetArg(args[n], XtNwidth, 700); n++;
    XtSetArg(args[n], XtNheight, 480); n++;
    XtSetArg(args[n], XtNborderWidth, 2); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    sdm_canvas_widget = XtCreateManagedWidget("sdmCanvas", simpleWidgetClass, form, args, n);

    /* Metric buttons row */
    Widget metric_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert, sdm_canvas_widget); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    metric_box = XtCreateManagedWidget("metricBox", boxWidgetClass, form, args, n);

    /* Metric label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Y-axis:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("metricLabel", labelWidgetClass, metric_box, args, n);

    for (int i = 0; i < SDM_N_METRICS; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, sdm_metric_labels[i]); n++;
        sdm_metric_buttons[i] = XtCreateManagedWidget(sdm_metric_labels[i],
            commandWidgetClass, metric_box, args, n);
        XtAddCallback(sdm_metric_buttons[i], XtNcallback, sdm_metric_callback, (XtPointer)(long)i);
    }

    /* Options row (LogX, LogY, Settings) */
    Widget options_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert, metric_box); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    options_box = XtCreateManagedWidget("optionsBox", boxWidgetClass, form, args, n);

    /* Options label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Options:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("optLabel", labelWidgetClass, options_box, args, n);

    /* LogX toggle button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "LogX"); n++;
    button = XtCreateManagedWidget("logX", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sdm_logx_callback, NULL);

    /* LogY toggle button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "LogY"); n++;
    button = XtCreateManagedWidget("logY", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sdm_logy_callback, NULL);

    /* LogBin toggle button - log-spaced bins */
    n = 0;
    XtSetArg(args[n], XtNlabel, "LogBin"); n++;
    button = XtCreateManagedWidget("logBin", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sdm_logbin_callback, NULL);

    /* PDF toggle button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "PDF"); n++;
    button = XtCreateManagedWidget("pdf", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sdm_pdf_callback, NULL);

    /* Settings button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Settings"); n++;
    button = XtCreateManagedWidget("settings", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sdm_settings_button_callback, NULL);

    /* Time navigation row (if multi-timestep) */
    if (n_timesteps > 1) {
        Widget time_box;
        n = 0;
        XtSetArg(args[n], XtNfromVert, options_box); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
        XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
        XtSetArg(args[n], XtNleft, XawChainLeft); n++;
        time_box = XtCreateManagedWidget("timeBox", boxWidgetClass, form, args, n);

        /* Time label */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Time"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        XtCreateManagedWidget("timeText", labelWidgetClass, time_box, args, n);

        /* < > buttons */
        const char *time_labels[] = {"<", ">"};
        for (int i = 0; i < 2; i++) {
            n = 0;
            XtSetArg(args[n], XtNlabel, time_labels[i]); n++;
            button = XtCreateManagedWidget(time_labels[i], commandWidgetClass, time_box, args, n);
            XtAddCallback(button, XtNcallback, sdm_time_nav_callback, (XtPointer)(long)i);
        }

        /* Time index display */
        char ts_text[32];
        snprintf(ts_text, sizeof(ts_text), "%d/%d", current_timestep + 1, n_timesteps);
        n = 0;
        XtSetArg(args[n], XtNlabel, ts_text); n++;
        XtSetArg(args[n], XtNwidth, 60); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        time_label = XtCreateManagedWidget("timeLabel", labelWidgetClass, time_box, args, n);
    }

    XtRealizeWidget(toplevel);

    /* Get canvas window */
    sdm_canvas = XtWindow(sdm_canvas_widget);

    /* Register expose handler */
    XtAddEventHandler(sdm_canvas_widget, ExposureMask, False, sdm_canvas_expose_callback, NULL);

    /* Enable keyboard events on canvas */
    XSelectInput(display, sdm_canvas, ExposureMask | KeyPressMask);
}

/* ========== SBM Mode GUI and Rendering ========== */

static const char *sbm_metric_labels[SBM_N_METRICS] = {
    "QC Mass", "QI Mass", "Total Mass", "QC Num", "QI Num", "Total Num"
};

static const char *sbm_metric_ylabels[SBM_N_METRICS] = {
    "Cloud water mass (kg/m3)",
    "Ice mass (kg/m3)",
    "Total mass (kg/m3)",
    "Cloud water number (#/m3)",
    "Ice number (#/m3)",
    "Total number (#/m3)"
};

static const char *sbm_metric_titles[SBM_N_METRICS] = {
    "SBM Droplet Size Distribution - Cloud Water Mass",
    "SBM Droplet Size Distribution - Ice Mass",
    "SBM Droplet Size Distribution - Total Mass (Cloud + Ice)",
    "SBM Droplet Size Distribution - Cloud Water Number",
    "SBM Droplet Size Distribution - Ice Number",
    "SBM Droplet Size Distribution - Total Number (Cloud + Ice)"
};

/* Read bin_info.txt to get SBM bin coordinates */
int read_sbm_bin_info(SBMData *sbm, const char *plotfile_dir) {
    char filepath[MAX_PATH];
    char line[MAX_LINE];

    snprintf(filepath, MAX_PATH, "%s/%s", plotfile_dir, SBM_BIN_INFO_FILE);

    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", filepath);
        return -1;
    }

    sbm->n_bins = 0;
    int in_hydrometeor_section = 0;
    int in_aerosol_section = 0;

    while (fgets(line, sizeof(line), fp) && sbm->n_bins < MAX_SBM_BINS) {
        /* Skip empty lines and comment-only lines */
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;

        if (*p == '\0') continue;

        /* Check for section headers */
        if (strstr(line, "Hydrometeor bin coordinates")) {
            in_hydrometeor_section = 1;
            in_aerosol_section = 0;
            continue;
        }
        if (strstr(line, "Aerosol bin coordinates")) {
            in_hydrometeor_section = 0;
            in_aerosol_section = 1;
            continue;
        }

        /* Skip comment lines */
        if (*p == '#') continue;

        /* Only parse hydrometeor section data */
        if (!in_hydrometeor_section || in_aerosol_section) continue;

        /* Parse: bin_index  diameter_um  radius_um  mass_kg ... */
        int bin_index;
        double diameter, radius;
        if (sscanf(line, "%d %lf %lf", &bin_index, &diameter, &radius) == 3) {
            if (sbm->n_bins < MAX_SBM_BINS) {
                sbm->bin_diameter_um[sbm->n_bins] = diameter;
                sbm->bin_radius_um[sbm->n_bins] = radius;
                sbm->n_bins++;
            }
        }
    }

    fclose(fp);

    if (sbm->n_bins == 0) {
        fprintf(stderr, "Error: No bins found in %s\n", filepath);
        return -1;
    }

    printf("SBM: Read %d bins from bin_info.txt\n", sbm->n_bins);
    return sbm->n_bins;
}

/* Compute SBM bin values by summing variable data across the domain */
int compute_sbm_values(SBMData *sbm, const char *plotfile_dir) {
    /* We need a temporary PlotfileData to read the plotfile */
    PlotfileData pf = {0};
    strncpy(pf.plotfile_dir, plotfile_dir, MAX_PATH - 1);

    if (read_header(&pf) < 0) {
        fprintf(stderr, "Error: Failed to read plotfile header\n");
        return -1;
    }

    if (read_cell_h(&pf) < 0) {
        fprintf(stderr, "Error: Failed to read cell header\n");
        return -1;
    }

    /* Compute domain volume */
    sbm->domain_volume = (pf.prob_hi[0] - pf.prob_lo[0]) *
                         (pf.prob_hi[1] - pf.prob_lo[1]) *
                         (pf.prob_hi[2] - pf.prob_lo[2]);

    /* Clear bin values */
    for (int i = 0; i < sbm->n_bins; i++) {
        sbm->bin_values[i] = 0.0;
    }

    /* Determine variable prefix based on metric */
    const char *prefix1 = NULL;
    const char *prefix2 = NULL;  /* For total metrics */

    switch (sbm->current_metric) {
        case SBM_METRIC_QC_MASS:
            prefix1 = "sbm_qc_bin_mass_";
            break;
        case SBM_METRIC_QI_MASS:
            prefix1 = "sbm_qi_bin_mass_";
            break;
        case SBM_METRIC_TOTAL_MASS:
            prefix1 = "sbm_qc_bin_mass_";
            prefix2 = "sbm_qi_bin_mass_";
            break;
        case SBM_METRIC_QC_NUM:
            prefix1 = "sbm_qc_bin_num_";
            break;
        case SBM_METRIC_QI_NUM:
            prefix1 = "sbm_qi_bin_num_";
            break;
        case SBM_METRIC_TOTAL_NUM:
            prefix1 = "sbm_qc_bin_num_";
            prefix2 = "sbm_qi_bin_num_";
            break;
        default:
            prefix1 = "sbm_qc_bin_mass_";
            break;
    }

    /* Read and sum each bin */
    for (int bin = 0; bin < sbm->n_bins; bin++) {
        char varname[64];
        double total = 0.0;

        /* First component */
        snprintf(varname, sizeof(varname), "%s%02d", prefix1, bin);
        int var_idx = find_variable_index(&pf, varname);

        if (var_idx >= 0) {
            if (read_variable_data(&pf, var_idx) == 0) {
                /* Sum all grid cells */
                int total_size = pf.grid_dims[0] * pf.grid_dims[1] * pf.grid_dims[2];
                for (int i = 0; i < total_size; i++) {
                    total += pf.data[i];
                }
            }
        }

        /* Second component (for total metrics) */
        if (prefix2) {
            snprintf(varname, sizeof(varname), "%s%02d", prefix2, bin);
            var_idx = find_variable_index(&pf, varname);

            if (var_idx >= 0) {
                if (read_variable_data(&pf, var_idx) == 0) {
                    int total_size = pf.grid_dims[0] * pf.grid_dims[1] * pf.grid_dims[2];
                    for (int i = 0; i < total_size; i++) {
                        total += pf.data[i];
                    }
                }
            }
        }

        sbm->bin_values[bin] = total;
    }

    /* Cleanup temporary plotfile data */
    if (pf.data) {
        free(pf.data);
        pf.data = NULL;
    }

    return 0;
}

/* Compute SBM histogram into HistogramData struct */
void compute_sbm_histogram(SBMData *sbm, HistogramData *hist) {
    if (!sbm || sbm->n_bins <= 0) return;

    /* Allocate/reallocate histogram arrays */
    if (hist->bin_counts) free(hist->bin_counts);
    if (hist->bin_centers) free(hist->bin_centers);
    if (hist->bin_edges) free(hist->bin_edges);

    hist->bin_counts = (double *)malloc(sbm->n_bins * sizeof(double));
    hist->bin_centers = (double *)malloc(sbm->n_bins * sizeof(double));
    hist->bin_edges = (double *)malloc((sbm->n_bins + 1) * sizeof(double));
    hist->n_bins = sbm->n_bins;

    /* Copy bin values and centers, compute max and total */
    double count_max = 0;
    double total_sum = 0;
    for (int i = 0; i < sbm->n_bins; i++) {
        hist->bin_counts[i] = sbm->bin_values[i];
        hist->bin_centers[i] = sbm->bin_radius_um[i];
        if (sbm->bin_values[i] > count_max) count_max = sbm->bin_values[i];
        total_sum += sbm->bin_values[i];
    }
    hist->count_max = (count_max > 0) ? count_max : 1.0;
    hist->total = total_sum;

    /* Compute bin edges using geometric mean: edge[i] = sqrt(center[i-1] * center[i]) */
    /* For log-spaced bins, this gives proper boundaries */
    for (int i = 1; i < sbm->n_bins; i++) {
        hist->bin_edges[i] = sqrt(sbm->bin_radius_um[i-1] * sbm->bin_radius_um[i]);
    }
    /* Extrapolate first and last edges using the same log ratio */
    if (sbm->n_bins >= 2) {
        double ratio = sbm->bin_radius_um[1] / sbm->bin_radius_um[0];
        hist->bin_edges[0] = sbm->bin_radius_um[0] / sqrt(ratio);
        hist->bin_edges[sbm->n_bins] = sbm->bin_radius_um[sbm->n_bins-1] * sqrt(ratio);
    } else {
        hist->bin_edges[0] = sbm->bin_radius_um[0] * 0.5;
        hist->bin_edges[1] = sbm->bin_radius_um[0] * 2.0;
    }

    /* Set range from bin edges */
    hist->bin_min = hist->bin_edges[0];
    hist->bin_max = hist->bin_edges[sbm->n_bins];

    /* Apply xlim overrides if set */
    if (sbm->xlim_min > 0) hist->bin_min = sbm->xlim_min;
    if (sbm->xlim_max > 0) hist->bin_max = sbm->xlim_max;

    /* Compute statistics (weighted by bin values) - use original values */
    double total_weight = 0, sum_r = 0, sum_r2 = 0;
    for (int i = 0; i < sbm->n_bins; i++) {
        double w = sbm->bin_values[i];
        if (w > 0) {
            double r = sbm->bin_radius_um[i];
            total_weight += w;
            sum_r += r * w;
            sum_r2 += r * r * w;
        }
    }

    double mean = (total_weight > 0) ? sum_r / total_weight : 0;
    double variance = (total_weight > 0) ? (sum_r2 / total_weight) - (mean * mean) : 0;
    double std = (variance > 0) ? sqrt(variance) : 0;

    double sum_third = 0, sum_fourth = 0;
    for (int i = 0; i < sbm->n_bins; i++) {
        double w = sbm->bin_values[i];
        if (w > 0) {
            double diff = sbm->bin_radius_um[i] - mean;
            double d2 = diff * diff;
            sum_third += d2 * diff * w;
            sum_fourth += d2 * d2 * w;
        }
    }

    double skewness = 0, kurtosis = 0;
    if (std > 0 && total_weight > 0) {
        double std3 = std * std * std;
        skewness = (sum_third / total_weight) / std3;
        kurtosis = (sum_fourth / total_weight) / (std * std * std * std) - 3.0;
    }

    hist->mean = mean;
    hist->std = std;
    hist->skewness = skewness;
    hist->kurtosis = kurtosis;

    /* Apply PDF normalization if enabled: PDF = value / (total_sum * bin_width) */
    /* This makes the integral over radius equal to 1 */
    if (sbm->pdf_mode && total_weight > 0) {
        count_max = 0;
        for (int i = 0; i < sbm->n_bins; i++) {
            double bin_width = hist->bin_edges[i + 1] - hist->bin_edges[i];
            if (bin_width > 0) {
                hist->bin_counts[i] = sbm->bin_values[i] / (total_weight * bin_width);
            } else {
                hist->bin_counts[i] = 0;
            }
            if (hist->bin_counts[i] > count_max) count_max = hist->bin_counts[i];
        }
        hist->count_max = (count_max > 0) ? count_max : 1.0;
    }

    /* Set title */
    if (sbm->pdf_mode) {
        snprintf(hist->title, sizeof(hist->title), "%s (PDF)", sbm_metric_titles[sbm->current_metric]);
    } else {
        snprintf(hist->title, sizeof(hist->title), "%s", sbm_metric_titles[sbm->current_metric]);
    }
    hist->xlabel[0] = '\0';  /* No cutoff info for SBM */
}

/* Render SBM histogram to the SBM canvas */
void render_sbm_histogram(SBMData *sbm) {
    if (!sbm || !sbm_canvas || !display) return;

    /* Ensure histogram data exists */
    if (!sbm_hist_data) {
        sbm_hist_data = (HistogramData *)calloc(1, sizeof(HistogramData));
    }

    compute_sbm_histogram(sbm, sbm_hist_data);

    /* Get canvas dimensions */
    Dimension width, height;
    XtVaGetValues(sbm_canvas_widget, XtNwidth, &width, XtNheight, &height, NULL);

    /* Use PDF label when in PDF mode */
    const char *ylabel;
    if (sbm->pdf_mode) {
        ylabel = "Probability density (1/um)";
    } else {
        ylabel = sbm_metric_ylabels[sbm->current_metric];
    }

    GC plot_gc = XCreateGC(display, sbm_canvas, 0, NULL);
    if (font) XSetFont(display, plot_gc, font->fid);
    draw_sdm_histogram(display, sbm_canvas, plot_gc, sbm_hist_data,
                       width, height, sbm->log_x, sbm->log_y, ylabel);
    XFreeGC(display, plot_gc);
}

/* Update SBM info label */
void update_sbm_info_label(SBMData *sbm, const char *plotfile_dir) {
    if (!sbm_info_label || !sbm) return;
    char text[512];
    char total_str[64];
    const char *basename = strrchr(plotfile_dir, '/');
    basename = basename ? basename + 1 : plotfile_dir;

    const char *pdf_str = sbm->pdf_mode ? " (PDF)" : "";

    /* Format total based on magnitude */
    double total = sbm_hist_data ? sbm_hist_data->total : 0;
    if (total >= 1e9) {
        snprintf(total_str, sizeof(total_str), "%.3e", total);
    } else if (total >= 1e6) {
        snprintf(total_str, sizeof(total_str), "%.2fM", total / 1e6);
    } else if (total >= 1e3) {
        snprintf(total_str, sizeof(total_str), "%.2fK", total / 1e3);
    } else if (total >= 1.0) {
        snprintf(total_str, sizeof(total_str), "%.2f", total);
    } else {
        snprintf(total_str, sizeof(total_str), "%.3e", total);
    }

    if (n_timesteps > 1) {
        snprintf(text, sizeof(text), "SBM: %s  |  Metric: %s%s  |  Total: %s  |  Step %d/%d",
                 basename, sbm_metric_labels[sbm->current_metric], pdf_str, total_str,
                 current_timestep + 1, n_timesteps);
    } else {
        snprintf(text, sizeof(text), "SBM: %s  |  Metric: %s%s  |  Total: %s",
                 basename, sbm_metric_labels[sbm->current_metric], pdf_str, total_str);
    }

    Arg args[1];
    XtSetArg(args[0], XtNlabel, text);
    XtSetValues(sbm_info_label, args, 1);
}

/* SBM metric button callback */
void sbm_metric_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int metric = (int)(long)client_data;
    if (global_sbm && metric >= 0 && metric < SBM_N_METRICS) {
        global_sbm->current_metric = metric;
        compute_sbm_values(global_sbm, timestep_paths[current_timestep]);
        render_sbm_histogram(global_sbm);
        update_sbm_info_label(global_sbm, timestep_paths[current_timestep]);
    }
}

/* SBM LogX toggle callback */
void sbm_logx_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_sbm) return;
    global_sbm->log_x = !global_sbm->log_x;
    render_sbm_histogram(global_sbm);
}

/* SBM LogY toggle callback */
void sbm_logy_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_sbm) return;
    global_sbm->log_y = !global_sbm->log_y;
    render_sbm_histogram(global_sbm);
}

/* SBM PDF toggle callback */
void sbm_pdf_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_sbm) return;
    global_sbm->pdf_mode = !global_sbm->pdf_mode;
    render_sbm_histogram(global_sbm);
    update_sbm_info_label(global_sbm, timestep_paths[current_timestep]);
}

/* SBM Settings apply callback */
void sbm_settings_apply_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Arg args[1];
    String xlim_min_str, xlim_max_str;

    if (sbm_settings_text_xlim_min) {
        XtSetArg(args[0], XtNstring, &xlim_min_str);
        XtGetValues(sbm_settings_text_xlim_min, args, 1);
        if (xlim_min_str && strlen(xlim_min_str) > 0)
            global_sbm->xlim_min = atof(xlim_min_str);
        else
            global_sbm->xlim_min = 0;
    }

    if (sbm_settings_text_xlim_max) {
        XtSetArg(args[0], XtNstring, &xlim_max_str);
        XtGetValues(sbm_settings_text_xlim_max, args, 1);
        if (xlim_max_str && strlen(xlim_max_str) > 0)
            global_sbm->xlim_max = atof(xlim_max_str);
        else
            global_sbm->xlim_max = 0;
    }

    /* Re-render histogram with new settings */
    render_sbm_histogram(global_sbm);
}

/* SBM Settings close callback */
void sbm_settings_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (sbm_dialog_shell) {
        XtPopdown(sbm_dialog_shell);
        XtDestroyWidget(sbm_dialog_shell);
        sbm_dialog_shell = NULL;
    }
    sbm_dialog_active = 0;
    sbm_active_text_widget = NULL;
    sbm_settings_text_xlim_min = NULL;
    sbm_settings_text_xlim_max = NULL;
}

/* SBM Settings text widget click-to-focus event handler */
void sbm_text_click_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != ButtonPress || !sbm_dialog_shell) {
        *continue_dispatch = True;
        return;
    }

    /* Set keyboard focus to clicked widget using the global dialog shell */
    XtSetKeyboardFocus(sbm_dialog_shell, w);
    Time time_val = CurrentTime;
    XtCallAcceptFocus(w, &time_val);

    /* Also set X input focus directly */
    XSetInputFocus(display, XtWindow(w), RevertToParent, CurrentTime);

    /* Update active widget tracking */
    sbm_active_text_widget = w;
    if (w == sbm_settings_text_xlim_min) {
        sbm_active_field = 0;
    } else if (w == sbm_settings_text_xlim_max) {
        sbm_active_field = 1;
    }

    *continue_dispatch = True;  /* Allow normal text widget processing */
}

/* SBM Settings button callback - open popup dialog */
void sbm_settings_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Arg args[10];
    int n;
    Widget form_w, label, button_w;
    char xlim_min_str[64], xlim_max_str[64];

    /* Format current values */
    if (global_sbm && global_sbm->xlim_min > 0)
        snprintf(xlim_min_str, sizeof(xlim_min_str), "%.4f", global_sbm->xlim_min);
    else
        xlim_min_str[0] = '\0';

    if (global_sbm && global_sbm->xlim_max > 0)
        snprintf(xlim_max_str, sizeof(xlim_max_str), "%.4f", global_sbm->xlim_max);
    else
        xlim_max_str[0] = '\0';

    n = 0;
    XtSetArg(args[n], XtNtitle, "SBM Settings"); n++;
    sbm_dialog_shell = XtCreatePopupShell("sbmSettings", transientShellWidgetClass, toplevel, args, n);

    n = 0;
    form_w = XtCreateManagedWidget("form", formWidgetClass, sbm_dialog_shell, args, n);

    /* Title label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Histogram settings:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    label = XtCreateManagedWidget("title", labelWidgetClass, form_w, args, n);

    /* X-lim min label */
    n = 0;
    XtSetArg(args[n], XtNfromVert, label); n++;
    XtSetArg(args[n], XtNlabel, "X-min (um):"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget xlim_min_label = XtCreateManagedWidget("xlimMinLabel", labelWidgetClass, form_w, args, n);

    /* X-lim min text input */
    n = 0;
    XtSetArg(args[n], XtNfromVert, label); n++;
    XtSetArg(args[n], XtNfromHoriz, xlim_min_label); n++;
    XtSetArg(args[n], XtNwidth, 120); n++;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetArg(args[n], XtNstring, xlim_min_str); n++;
    sbm_settings_text_xlim_min = XtCreateManagedWidget("xlimMinInput", asciiTextWidgetClass, form_w, args, n);
    XtAddEventHandler(sbm_settings_text_xlim_min, ButtonPressMask, False, sbm_text_click_handler, NULL);

    /* X-lim max label */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_min_label); n++;
    XtSetArg(args[n], XtNlabel, "X-max (um):"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    Widget xlim_max_label = XtCreateManagedWidget("xlimMaxLabel", labelWidgetClass, form_w, args, n);

    /* X-lim max text input */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_min_label); n++;
    XtSetArg(args[n], XtNfromHoriz, xlim_max_label); n++;
    XtSetArg(args[n], XtNwidth, 120); n++;
    XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
    XtSetArg(args[n], XtNstring, xlim_max_str); n++;
    sbm_settings_text_xlim_max = XtCreateManagedWidget("xlimMaxInput", asciiTextWidgetClass, form_w, args, n);
    XtAddEventHandler(sbm_settings_text_xlim_max, ButtonPressMask, False, sbm_text_click_handler, NULL);

    /* Apply button */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_max_label); n++;
    XtSetArg(args[n], XtNlabel, "Apply"); n++;
    button_w = XtCreateManagedWidget("apply", commandWidgetClass, form_w, args, n);
    XtAddCallback(button_w, XtNcallback, sbm_settings_apply_callback, NULL);

    /* Close button */
    n = 0;
    XtSetArg(args[n], XtNfromVert, xlim_max_label); n++;
    XtSetArg(args[n], XtNfromHoriz, button_w); n++;
    XtSetArg(args[n], XtNlabel, "Close"); n++;
    button_w = XtCreateManagedWidget("close", commandWidgetClass, form_w, args, n);
    XtAddCallback(button_w, XtNcallback, sbm_settings_close_callback, NULL);

    XtRealizeWidget(sbm_dialog_shell);
    XtPopup(sbm_dialog_shell, XtGrabNonexclusive);

    /* Set keyboard focus to xlim_min text input */
    XtSetKeyboardFocus(sbm_dialog_shell, sbm_settings_text_xlim_min);
    XSync(display, False);
    Time time_val = CurrentTime;
    XtCallAcceptFocus(sbm_settings_text_xlim_min, &time_val);

    /* Also set X input focus to the text widget's window */
    XSetInputFocus(display, XtWindow(sbm_settings_text_xlim_min), RevertToParent, CurrentTime);

    sbm_dialog_active = 1;
    sbm_active_text_widget = sbm_settings_text_xlim_min;
    sbm_active_field = 0;
}

/* SBM timestep switch */
void sbm_switch_timestep(SBMData *sbm, int new_timestep) {
    if (new_timestep < 0 || new_timestep >= n_timesteps) return;
    current_timestep = new_timestep;

    /* Re-read bin info from new timestep (in case bins differ) */
    read_sbm_bin_info(sbm, timestep_paths[current_timestep]);

    /* Recompute values from new timestep */
    compute_sbm_values(sbm, timestep_paths[current_timestep]);

    render_sbm_histogram(sbm);
    update_sbm_info_label(sbm, timestep_paths[current_timestep]);
    update_time_label();
}

/* SBM time navigation button callback */
void sbm_time_nav_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int dir = (int)(long)client_data;
    if (!global_sbm || n_timesteps <= 1) return;

    int new_ts;
    if (dir == 0) {  /* prev */
        new_ts = current_timestep - 1;
        if (new_ts < 0) new_ts = n_timesteps - 1;
    } else {  /* next */
        new_ts = current_timestep + 1;
        if (new_ts >= n_timesteps) new_ts = 0;
    }
    sbm_switch_timestep(global_sbm, new_ts);
}

/* SBM canvas expose handler */
void sbm_canvas_expose_callback(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (event->type != Expose) return;
    if (global_sbm) {
        render_sbm_histogram(global_sbm);
    }
}

/* Initialize SBM GUI */
void init_sbm_gui(SBMData *sbm, const char *plotfile_dir, int argc, char **argv) {
    Arg args[20];
    int n;
    Widget button;

    global_sbm = sbm;

    toplevel = XtAppInitialize(NULL, "PLTView-SBM", NULL, 0, &argc, argv, NULL, NULL, 0);
    display = XtDisplay(toplevel);
    screen = DefaultScreen(display);

    /* Load font */
    font = XLoadQueryFont(display, "fixed");
    if (!font) font = XLoadQueryFont(display, "*");

    /* Main form */
    n = 0;
    XtSetArg(args[n], XtNwidth, 750); n++;
    XtSetArg(args[n], XtNheight, 600); n++;
    form = XtCreateManagedWidget("form", formWidgetClass, toplevel, args, n);

    /* Info label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "SBM - Loading..."); n++;
    XtSetArg(args[n], XtNwidth, 730); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    sbm_info_label = XtCreateManagedWidget("info", labelWidgetClass, form, args, n);

    /* Histogram canvas */
    n = 0;
    XtSetArg(args[n], XtNfromVert, sbm_info_label); n++;
    XtSetArg(args[n], XtNwidth, 700); n++;
    XtSetArg(args[n], XtNheight, 480); n++;
    XtSetArg(args[n], XtNborderWidth, 2); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    sbm_canvas_widget = XtCreateManagedWidget("sbmCanvas", simpleWidgetClass, form, args, n);

    /* Metric buttons row */
    Widget metric_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert, sbm_canvas_widget); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    metric_box = XtCreateManagedWidget("metricBox", boxWidgetClass, form, args, n);

    /* Metric label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Y-axis:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("metricLabel", labelWidgetClass, metric_box, args, n);

    for (int i = 0; i < SBM_N_METRICS; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, sbm_metric_labels[i]); n++;
        sbm_metric_buttons[i] = XtCreateManagedWidget(sbm_metric_labels[i],
            commandWidgetClass, metric_box, args, n);
        XtAddCallback(sbm_metric_buttons[i], XtNcallback, sbm_metric_callback, (XtPointer)(long)i);
    }

    /* Options row (LogX, LogY, PDF) */
    Widget options_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert, metric_box); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    options_box = XtCreateManagedWidget("optionsBox", boxWidgetClass, form, args, n);

    /* Options label */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Options:"); n++;
    XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("optLabel", labelWidgetClass, options_box, args, n);

    /* LogX toggle button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "LogX"); n++;
    button = XtCreateManagedWidget("logX", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sbm_logx_callback, NULL);

    /* LogY toggle button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "LogY"); n++;
    button = XtCreateManagedWidget("logY", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sbm_logy_callback, NULL);

    /* PDF toggle button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "PDF"); n++;
    button = XtCreateManagedWidget("pdf", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sbm_pdf_callback, NULL);

    /* Settings button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Settings"); n++;
    button = XtCreateManagedWidget("settings", commandWidgetClass, options_box, args, n);
    XtAddCallback(button, XtNcallback, sbm_settings_button_callback, NULL);

    /* Time navigation row (if multi-timestep) */
    if (n_timesteps > 1) {
        Widget time_box;
        n = 0;
        XtSetArg(args[n], XtNfromVert, options_box); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
        XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
        XtSetArg(args[n], XtNleft, XawChainLeft); n++;
        time_box = XtCreateManagedWidget("timeBox", boxWidgetClass, form, args, n);

        /* Time label */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Time"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        XtCreateManagedWidget("timeText", labelWidgetClass, time_box, args, n);

        /* < > buttons */
        const char *time_labels[] = {"<", ">"};
        for (int i = 0; i < 2; i++) {
            n = 0;
            XtSetArg(args[n], XtNlabel, time_labels[i]); n++;
            button = XtCreateManagedWidget(time_labels[i], commandWidgetClass, time_box, args, n);
            XtAddCallback(button, XtNcallback, sbm_time_nav_callback, (XtPointer)(long)i);
        }

        /* Time index display */
        char ts_text[32];
        snprintf(ts_text, sizeof(ts_text), "%d/%d", current_timestep + 1, n_timesteps);
        n = 0;
        XtSetArg(args[n], XtNlabel, ts_text); n++;
        XtSetArg(args[n], XtNwidth, 60); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        time_label = XtCreateManagedWidget("timeLabel", labelWidgetClass, time_box, args, n);
    }

    XtRealizeWidget(toplevel);

    /* Get canvas window */
    sbm_canvas = XtWindow(sbm_canvas_widget);

    /* Register expose handler */
    XtAddEventHandler(sbm_canvas_widget, ExposureMask, False, sbm_canvas_expose_callback, NULL);

    /* Enable keyboard events on canvas */
    XSelectInput(display, sbm_canvas, ExposureMask | KeyPressMask);
}

/* ============================================================
 * 1D Profile viewer implementation
 * ============================================================ */

/* Column name tables for each file type */
static const char *profile_col_names_surf[]    = {"time","u_star","t_star","olen",NULL};
static const char *profile_col_names_mean[]    = {"time","z","u","v","w","rho","theta","ksgs","Kmv","Khv","qv","qc","qr","qi","qs","qg",NULL};
static const char *profile_col_names_flux[]    = {"time","z","u'u'","u'v'","u'w'","v'v'","v'w'","w'w'","u'th'","v'th'","w'th'","th'th'","uiuiu'","uiuiv'","uiuiw'","p'u'","p'v'","p'w'","w'qv'","w'qc'","w'qr'","w'thv'",NULL};
static const char *profile_col_names_subgrid[] = {"time","z","tau11","tau12","tau13","tau22","tau23","tau33","sgshfx","sgsq1fx","sgsq2fx","sgsdiss",NULL};
static const char *profile_file_labels[N_PROFILE_FILES] = {"surf","mean","flux","subgrid"};

int read_profile_file(ProfileFile *pf, const char *path, int file_type) {
    FILE *fp = fopen(path, "r");
    if (!fp) return -1;

    /* Set column names from table */
    const char **names = NULL;
    switch (file_type) {
        case PROFILE_FILE_SURF:    names = profile_col_names_surf;    pf->has_z = 0; break;
        case PROFILE_FILE_MEAN:    names = profile_col_names_mean;    pf->has_z = 1; break;
        case PROFILE_FILE_FLUX:    names = profile_col_names_flux;    pf->has_z = 1; break;
        case PROFILE_FILE_SUBGRID: names = profile_col_names_subgrid; pf->has_z = 1; break;
        default: fclose(fp); return -1;
    }
    pf->ncols = 0;
    while (names[pf->ncols]) {
        strncpy(pf->col_names[pf->ncols], names[pf->ncols], 31);
        pf->col_names[pf->ncols][31] = '\0';
        pf->ncols++;
    }
    strncpy(pf->filename, profile_file_labels[file_type], 255);

    /* Allocate data buffer */
    int capacity = 4096;
    pf->data = (double *)malloc(capacity * pf->ncols * sizeof(double));
    if (!pf->data) { fclose(fp); return -1; }
    pf->nrows = 0;

    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        /* Skip header lines (start with non-numeric after whitespace) */
        char *p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\r') continue;
        /* Try to detect header: if first token can't be parsed as a double, skip */
        char *end;
        strtod(p, &end);
        if (end == p) continue;  /* Not a number, skip header row */

        /* Parse pf->ncols doubles */
        if (pf->nrows >= capacity) {
            if (capacity >= MAX_PROFILE_ROWS) break;
            capacity *= 2;
            if (capacity > MAX_PROFILE_ROWS) capacity = MAX_PROFILE_ROWS;
            double *tmp = (double *)realloc(pf->data, capacity * pf->ncols * sizeof(double));
            if (!tmp) break;
            pf->data = tmp;
        }

        double *row = pf->data + pf->nrows * pf->ncols;
        char *tok = strtok(line, " \t\n\r");
        int col = 0;
        while (tok && col < pf->ncols) {
            row[col++] = strtod(tok, NULL);
            tok = strtok(NULL, " \t\n\r");
        }
        if (col < pf->ncols) {
            /* Partial row: zero-fill */
            while (col < pf->ncols) row[col++] = 0.0;
        }
        pf->nrows++;
    }
    fclose(fp);

    if (pf->nrows == 0) {
        free(pf->data); pf->data = NULL;
        return -1;
    }

    /* Build unique times array and determine nz */
    if (pf->has_z) {
        /* Count unique times (time is column 0) */
        int times_cap = 4096;
        pf->times = (double *)malloc(times_cap * sizeof(double));
        pf->ntimes = 0;
        pf->z_min = pf->data[1];
        pf->z_max = pf->data[1];
        for (int r = 0; r < pf->nrows; r++) {
            double t = pf->data[r * pf->ncols + 0];
            double z = pf->data[r * pf->ncols + 1];
            if (z < pf->z_min) pf->z_min = z;
            if (z > pf->z_max) pf->z_max = z;
            /* Check if t already seen */
            int found = 0;
            for (int i = 0; i < pf->ntimes; i++) {
                if (fabs(pf->times[i] - t) < 1e-15 * (1.0 + fabs(t))) { found = 1; break; }
            }
            if (!found) {
                if (pf->ntimes >= times_cap) {
                    times_cap *= 2;
                    double *tmp = (double *)realloc(pf->times, times_cap * sizeof(double));
                    if (!tmp) break;
                    pf->times = tmp;
                }
                pf->times[pf->ntimes++] = t;
            }
        }
        /* nz = rows per timestep (assume uniform) */
        if (pf->ntimes > 0)
            pf->nz = pf->nrows / pf->ntimes;
        else
            pf->nz = pf->nrows;
    } else {
        /* surf: all rows are individual timesteps */
        pf->ntimes = pf->nrows;
        pf->nz = 1;
        pf->times = (double *)malloc(pf->nrows * sizeof(double));
        for (int r = 0; r < pf->nrows; r++)
            pf->times[r] = pf->data[r * pf->ncols + 0];
    }

    printf("Profile: loaded %s (%d rows, %d cols, %d timesteps, %d z-levels)\n",
           pf->filename, pf->nrows, pf->ncols, pf->ntimes, pf->nz);
    return 0;
}

void free_profile_file(ProfileFile *pf) {
    if (pf->data)  { free(pf->data);  pf->data  = NULL; }
    if (pf->times) { free(pf->times); pf->times = NULL; }
    pf->nrows = 0; pf->ntimes = 0;
}

/* Helper: nice tick labels */
void profile_fmt_val(char *buf, int bufsz, double v) {
    double av = fabs(v);
    if (av == 0.0)                    snprintf(buf, bufsz, "0");
    else if (av >= 1e4 || av < 1e-3) snprintf(buf, bufsz, "%.2e", v);
    else if (av >= 100)               snprintf(buf, bufsz, "%.1f", v);
    else if (av >= 1)                 snprintf(buf, bufsz, "%.3g", v);
    else                              snprintf(buf, bufsz, "%.3g", v);
}

void draw_profile_plot(Display *dpy, Window win, GC plot_gc, ProfileData *pd,
                       int file_idx, int col_idx, int time_idx,
                       int width, int height, int log_x) {
    ProfileFile *pf = &pd->files[file_idx];

    XSetForeground(dpy, plot_gc, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, plot_gc, 0, 0, width, height);
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    if (font) XSetFont(dpy, plot_gc, font->fid);

    if (!pf || pf->nrows == 0 || col_idx < 0 || col_idx >= pf->ncols) {
        XDrawString(dpy, win, plot_gc, width/2-20, height/2, "No data", 7);
        XFlush(dpy); return;
    }

    int plot_left = 75, plot_right = width - 20;
    int plot_top = 30, plot_bottom = height - 45;
    int plot_w = plot_right - plot_left;
    int plot_h = plot_bottom - plot_top;
    if (plot_w <= 0 || plot_h <= 0) return;

    char label[64];

    if (pf->has_z) {
        double t_target = (time_idx < pf->ntimes) ? pf->times[time_idx] : pf->times[0];

        /* Collect (val, z) for this timestep */
        int n = 0;
        double vals_buf[8192], z_buf[8192];
        for (int r = 0; r < pf->nrows && n < 8192; r++) {
            double t = pf->data[r * pf->ncols + 0];
            if (fabs(t - t_target) > 1e-10 * (1.0 + fabs(t_target))) continue;
            z_buf[n]    = pf->data[r * pf->ncols + 1];
            vals_buf[n] = pf->data[r * pf->ncols + col_idx];
            n++;
        }
        if (n == 0) {
            XDrawString(dpy, win, plot_gc, 10, height/2, "No data for this timestep", 25);
            XFlush(dpy); return;
        }

        double vmin = vals_buf[0], vmax = vals_buf[0];
        for (int i = 1; i < n; i++) {
            if (vals_buf[i] < vmin) vmin = vals_buf[i];
            if (vals_buf[i] > vmax) vmax = vals_buf[i];
        }
        if (vmin == vmax) { vmin -= 0.5; vmax += 0.5; }

        double zmin = pf->z_min, zmax = pf->z_max;
        if (zmin == zmax) { zmax = zmin + 1.0; }

        int use_log = log_x && vmin > 0.0;
        double lv0 = use_log ? log10(vmin) : 0;
        double lv1 = use_log ? log10(vmax) : 0;

        /* Axes */
        XDrawLine(dpy, win, plot_gc, plot_left, plot_bottom, plot_right, plot_bottom);
        XDrawLine(dpy, win, plot_gc, plot_left, plot_top,    plot_left,  plot_bottom);
        if (!use_log && vmin < 0.0 && vmax > 0.0) {
            int x0 = plot_left + (int)((-vmin) / (vmax - vmin) * plot_w);
            XSetForeground(dpy, plot_gc, 0xBBBBBB);
            XDrawLine(dpy, win, plot_gc, x0, plot_top, x0, plot_bottom);
            XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
        }

        /* X ticks */
        for (int i = 0; i <= 5; i++) {
            double f = (double)i / 5;
            double val = use_log ? pow(10.0, lv0 + f*(lv1-lv0)) : vmin + f*(vmax-vmin);
            int xp = plot_left + (int)(f * plot_w);
            XDrawLine(dpy, win, plot_gc, xp, plot_bottom, xp, plot_bottom+3);
            profile_fmt_val(label, sizeof(label), val);
            int lw = font ? XTextWidth(font, label, strlen(label)) : 40;
            XDrawString(dpy, win, plot_gc, xp - lw/2, plot_bottom+16, label, strlen(label));
        }
        /* Y (z) ticks */
        for (int i = 0; i <= 6; i++) {
            double z = zmin + (zmax-zmin) * i / 6;
            int yp = plot_bottom - (int)((z-zmin)/(zmax-zmin) * plot_h);
            XDrawLine(dpy, win, plot_gc, plot_left-3, yp, plot_left, yp);
            profile_fmt_val(label, sizeof(label), z);
            int lw = font ? XTextWidth(font, label, strlen(label)) : 40;
            XDrawString(dpy, win, plot_gc, plot_left-lw-5, yp+4, label, strlen(label));
        }
        /* Axis labels */
        int xlw = font ? XTextWidth(font, pf->col_names[col_idx], strlen(pf->col_names[col_idx])) : 40;
        XDrawString(dpy, win, plot_gc, plot_left+(plot_w-xlw)/2, plot_bottom+32, pf->col_names[col_idx], strlen(pf->col_names[col_idx]));
        XDrawString(dpy, win, plot_gc, 2, plot_top+10, "z(m)", 4);
        char tlabel[64];
        snprintf(tlabel, sizeof(tlabel), "t=%.4g s  [%d/%d]", t_target, time_idx+1, pf->ntimes);
        XDrawString(dpy, win, plot_gc, plot_left+4, plot_top-2, tlabel, strlen(tlabel));

        /* Profile line */
        XSetForeground(dpy, plot_gc, 0x0000CC);
        int px0 = -1, py0 = -1;
        for (int i = 0; i < n; i++) {
            double fx, fy;
            if (use_log) {
                if (vals_buf[i] <= 0.0) { px0=-1; continue; }
                fx = (log10(vals_buf[i]) - lv0) / (lv1 - lv0);
            } else {
                fx = (vals_buf[i] - vmin) / (vmax - vmin);
            }
            fy = (z_buf[i] - zmin) / (zmax - zmin);
            int px = plot_left + (int)(fx * plot_w);
            int py = plot_bottom - (int)(fy * plot_h);
            if (px < plot_left)   px = plot_left;
            if (px > plot_right)  px = plot_right;
            if (py < plot_top)    py = plot_top;
            if (py > plot_bottom) py = plot_bottom;
            if (px0 >= 0) XDrawLine(dpy, win, plot_gc, px0, py0, px, py);
            XFillRectangle(dpy, win, plot_gc, px-1, py-1, 3, 3);
            px0 = px; py0 = py;
        }

    } else {
        /* Surf: time-series. X=time, Y=variable */
        if (col_idx == 0) col_idx = 1;
        if (col_idx >= pf->ncols) { XDrawString(dpy,win,plot_gc,10,height/2,"Select variable",15); XFlush(dpy); return; }

        double tmin = pf->times[0], tmax = pf->times[pf->ntimes-1];
        if (tmin == tmax) tmax = tmin+1.0;
        double vmin = pf->data[col_idx], vmax = pf->data[col_idx];
        for (int r = 0; r < pf->nrows; r++) {
            double v = pf->data[r*pf->ncols+col_idx];
            if (v < vmin) vmin=v; if (v > vmax) vmax=v;
        }
        if (vmin == vmax) { vmin -= 0.5; vmax += 0.5; }

        XDrawLine(dpy, win, plot_gc, plot_left, plot_bottom, plot_right, plot_bottom);
        XDrawLine(dpy, win, plot_gc, plot_left, plot_top,    plot_left,  plot_bottom);
        for (int i = 0; i <= 5; i++) {
            double t = tmin + (tmax-tmin)*i/5;
            int xp = plot_left + (int)((double)i/5 * plot_w);
            XDrawLine(dpy, win, plot_gc, xp, plot_bottom, xp, plot_bottom+3);
            profile_fmt_val(label, sizeof(label), t);
            int lw = font ? XTextWidth(font, label, strlen(label)) : 40;
            XDrawString(dpy, win, plot_gc, xp-lw/2, plot_bottom+16, label, strlen(label));
        }
        for (int i = 0; i <= 5; i++) {
            double v = vmin + (vmax-vmin)*i/5;
            int yp = plot_bottom - (int)((double)i/5 * plot_h);
            XDrawLine(dpy, win, plot_gc, plot_left-3, yp, plot_left, yp);
            profile_fmt_val(label, sizeof(label), v);
            int lw = font ? XTextWidth(font, label, strlen(label)) : 40;
            XDrawString(dpy, win, plot_gc, plot_left-lw-5, yp+4, label, strlen(label));
        }
        const char *xl = "time (s)";
        int xlw = font ? XTextWidth(font, xl, strlen(xl)) : 50;
        XDrawString(dpy, win, plot_gc, plot_left+(plot_w-xlw)/2, plot_bottom+32, xl, strlen(xl));
        XDrawString(dpy, win, plot_gc, 2, plot_top+10, pf->col_names[col_idx], strlen(pf->col_names[col_idx]));

        XSetForeground(dpy, plot_gc, 0x0000CC);
        int px0=-1, py0=-1;
        for (int r = 0; r < pf->nrows; r++) {
            double t = pf->data[r*pf->ncols+0];
            double v = pf->data[r*pf->ncols+col_idx];
            int px = plot_left + (int)((t-tmin)/(tmax-tmin)*plot_w);
            int py = plot_bottom - (int)((v-vmin)/(vmax-vmin)*plot_h);
            if (px<plot_left) px=plot_left; if (px>plot_right) px=plot_right;
            if (py<plot_top)  py=plot_top;  if (py>plot_bottom) py=plot_bottom;
            if (px0>=0) XDrawLine(dpy, win, plot_gc, px0, py0, px, py);
            XFillRectangle(dpy, win, plot_gc, px-1, py-1, 3, 3);
            px0=px; py0=py;
        }
    }
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XFlush(dpy);
}

void draw_profile_contour(Display *dpy, Window win, GC plot_gc, ProfileData *pd,
                          int file_idx, int col_idx, int width, int height) {
    ProfileFile *pf = &pd->files[file_idx];

    XSetForeground(dpy, plot_gc, WhitePixel(dpy, screen));
    XFillRectangle(dpy, win, plot_gc, 0, 0, width, height);
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    if (font) XSetFont(dpy, plot_gc, font->fid);

    if (!pf || !pf->has_z || pf->ntimes < 1 || pf->nz < 1 ||
        col_idx < 2 || col_idx >= pf->ncols) {
        const char *msg = "No profile data for contour";
        XDrawString(dpy, win, plot_gc, width/2-60, height/2, msg, strlen(msg));
        XFlush(dpy); return;
    }

    int plot_left = 75, plot_right = width - 95;  /* leave ~90px for colorbar + labels */
    int plot_top = 30, plot_bottom = height - 45;
    int plot_w = plot_right - plot_left;
    int plot_h = plot_bottom - plot_top;
    if (plot_w <= 0 || plot_h <= 0) return;

    /* Find global min/max of the variable across all timesteps */
    double vmin = 1e300, vmax = -1e300;
    for (int r = 0; r < pf->nrows; r++) {
        double v = pf->data[r * pf->ncols + col_idx];
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    if (vmin == vmax) { vmin -= 0.5; vmax += 0.5; }

    double zmin = pf->z_min, zmax = pf->z_max;
    if (zmin == zmax) zmax = zmin + 1.0;
    double tmin = pf->times[0], tmax = pf->times[pf->ntimes-1];
    if (tmin == tmax) tmax = tmin + 1.0;

    /* Render pixel by pixel using XImage */
    XImage *img = XCreateImage(dpy, DefaultVisual(dpy, screen), DefaultDepth(dpy, screen),
                               ZPixmap, 0, NULL, plot_w, plot_h, 32, 0);
    if (!img) goto draw_axes_only;
    img->data = (char *)malloc(img->bytes_per_line * plot_h);
    if (!img->data) { XDestroyImage(img); goto draw_axes_only; }

    /* Build lookup: for each pixel column, which time index? */
    /* For each pixel row, which z index (nearest neighbor)? */
    for (int py = 0; py < plot_h; py++) {
        /* z from top=zmax to bottom=zmin */
        double z_frac = 1.0 - (double)py / (plot_h - 1);
        double z_target = zmin + z_frac * (zmax - zmin);

        for (int px = 0; px < plot_w; px++) {
            /* map px to time index */
            double t_frac = (double)px / (plot_w - 1);
            double t_target = tmin + t_frac * (tmax - tmin);

            /* Find nearest time index */
            int ti = 0;
            double best_dt = fabs(pf->times[0] - t_target);
            for (int i = 1; i < pf->ntimes; i++) {
                double dt = fabs(pf->times[i] - t_target);
                if (dt < best_dt) { best_dt = dt; ti = i; }
            }

            /* Find nearest z row within this timestep */
            /* Rows for timestep ti start at ti*nz (assuming uniform layout) */
            int row_start = ti * pf->nz;
            int zi = 0;
            double best_dz = 1e300;
            for (int i = 0; i < pf->nz && (row_start + i) < pf->nrows; i++) {
                double z = pf->data[(row_start + i) * pf->ncols + 1];
                double dz = fabs(z - z_target);
                if (dz < best_dz) { best_dz = dz; zi = i; }
            }

            int row = row_start + zi;
            double v = (row < pf->nrows) ? pf->data[row * pf->ncols + col_idx] : 0.0;
            double t = (v - vmin) / (vmax - vmin);
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;

            RGB rgb = viridis_colormap(t);
            unsigned long pixel = ((unsigned long)(rgb.r) << 16) |
                                  ((unsigned long)(rgb.g) << 8)  |
                                  ((unsigned long)(rgb.b));
            XPutPixel(img, px, py, pixel);
        }
    }
    XPutImage(dpy, win, plot_gc, img, 0, 0, plot_left, plot_top, plot_w, plot_h);
    free(img->data); img->data = NULL;
    XDestroyImage(img);

draw_axes_only:
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, plot_gc, plot_left, plot_top, plot_w, plot_h);

    /* X (time) ticks */
    char label[64];
    for (int i = 0; i <= 5; i++) {
        double t = tmin + (tmax-tmin)*i/5;
        int xp = plot_left + (int)((double)i/5 * plot_w);
        XDrawLine(dpy, win, plot_gc, xp, plot_bottom, xp, plot_bottom+3);
        profile_fmt_val(label, sizeof(label), t);
        int lw = font ? XTextWidth(font, label, strlen(label)) : 40;
        XDrawString(dpy, win, plot_gc, xp-lw/2, plot_bottom+16, label, strlen(label));
    }
    /* Y (z) ticks */
    for (int i = 0; i <= 6; i++) {
        double z = zmin + (zmax-zmin)*i/6;
        int yp = plot_bottom - (int)((z-zmin)/(zmax-zmin)*plot_h);
        XDrawLine(dpy, win, plot_gc, plot_left-3, yp, plot_left, yp);
        profile_fmt_val(label, sizeof(label), z);
        int lw = font ? XTextWidth(font, label, strlen(label)) : 40;
        XDrawString(dpy, win, plot_gc, plot_left-lw-5, yp+4, label, strlen(label));
    }
    const char *xl = "time (s)";
    int xlw = font ? XTextWidth(font, xl, strlen(xl)) : 50;
    XDrawString(dpy, win, plot_gc, plot_left+(plot_w-xlw)/2, plot_bottom+32, xl, strlen(xl));
    XDrawString(dpy, win, plot_gc, 2, plot_top+10, "z(m)", 4);

    /* Colorbar: right side strip (cb_x is 8px right of plot area) */
    int cb_x = plot_right + 8, cb_w = 14;
    for (int py = 0; py < plot_h; py++) {
        double t = 1.0 - (double)py / (plot_h - 1);
        RGB rgb = viridis_colormap(t);
        unsigned long px_col = ((unsigned long)(rgb.r)<<16)|((unsigned long)(rgb.g)<<8)|(rgb.b);
        XSetForeground(dpy, plot_gc, px_col);
        XFillRectangle(dpy, win, plot_gc, cb_x, plot_top + py, cb_w, 1);
    }
    XSetForeground(dpy, plot_gc, BlackPixel(dpy, screen));
    XDrawRectangle(dpy, win, plot_gc, cb_x, plot_top, cb_w, plot_h);
    /* Colorbar tick labels: max, mid, min */
    int lx = cb_x + cb_w + 4;
    profile_fmt_val(label, sizeof(label), vmax);
    XDrawLine(dpy, win, plot_gc, cb_x, plot_top, cb_x + cb_w + 3, plot_top);
    XDrawString(dpy, win, plot_gc, lx, plot_top + 10, label, strlen(label));
    profile_fmt_val(label, sizeof(label), (vmin + vmax) * 0.5);
    int mid_y = plot_top + plot_h / 2;
    XDrawLine(dpy, win, plot_gc, cb_x, mid_y, cb_x + cb_w + 3, mid_y);
    XDrawString(dpy, win, plot_gc, lx, mid_y + 4, label, strlen(label));
    profile_fmt_val(label, sizeof(label), vmin);
    XDrawLine(dpy, win, plot_gc, cb_x, plot_bottom, cb_x + cb_w + 3, plot_bottom);
    XDrawString(dpy, win, plot_gc, lx, plot_bottom + 4, label, strlen(label));

    /* Title */
    snprintf(label, sizeof(label), "%s  [contour]", pf->col_names[col_idx]);
    XDrawString(dpy, win, plot_gc, plot_left+4, plot_top-2, label, strlen(label));

    XFlush(dpy);
}

void render_profile_canvas(ProfileData *pd) {
    if (!pd || profile_canvas == 0) return;
    int file_idx = profile_current_file;
    if (!pd->loaded[file_idx]) {
        for (int i = 0; i < N_PROFILE_FILES; i++)
            if (pd->loaded[i]) { file_idx = i; break; }
    }
    unsigned int w = 0, h = 0;
    Window root; int x, y; unsigned int bw, depth;
    XGetGeometry(display, profile_canvas, &root, &x, &y, &w, &h, &bw, &depth);
    if ((int)w <= 0 || (int)h <= 0) return;
    GC gc2 = XCreateGC(display, profile_canvas, 0, NULL);
    if (font) XSetFont(display, gc2, font->fid);
    if (profile_contour_mode && pd->files[file_idx].has_z)
        draw_profile_contour(display, profile_canvas, gc2, pd, file_idx, profile_current_col, w, h);
    else
        draw_profile_plot(display, profile_canvas, gc2, pd, file_idx, profile_current_col,
                          profile_current_time_idx, w, h, profile_log_x);
    XFreeGC(display, gc2);
}

void update_profile_info_label(ProfileData *pd) {
    if (!pd || !profile_info_label) return;
    int fi = profile_current_file;
    if (!pd->loaded[fi]) { XtVaSetValues(profile_info_label, XtNlabel, "No data for this file", NULL); return; }
    ProfileFile *pf = &pd->files[fi];
    char buf[256];
    const char *mode = profile_contour_mode ? "contour" : "profile";
    if (pf->has_z && pf->ntimes > 0 && !profile_contour_mode) {
        double t = pf->times[profile_current_time_idx];
        snprintf(buf, sizeof(buf), "%s | %s | t=%.4g s [%d/%d] | %s",
                 pf->filename, pf->col_names[profile_current_col],
                 t, profile_current_time_idx+1, pf->ntimes, mode);
    } else {
        snprintf(buf, sizeof(buf), "%s | %s | %d timesteps | %s",
                 pf->filename, pf->col_names[profile_current_col], pf->ntimes, mode);
    }
    XtVaSetValues(profile_info_label, XtNlabel, buf, NULL);
}

void profile_canvas_expose_callback(Widget w, XtPointer client_data, XEvent *event, Boolean *cont) {
    (void)w; (void)client_data; (void)event; (void)cont;
    render_profile_canvas(global_profile);
}

void profile_rebuild_var_buttons(ProfileData *pd) {
    if (!profile_var_viewport) return;

    /* Destroy old inner box and replace with a fresh one */
    if (profile_var_box_widget) {
        XtDestroyWidget(profile_var_box_widget);
        profile_var_box_widget = NULL;
    }
    if (profile_var_buttons) { free(profile_var_buttons); profile_var_buttons = NULL; }
    profile_n_var_buttons = 0;

    int fi = profile_current_file;
    if (!pd || !pd->loaded[fi]) return;
    ProfileFile *pf = &pd->files[fi];

    /* New inner box */
    Arg args[4]; int n;
    n = 0;
    XtSetArg(args[n], XtNorientation, XtorientVertical); n++;
    profile_var_box_widget = XtCreateManagedWidget("varInnerBox", boxWidgetClass,
                                                    profile_var_viewport, args, n);

    int start = pf->has_z ? 2 : 1;
    int nbtns = pf->ncols - start;
    if (nbtns <= 0) { XtRealizeWidget(profile_var_box_widget); return; }

    profile_var_buttons = (Widget *)malloc(nbtns * sizeof(Widget));
    profile_n_var_buttons = nbtns;

    for (int i = 0; i < nbtns; i++) {
        int col = start + i;
        n = 0;
        XtSetArg(args[n], XtNlabel, pf->col_names[col]); n++;
        profile_var_buttons[i] = XtCreateManagedWidget(pf->col_names[col],
            commandWidgetClass, profile_var_box_widget, args, n);
        XtAddCallback(profile_var_buttons[i], XtNcallback, profile_var_callback, (XtPointer)(long)col);
    }
    XtRealizeWidget(profile_var_box_widget);
}

void profile_file_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)call_data;
    int fi = (int)(long)client_data;
    if (!global_profile || !global_profile->loaded[fi]) return;
    profile_current_file = fi;
    ProfileFile *pf = &global_profile->files[fi];
    profile_current_col = pf->has_z ? 2 : 1;
    if (profile_current_col >= pf->ncols) profile_current_col = 0;
    profile_current_time_idx = 0;
    profile_rebuild_var_buttons(global_profile);
    update_profile_info_label(global_profile);
    render_profile_canvas(global_profile);
}

void profile_var_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)call_data;
    profile_current_col = (int)(long)client_data;
    update_profile_info_label(global_profile);
    render_profile_canvas(global_profile);
}

void profile_time_nav_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)call_data;
    int dir = (int)(long)client_data;
    if (!global_profile) return;
    int fi = profile_current_file;
    if (!global_profile->loaded[fi]) return;
    ProfileFile *pf = &global_profile->files[fi];
    if (!pf->has_z) return;
    profile_current_time_idx += dir;
    if (profile_current_time_idx < 0) profile_current_time_idx = pf->ntimes - 1;
    if (profile_current_time_idx >= pf->ntimes) profile_current_time_idx = 0;
    update_profile_info_label(global_profile);
    render_profile_canvas(global_profile);
}

void profile_logx_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)client_data; (void)call_data;
    profile_log_x = !profile_log_x;
    render_profile_canvas(global_profile);
}

void profile_contour_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    (void)w; (void)client_data; (void)call_data;
    profile_contour_mode = !profile_contour_mode;
    update_profile_info_label(global_profile);
    render_profile_canvas(global_profile);
}

void init_profile_gui(ProfileData *pd, int argc, char **argv) {
    Arg args[20];
    int n;

    global_profile = pd;

    toplevel = XtAppInitialize(NULL, "PLTView-Profile", NULL, 0, &argc, argv, NULL, NULL, 0);
    display = XtDisplay(toplevel);
    screen  = DefaultScreen(display);

    font = XLoadQueryFont(display, "fixed");
    if (!font) font = XLoadQueryFont(display, "*");

    /* Compute var viewport width from longest column name */
    int vp_width = 120;
    if (font) {
        int fi = profile_current_file;
        if (pd->loaded[fi]) {
            ProfileFile *pf = &pd->files[fi];
            for (int i = 0; i < pf->ncols; i++) {
                int tw = XTextWidth(font, pf->col_names[i], strlen(pf->col_names[i])) + 24;
                if (tw > vp_width) vp_width = tw;
            }
        }
        /* also check all files for widest name */
        for (int fi2 = 0; fi2 < N_PROFILE_FILES; fi2++) {
            if (!pd->loaded[fi2]) continue;
            ProfileFile *pf = &pd->files[fi2];
            for (int i = 0; i < pf->ncols; i++) {
                int tw = XTextWidth(font, pf->col_names[i], strlen(pf->col_names[i])) + 24;
                if (tw > vp_width) vp_width = tw;
            }
        }
        if (vp_width > 160) vp_width = 160;
    }
    int canvas_w = 820, canvas_h = 520;

    /* Main form */
    n = 0;
    XtSetArg(args[n], XtNwidth,  canvas_w + vp_width + 20); n++;
    XtSetArg(args[n], XtNheight, canvas_h + 100); n++;
    form = XtCreateManagedWidget("form", formWidgetClass, toplevel, args, n);

    /* Info label at top spanning full width */
    n = 0;
    XtSetArg(args[n], XtNlabel, "1D Profiles - Loading..."); n++;
    XtSetArg(args[n], XtNwidth, canvas_w + vp_width); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNtop,   XawChainTop);  n++;
    XtSetArg(args[n], XtNleft,  XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    profile_info_label = XtCreateManagedWidget("info", labelWidgetClass, form, args, n);

    /* Variable viewport: left column, below info label */
    n = 0;
    XtSetArg(args[n], XtNfromVert,    profile_info_label); n++;
    XtSetArg(args[n], XtNwidth,       vp_width); n++;
    XtSetArg(args[n], XtNheight,      canvas_h); n++;
    XtSetArg(args[n], XtNallowVert,   True); n++;
    XtSetArg(args[n], XtNallowHoriz,  False); n++;
    XtSetArg(args[n], XtNforceBars,   True); n++;
    XtSetArg(args[n], XtNtop,         XawChainTop); n++;
    XtSetArg(args[n], XtNbottom,      XawChainBottom); n++;
    XtSetArg(args[n], XtNleft,        XawChainLeft); n++;
    profile_var_viewport = XtCreateManagedWidget("varViewport", viewportWidgetClass, form, args, n);

    /* Canvas: to the right of var viewport */
    n = 0;
    XtSetArg(args[n], XtNfromVert,    profile_info_label); n++;
    XtSetArg(args[n], XtNfromHoriz,   profile_var_viewport); n++;
    XtSetArg(args[n], XtNwidth,       canvas_w); n++;
    XtSetArg(args[n], XtNheight,      canvas_h); n++;
    XtSetArg(args[n], XtNborderWidth, 2); n++;
    XtSetArg(args[n], XtNtop,         XawChainTop); n++;
    XtSetArg(args[n], XtNbottom,      XawChainBottom); n++;
    XtSetArg(args[n], XtNleft,        XawChainLeft); n++;
    XtSetArg(args[n], XtNright,       XawChainRight); n++;
    profile_canvas_widget = XtCreateManagedWidget("profileCanvas", simpleWidgetClass, form, args, n);

    /* Controls row below both: file buttons + nav + options */
    Widget ctrl_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert,    profile_canvas_widget); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNbottom,      XawChainBottom); n++;
    XtSetArg(args[n], XtNleft,        XawChainLeft); n++;
    XtSetArg(args[n], XtNright,       XawChainRight); n++;
    ctrl_box = XtCreateManagedWidget("ctrlBox", boxWidgetClass, form, args, n);

    /* File label + buttons */
    n=0; XtSetArg(args[n], XtNlabel, "File:"); n++; XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("fileLabel", labelWidgetClass, ctrl_box, args, n);
    for (int i = 0; i < N_PROFILE_FILES; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, profile_file_labels[i]); n++;
        if (!pd->loaded[i]) { XtSetArg(args[n], XtNsensitive, False); n++; }
        profile_file_buttons[i] = XtCreateManagedWidget(profile_file_labels[i],
            commandWidgetClass, ctrl_box, args, n);
        XtAddCallback(profile_file_buttons[i], XtNcallback, profile_file_callback, (XtPointer)(long)i);
    }

    /* Separator label */
    n=0; XtSetArg(args[n], XtNlabel, " | "); n++; XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("sep1", labelWidgetClass, ctrl_box, args, n);

    /* Timestep nav */
    n=0; XtSetArg(args[n], XtNlabel, "t:"); n++; XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("tLabel", labelWidgetClass, ctrl_box, args, n);
    n=0; XtSetArg(args[n], XtNlabel, "<"); n++;
    Widget pb = XtCreateManagedWidget("prev", commandWidgetClass, ctrl_box, args, n);
    XtAddCallback(pb, XtNcallback, profile_time_nav_callback, (XtPointer)-1L);
    n=0; XtSetArg(args[n], XtNlabel, ">"); n++;
    Widget nb = XtCreateManagedWidget("next", commandWidgetClass, ctrl_box, args, n);
    XtAddCallback(nb, XtNcallback, profile_time_nav_callback, (XtPointer)1L);

    /* Separator */
    n=0; XtSetArg(args[n], XtNlabel, " | "); n++; XtSetArg(args[n], XtNborderWidth, 0); n++;
    XtCreateManagedWidget("sep2", labelWidgetClass, ctrl_box, args, n);

    /* LogX */
    n=0; XtSetArg(args[n], XtNlabel, "LogX"); n++;
    Widget lx = XtCreateManagedWidget("logX", commandWidgetClass, ctrl_box, args, n);
    XtAddCallback(lx, XtNcallback, profile_logx_callback, NULL);

    /* Contour toggle */
    n=0; XtSetArg(args[n], XtNlabel, "Contour"); n++;
    Widget ct = XtCreateManagedWidget("contour", commandWidgetClass, ctrl_box, args, n);
    XtAddCallback(ct, XtNcallback, profile_contour_callback, NULL);

    /* Realize */
    XtRealizeWidget(toplevel);

    profile_canvas = XtWindow(profile_canvas_widget);
    XtAddEventHandler(profile_canvas_widget, ExposureMask, False,
                      profile_canvas_expose_callback, NULL);
    XSelectInput(display, profile_canvas, ExposureMask | KeyPressMask);

    /* Build initial var buttons inside viewport */
    profile_rebuild_var_buttons(pd);
}

int main(int argc, char **argv) {
    PlotfileData pf = {0};
    Arg args[2];
    char check_path[MAX_PATH];
    const char *prefix = "plt";  /* Default prefix */

    init_map_layers_dir(argv[0]);

    /* Check for --sdm, --sbm, --profile flags */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--sdm") == 0) {
            sdm_mode = 1;
            for (int j = i; j < argc - 1; j++) argv[j] = argv[j + 1];
            argc--; i--;
        } else if (strcmp(argv[i], "--sbm") == 0) {
            sbm_mode = 1;
            for (int j = i; j < argc - 1; j++) argv[j] = argv[j + 1];
            argc--; i--;
        } else if (strcmp(argv[i], "--profile") == 0) {
            profile_mode = 1;
            for (int j = i; j < argc - 1; j++) argv[j] = argv[j + 1];
            argc--; i--;
        }
    }

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [--sdm|--sbm|--profile] <directory> [prefix]\n", argv[0]);
        fprintf(stderr, "  Single plotfile:    %s plt00100\n", argv[0]);
        fprintf(stderr, "  Multi-timestep:     %s /path/to/dir plt\n", argv[0]);
        fprintf(stderr, "  With prefix plt2d:  %s /path/to/dir plt2d\n", argv[0]);
        fprintf(stderr, "  SDM mode:           %s --sdm plt00100\n", argv[0]);
        fprintf(stderr, "  SDM multi-timestep: %s --sdm /path/to/dir plt\n", argv[0]);
        fprintf(stderr, "  SBM mode:           %s --sbm plt00100\n", argv[0]);
        fprintf(stderr, "  SBM multi-timestep: %s --sbm /path/to/dir plt\n", argv[0]);
        fprintf(stderr, "  1D Profiles:        %s --profile /path/to/dir\n", argv[0]);
        return 1;
    }

    /* Get prefix from argument if provided */
    if (argc >= 3) {
        prefix = argv[2];
    }

    /* Check if argument is a single plotfile or a directory containing plotfiles */
    if (sdm_mode) {
        /* SDM mode: check for super_droplets_moisture/Header */
        snprintf(check_path, MAX_PATH, "%s/%s/Header", argv[1], SDM_SUBDIR);
        FILE *fp = fopen(check_path, "r");

        if (fp) {
            /* Single plotfile with SDM data */
            fclose(fp);
            n_timesteps = 1;
            current_timestep = 0;
            timestep_paths[0] = strdup(argv[1]);
            timestep_numbers[0] = 0;
            printf("SDM single plotfile mode: %s\n", argv[1]);
        } else {
            /* Multi-timestep: try explicit prefix first, then auto-detect */
            int found = 0;
            if (argc >= 3) {
                /* Explicit prefix given */
                printf("Scanning for SDM plotfiles with prefix '%s'...\n", prefix);
                found = scan_sdm_timesteps(argv[1], prefix) > 0;
            }
            if (!found) {
                /* Auto-detect prefix from first directory with SDM data */
                DIR *autodir = opendir(argv[1]);
                if (autodir) {
                    struct dirent *ae;
                    char detected_prefix[128] = "";
                    while ((ae = readdir(autodir)) != NULL) {
                        snprintf(check_path, MAX_PATH, "%s/%s/%s/Header",
                                 argv[1], ae->d_name, SDM_SUBDIR);
                        FILE *af = fopen(check_path, "r");
                        if (af) {
                            fclose(af);
                            /* Extract prefix: everything before trailing digits */
                            const char *name = ae->d_name;
                            int len = strlen(name);
                            int end = len;
                            while (end > 0 && isdigit(name[end - 1])) end--;
                            if (end > 0 && end < len) {
                                strncpy(detected_prefix, name, end);
                                detected_prefix[end] = '\0';
                                break;
                            }
                        }
                    }
                    closedir(autodir);

                    if (detected_prefix[0]) {
                        printf("Auto-detected SDM prefix: '%s'\n", detected_prefix);
                        prefix = strdup(detected_prefix);
                        found = scan_sdm_timesteps(argv[1], prefix) > 0;
                    }
                }
            }
            if (!found) {
                fprintf(stderr, "Error: No plotfiles with SDM data found in %s\n", argv[1]);
                return 1;
            }
            current_timestep = 0;
            printf("SDM multi-timestep mode: %d timesteps found\n", n_timesteps);
        }
    } else if (sbm_mode) {
        /* SBM mode: check for bin_info.txt */
        snprintf(check_path, MAX_PATH, "%s/%s", argv[1], SBM_BIN_INFO_FILE);
        FILE *fp = fopen(check_path, "r");

        if (fp) {
            /* Single plotfile with SBM data */
            fclose(fp);
            n_timesteps = 1;
            current_timestep = 0;
            timestep_paths[0] = strdup(argv[1]);
            timestep_numbers[0] = 0;
            printf("SBM single plotfile mode: %s\n", argv[1]);
        } else {
            /* Multi-timestep: try explicit prefix first, then auto-detect */
            int found = 0;
            if (argc >= 3) {
                /* Explicit prefix given */
                printf("Scanning for SBM plotfiles with prefix '%s'...\n", prefix);
                found = scan_sbm_timesteps(argv[1], prefix) > 0;
            }
            if (!found) {
                /* Auto-detect prefix from first directory with SBM data */
                DIR *autodir = opendir(argv[1]);
                if (autodir) {
                    struct dirent *ae;
                    char detected_prefix[128] = "";
                    while ((ae = readdir(autodir)) != NULL) {
                        snprintf(check_path, MAX_PATH, "%s/%s/%s",
                                 argv[1], ae->d_name, SBM_BIN_INFO_FILE);
                        FILE *af = fopen(check_path, "r");
                        if (af) {
                            fclose(af);
                            /* Extract prefix: everything before trailing digits */
                            const char *name = ae->d_name;
                            int len = strlen(name);
                            int end = len;
                            while (end > 0 && isdigit(name[end - 1])) end--;
                            if (end > 0 && end < len) {
                                strncpy(detected_prefix, name, end);
                                detected_prefix[end] = '\0';
                                break;
                            }
                        }
                    }
                    closedir(autodir);

                    if (detected_prefix[0]) {
                        printf("Auto-detected SBM prefix: '%s'\n", detected_prefix);
                        prefix = strdup(detected_prefix);
                        found = scan_sbm_timesteps(argv[1], prefix) > 0;
                    }
                }
            }
            if (!found) {
                fprintf(stderr, "Error: No plotfiles with SBM data found in %s\n", argv[1]);
                return 1;
            }
            current_timestep = 0;
            printf("SBM multi-timestep mode: %d timesteps found\n", n_timesteps);
        }
    } else if (!profile_mode) {
        snprintf(check_path, MAX_PATH, "%s/Header", argv[1]);
        FILE *fp = fopen(check_path, "r");

        if (fp) {
            /* Single plotfile mode */
            fclose(fp);
            n_timesteps = 1;
            current_timestep = 0;
            timestep_paths[0] = strdup(argv[1]);
            timestep_numbers[0] = 0;
            strncpy(pf.plotfile_dir, argv[1], MAX_PATH - 1);
            printf("Single plotfile mode: %s\n", argv[1]);
        } else {
            /* Try multi-timestep mode - scan directory for plotfiles */
            printf("Scanning for plotfiles with prefix '%s'...\n", prefix);
            if (scan_timesteps(argv[1], prefix) <= 0) {
                fprintf(stderr, "Error: No valid plotfiles with prefix '%s' found in %s\n", prefix, argv[1]);
                return 1;
            }
            current_timestep = 0;
            strncpy(pf.plotfile_dir, timestep_paths[0], MAX_PATH - 1);
            printf("Multi-timestep mode: %d timesteps found\n", n_timesteps);
        }
    }

    /* SDM mode: particle histogram viewer */
    if (sdm_mode) {
        ParticleData pd = {0};
        pd.current_metric = SDM_METRIC_PARTICLE_COUNT;

        if (read_sdm_header(&pd, timestep_paths[current_timestep]) < 0) {
            fprintf(stderr, "Error: Failed to read SDM header\n");
            return 1;
        }

        pd.domain_volume = compute_domain_volume(timestep_paths[current_timestep]);

        if (read_sdm_data(&pd, timestep_paths[current_timestep]) < 0) {
            fprintf(stderr, "Error: Failed to read SDM data\n");
            return 1;
        }

        init_sdm_gui(&pd, timestep_paths[current_timestep], argc, argv);
        update_sdm_info_label(&pd, timestep_paths[current_timestep]);
        render_sdm_histogram(&pd);

        printf("\nSDM Mode Controls:\n");
        printf("  Click metric buttons to change y-axis\n");
        printf("  Click LogX/LogY to toggle log scale\n");
        printf("  Click Settings to set cutoff radius and bin width\n");
        if (n_timesteps > 1) {
            printf("  Click </> or use Left/Right arrow keys to navigate timesteps\n");
        }
        printf("\n");

        /* SDM event loop */
        XtAppContext app_context = XtWidgetToApplicationContext(toplevel);
        while (1) {
            XEvent event;
            XtAppNextEvent(app_context, &event);

            if (event.type == Expose) {
                if (event.xexpose.window == sdm_canvas && global_pd) {
                    render_sdm_histogram(global_pd);
                    if (!initial_focus_set) {
                        XSetInputFocus(display, sdm_canvas, RevertToParent, CurrentTime);
                        initial_focus_set = 1;
                    }
                }
            } else if (event.type == KeyPress && global_pd) {
                /* Handle keyboard input for SDM settings dialog */
                if (sdm_dialog_active && sdm_active_text_widget) {
                    char buf[32];
                    KeySym keysym;
                    int len = XLookupString(&event.xkey, buf, sizeof(buf) - 1, &keysym, NULL);

                    String current_value;
                    Arg kargs[1];
                    XtSetArg(kargs[0], XtNstring, &current_value);
                    XtGetValues(sdm_active_text_widget, kargs, 1);

                    char new_value[256];
                    strncpy(new_value, current_value ? current_value : "", sizeof(new_value) - 1);
                    new_value[sizeof(new_value) - 1] = '\0';
                    size_t current_len = strlen(new_value);

                    if (keysym == XK_BackSpace || keysym == XK_Delete) {
                        if (current_len > 0) {
                            new_value[current_len - 1] = '\0';
                            XtSetArg(kargs[0], XtNstring, new_value);
                            XtSetValues(sdm_active_text_widget, kargs, 1);
                        }
                    } else if (keysym == XK_Tab) {
                        /* Switch between cutoff and binwidth fields */
                        if (sdm_active_field == 0 && sdm_settings_text_binwidth) {
                            sdm_active_text_widget = sdm_settings_text_binwidth;
                            sdm_active_field = 1;
                        } else if (sdm_settings_text_cutoff) {
                            sdm_active_text_widget = sdm_settings_text_cutoff;
                            sdm_active_field = 0;
                        }
                    } else if (keysym == XK_Return || keysym == XK_KP_Enter) {
                        sdm_settings_apply_callback(NULL, NULL, NULL);
                    } else if (keysym == XK_Escape) {
                        sdm_settings_close_callback(NULL, NULL, NULL);
                    } else if (len > 0 && isprint((unsigned char)buf[0])) {
                        if (current_len + (size_t)len < sizeof(new_value) - 1) {
                            buf[len] = '\0';
                            strcat(new_value, buf);
                            XtSetArg(kargs[0], XtNstring, new_value);
                            XtSetValues(sdm_active_text_widget, kargs, 1);
                        }
                    }
                    continue;
                }
                if (event.xkey.window == sdm_canvas) {
                    KeySym key = XLookupKeysym(&event.xkey, 0);
                    if (key == XK_Right && n_timesteps > 1) {
                        int new_ts = current_timestep + 1;
                        if (new_ts >= n_timesteps) new_ts = 0;
                        sdm_switch_timestep(global_pd, new_ts);
                        update_time_label();
                        continue;
                    } else if (key == XK_Left && n_timesteps > 1) {
                        int new_ts = current_timestep - 1;
                        if (new_ts < 0) new_ts = n_timesteps - 1;
                        sdm_switch_timestep(global_pd, new_ts);
                        update_time_label();
                        continue;
                    }
                }
            }

            XtDispatchEvent(&event);
        }

        /* Cleanup */
        if (pd.radius) free(pd.radius);
        if (pd.multiplicity) free(pd.multiplicity);
        if (pd.mass) free(pd.mass);
        if (sdm_hist_data) {
            if (sdm_hist_data->bin_counts) free(sdm_hist_data->bin_counts);
            if (sdm_hist_data->bin_centers) free(sdm_hist_data->bin_centers);
            if (sdm_hist_data->bin_edges) free(sdm_hist_data->bin_edges);
            free(sdm_hist_data);
        }
        return 0;
    }

    /* SBM mode: spectral bin microphysics histogram viewer */
    if (sbm_mode) {
        SBMData sbm = {0};
        sbm.current_metric = SBM_METRIC_QC_MASS;

        if (read_sbm_bin_info(&sbm, timestep_paths[current_timestep]) < 0) {
            fprintf(stderr, "Error: No bin_info.txt found\n");
            return 1;
        }

        if (compute_sbm_values(&sbm, timestep_paths[current_timestep]) < 0) {
            fprintf(stderr, "Error: Failed to read SBM data\n");
            return 1;
        }

        init_sbm_gui(&sbm, timestep_paths[current_timestep], argc, argv);
        update_sbm_info_label(&sbm, timestep_paths[current_timestep]);
        render_sbm_histogram(&sbm);

        printf("\nSBM Mode Controls:\n");
        printf("  Click metric buttons to change y-axis\n");
        printf("  Click LogX/LogY to toggle log scale\n");
        if (n_timesteps > 1) {
            printf("  Click </> or use Left/Right arrow keys to navigate timesteps\n");
        }
        printf("\n");

        /* SBM event loop */
        XtAppContext app_context = XtWidgetToApplicationContext(toplevel);
        while (1) {
            XEvent event;
            XtAppNextEvent(app_context, &event);

            if (event.type == Expose) {
                if (event.xexpose.window == sbm_canvas && global_sbm) {
                    render_sbm_histogram(global_sbm);
                    if (!initial_focus_set) {
                        XSetInputFocus(display, sbm_canvas, RevertToParent, CurrentTime);
                        initial_focus_set = 1;
                    }
                }
            } else if (event.type == KeyPress && global_sbm) {
                if (event.xkey.window == sbm_canvas) {
                    KeySym key = XLookupKeysym(&event.xkey, 0);
                    if (key == XK_Right && n_timesteps > 1) {
                        int new_ts = current_timestep + 1;
                        if (new_ts >= n_timesteps) new_ts = 0;
                        sbm_switch_timestep(global_sbm, new_ts);
                        update_time_label();
                        continue;
                    } else if (key == XK_Left && n_timesteps > 1) {
                        int new_ts = current_timestep - 1;
                        if (new_ts < 0) new_ts = n_timesteps - 1;
                        sbm_switch_timestep(global_sbm, new_ts);
                        update_time_label();
                        continue;
                    }
                }
            }

            XtDispatchEvent(&event);
        }

        /* Cleanup */
        if (sbm_hist_data) {
            if (sbm_hist_data->bin_counts) free(sbm_hist_data->bin_counts);
            if (sbm_hist_data->bin_centers) free(sbm_hist_data->bin_centers);
            if (sbm_hist_data->bin_edges) free(sbm_hist_data->bin_edges);
            free(sbm_hist_data);
        }
        return 0;
    }

    /* 1D Profile mode */
    if (profile_mode) {
        ProfileData *pd = (ProfileData *)calloc(1, sizeof(ProfileData));
        if (!pd) { fprintf(stderr, "Error: out of memory\n"); return 1; }
        strncpy(pd->dir, argv[1], MAX_PATH - 1);

        /* Try to load each profile file: exact match first, then substring match */
        const char *keywords[N_PROFILE_FILES] = {"surf", "mean", "flux", "subgrid"};
        char matched_paths[N_PROFILE_FILES][MAX_PATH];
        int any_loaded = 0;
        memset(matched_paths, 0, sizeof(matched_paths));

        for (int fi = 0; fi < N_PROFILE_FILES; fi++) {
            char exact[MAX_PATH];
            snprintf(exact, sizeof(exact), "%s/%s", pd->dir, keywords[fi]);
            if (read_profile_file(&pd->files[fi], exact, fi) == 0) {
                pd->loaded[fi] = 1;
                any_loaded = 1;
                strncpy(matched_paths[fi], exact, MAX_PATH - 1);
                continue;
            }
            /* Substring scan of directory entries */
            DIR *dp = opendir(pd->dir);
            if (dp) {
                struct dirent *de;
                while ((de = readdir(dp)) != NULL) {
                    if (strstr(de->d_name, keywords[fi])) {
                        char candidate[MAX_PATH];
                        snprintf(candidate, sizeof(candidate), "%s/%s", pd->dir, de->d_name);
                        if (read_profile_file(&pd->files[fi], candidate, fi) == 0) {
                            pd->loaded[fi] = 1;
                            any_loaded = 1;
                            strncpy(matched_paths[fi], candidate, MAX_PATH - 1);
                            printf("Profile: matched '%s' -> %s\n", keywords[fi], de->d_name);
                            break;
                        }
                    }
                }
                closedir(dp);
            }
        }

        if (!any_loaded) {
            /* Show X11 error popup before exiting */
            Display *err_dpy = XOpenDisplay(NULL);
            if (err_dpy) {
                Widget err_top = XtAppInitialize(NULL, "PLTViewErr", NULL, 0,
                                                  &argc, argv, NULL, NULL, 0);
                Display *ed = XtDisplay(err_top);
                int es = DefaultScreen(ed);
                XFontStruct *ef = XLoadQueryFont(ed, "fixed");

                Widget err_form = XtCreateManagedWidget("form", formWidgetClass, err_top, NULL, 0);

                Arg ea[8]; int en;
                char errmsg[512];
                snprintf(errmsg, sizeof(errmsg),
                    "No profile files found in:\n%s\n\n"
                    "Expected files containing: surf, mean, flux, subgrid\n"
                    "(e.g. abl_surf, abl_mean, ...)\n\n"
                    "Please rename or use erf.data_log to set names.", argv[1]);

                en = 0;
                XtSetArg(ea[en], XtNlabel, errmsg); en++;
                XtSetArg(ea[en], XtNwidth, 420); en++;
                XtSetArg(ea[en], XtNborderWidth, 0); en++;
                Widget emsg = XtCreateManagedWidget("msg", labelWidgetClass, err_form, ea, en);

                en = 0;
                XtSetArg(ea[en], XtNfromVert, emsg); en++;
                XtSetArg(ea[en], XtNlabel, "OK"); en++;
                Widget eok = XtCreateManagedWidget("ok", commandWidgetClass, err_form, ea, en);
                (void)eok;

                XtRealizeWidget(err_top);
                /* Simple event loop: close on OK click or any key */
                XtAppContext ectx = XtWidgetToApplicationContext(err_top);
                XEvent ev;
                int done = 0;
                while (!done) {
                    XtAppNextEvent(ectx, &ev);
                    if (ev.type == ButtonPress) done = 1;
                    else XtDispatchEvent(&ev);
                }
                (void)ef; (void)es;
                XtDestroyWidget(err_top);
                XCloseDisplay(err_dpy);
            }
            fprintf(stderr, "Error: No profile files (surf/mean/flux/subgrid) found in %s\n", argv[1]);
            free(pd); return 1;
        }

        /* Set initial file to first loaded profile file (prefer mean) */
        profile_current_file = PROFILE_FILE_MEAN;
        if (!pd->loaded[profile_current_file]) {
            for (int fi = 0; fi < N_PROFILE_FILES; fi++) {
                if (pd->loaded[fi]) { profile_current_file = fi; break; }
            }
        }
        ProfileFile *init_pf = &pd->files[profile_current_file];
        profile_current_col = init_pf->has_z ? 2 : 1;
        if (profile_current_col >= init_pf->ncols) profile_current_col = 0;

        init_profile_gui(pd, argc, argv);
        update_profile_info_label(pd);
        render_profile_canvas(pd);

        printf("\n1D Profile Controls:\n");
        printf("  Click [surf/mean/flux/subgrid] to switch file\n");
        printf("  Click variable name button to select column\n");
        printf("  Click </> to step through timesteps\n");
        printf("  Click LogX to toggle log x-axis\n\n");

        XtAppContext app_context = XtWidgetToApplicationContext(toplevel);
        while (1) {
            XEvent event;
            XtAppNextEvent(app_context, &event);
            if (event.type == Expose) {
                if (event.xexpose.window == profile_canvas) {
                    render_profile_canvas(pd);
                    if (!initial_focus_set) {
                        XSetInputFocus(display, profile_canvas, RevertToParent, CurrentTime);
                        initial_focus_set = 1;
                    }
                }
            } else if (event.type == KeyPress &&
                       event.xkey.window == profile_canvas) {
                KeySym key = XLookupKeysym(&event.xkey, 0);
                ProfileFile *cpf = &pd->files[profile_current_file];
                if ((key == XK_Right || key == XK_n) && pd->loaded[profile_current_file] && cpf->has_z) {
                    profile_current_time_idx++;
                    if (profile_current_time_idx >= cpf->ntimes) profile_current_time_idx = 0;
                    update_profile_info_label(pd);
                    render_profile_canvas(pd);
                    continue;
                } else if ((key == XK_Left || key == XK_p) && pd->loaded[profile_current_file] && cpf->has_z) {
                    profile_current_time_idx--;
                    if (profile_current_time_idx < 0) profile_current_time_idx = cpf->ntimes - 1;
                    update_profile_info_label(pd);
                    render_profile_canvas(pd);
                    continue;
                }
            }
            XtDispatchEvent(&event);
        }

        for (int fi = 0; fi < N_PROFILE_FILES; fi++)
            free_profile_file(&pd->files[fi]);
        free(pd);
        return 0;
    }

    if (read_header(&pf) < 0) return 1;

    /* Initialize to first level */
    pf.current_level = 0;

    if (read_cell_h(&pf) < 0) return 1;

    /* Load first variable */
    pf.current_var = 0;
    pf.slice_axis = 2;  /* Z */
    pf.slice_idx = 0;  /* Start at first layer */
    pf.colormap = 0;  /* viridis */

    read_variable_data(&pf, 0);

    /* Initialize GUI */
    init_gui(&pf, argc, argv);

    update_layer_label(&pf);
    update_time_label();
    update_info_label(&pf);
    render_slice(&pf);

    printf("\nGUI Controls:\n");
    printf("  Click variable buttons to change variable\n");
    printf("  Click X/Y/Z buttons to switch axis\n");
    if (pf.n_levels > 1) {
        printf("  Click Level 0/Level 1/... buttons to switch level\n");
    }
    printf("  Click Colormap button to select colormap (or use keyboard 1-8)\n");
    printf("  Click v/^ buttons to navigate layers (or use keyboard Up/Down arrows)\n");
    if (n_timesteps > 1) {
        printf("  Click </> buttons to navigate timesteps (or use keyboard Left/Right arrows)\n");
    }
    printf("\n");
    
    /* Main event loop with expose and keyboard handling */
    XtAppContext app_context = XtWidgetToApplicationContext(toplevel);
    while (1) {
        XEvent event;
        XtAppNextEvent(app_context, &event);
        
        /* Handle expose events */
        if (event.type == Expose) {
            if (event.xexpose.window == canvas && global_pf && global_pf->data) {
                render_slice(global_pf);
                /* Set keyboard focus on first expose - needed for remote X11 */
                if (!initial_focus_set) {
                    XSetInputFocus(display, canvas, RevertToParent, CurrentTime);
                    initial_focus_set = 1;
                }
            } else if (event.xexpose.window == colorbar && global_pf) {
                draw_colorbar(current_vmin, current_vmax, global_pf->colormap,
                              global_pf->variables[global_pf->current_var]);
            }
        }
        /* Handle keyboard events */
        else if (event.type == KeyPress && global_pf) {
            /* Let Xaw text widgets handle input when a dialog is active */
            if (dialog_active) {
                XtDispatchEvent(&event);
                continue;
            }
            /* Only process keyboard shortcuts if the event is from the main canvas window */
            /* This prevents dialog text input from triggering colormap changes */
            if (event.xkey.window != canvas) {
                XtDispatchEvent(&event);
                continue;
            }
            
            KeySym key = XLookupKeysym(&event.xkey, 0);
            int max_idx = global_pf->grid_dims[global_pf->slice_axis] - 1;
            int changed = 0;
            
            if (key == XK_plus || key == XK_equal || key == XK_Up) {
                if (global_pf->slice_idx < max_idx) {
                    global_pf->slice_idx++;
                    changed = 1;
                }
            } else if (key == XK_minus || key == XK_underscore || key == XK_Down) {
                if (global_pf->slice_idx > 0) {
                    global_pf->slice_idx--;
                    changed = 1;
                }
            } else if (key >= XK_1 && key <= XK_8) {
                /* Switch colormap with 1-8 keys */
                global_pf->colormap = key - XK_1;
                changed = 1;
            } else if (key == XK_Right && n_timesteps > 1) {
                /* Next timestep */
                int new_timestep = current_timestep + 1;
                if (new_timestep >= n_timesteps) new_timestep = 0;
                switch_timestep(global_pf, new_timestep);
                continue;  /* switch_timestep handles all updates */
            } else if (key == XK_Left && n_timesteps > 1) {
                /* Previous timestep */
                int new_timestep = current_timestep - 1;
                if (new_timestep < 0) new_timestep = n_timesteps - 1;
                switch_timestep(global_pf, new_timestep);
                continue;  /* switch_timestep handles all updates */
            } else if (key == XK_r || key == XK_0) {
                /* Reset zoom */
                if (zoom_level > 1.0) {
                    zoom_reset();
                    changed = 1;
                }
            }

            if (changed) {
                update_layer_label(global_pf);
                update_info_label(global_pf);
                render_slice(global_pf);
                update_distribution_histogram(-1);  /* Auto-update distribution popup */
            }
        }

        XtDispatchEvent(&event);
    }

    cleanup(&pf);
    return 0;
}
