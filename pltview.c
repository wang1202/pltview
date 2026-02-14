/*
 * pltview.c - Fast AMReX plotfile viewer in C
 * Similar to ncview, using X11/Athena Widgets for GUI
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
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

#define MAX_VARS 128
#define MAX_BOXES 1024
#define MAX_PATH 512
#define MAX_LINE 1024
#define MAX_TIMESTEPS 1024
#define MAX_LEVELS 10

/* Data structures */
typedef struct {
    int lo[3];
    int hi[3];
    char filename[64];
} Box;

/* Per-level data storage for multi-level overlay rendering */
typedef struct {
    int grid_dims[3];       /* Grid dimensions for this level */
    int level_lo[3];        /* Lower index bounds in level's coordinates */
    int level_hi[3];        /* Upper index bounds in level's coordinates */
    Box boxes[MAX_BOXES];   /* Box definitions for this level */
    int n_boxes;            /* Number of boxes at this level */
    double *data;           /* Variable data for this level */
    int loaded;             /* Flag: 1 if data is loaded, 0 otherwise */
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
} PlotfileData;

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
    double cutoff_radius;   /* Cutoff in um (0 = no cutoff) */
    double custom_bin_width; /* Custom bin width in um (0 = auto/Sturges) */
    /* Per-grid info from particle Header */
    int n_grids;
    int grid_file_num[MAX_BOXES];
    int grid_count[MAX_BOXES];
    long grid_offset[MAX_BOXES];
} ParticleData;

/* X11 globals */
Display *display;
Widget toplevel, form, canvas_widget, var_box, info_label;
Widget axis_box, nav_box, colorbar_widget, layer_label;
Window canvas, colorbar;
GC gc, text_gc, colorbar_gc;
XImage *ximage;
int screen;
unsigned long *pixel_data;
int canvas_width = 800;
int canvas_height = 600;
Pixmap pixmap, colorbar_pixmap;
XFontStruct *font;
double current_vmin = 0, current_vmax = 1;

/* Current slice rendering info for mouse interaction */
double *current_slice_data = NULL;
int slice_width = 0, slice_height = 0;
int render_offset_x = 0, render_offset_y = 0;
int render_width = 0, render_height = 0;
char hover_value_text[256] = "";
int initial_focus_set = 0;  /* Flag for setting keyboard focus on first expose */
int dialog_active = 0;  /* Flag to track when a dialog is open */
Widget active_text_widget = NULL;  /* Text widget in active dialog */

/* Custom colorbar range */
int use_custom_range = 0;  /* Flag for using custom min/max */
double custom_vmin = 0.0;
double custom_vmax = 1.0;

/* Multi-timestep support */
char *timestep_paths[MAX_TIMESTEPS];  /* Array of plotfile paths */
int timestep_numbers[MAX_TIMESTEPS];   /* Numerical values for sorting */
int timestep_levels[MAX_TIMESTEPS];    /* Number of levels at each timestep */
int n_timesteps = 0;                   /* Number of timesteps found */
int current_timestep = 0;              /* Current timestep index */
int max_levels_all_timesteps = 1;      /* Max levels across all timesteps */
Widget time_label = NULL;              /* Time step display label */

/* Data structure for histogram expose handler (forward declaration for SDM) */
typedef struct {
    double *bin_counts;
    double *bin_centers;
    int n_bins;
    double count_max;
    double bin_min, bin_max;
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
int sdm_dialog_active = 0;
Widget sdm_active_text_widget = NULL;
int sdm_active_field = 0;  /* 0=cutoff, 1=binwidth */
Widget sdm_dialog_shell = NULL;
Widget overlay_button = NULL;  /* Overlay toggle button */
Widget map_dialog_shell = NULL;
Widget map_unavailable_shell = NULL;
int map_color_option = 0;  /* 0=black, 1=red, 2=gray, 3=white */
unsigned long map_color_pixel = 0;
int map_coastlines_enabled = 1;
double map_last_lon_min = 0.0, map_last_lon_max = 0.0, map_last_lat_min = 0.0, map_last_lat_max = 0.0;
int map_has_bounds = 0;
int map_auto_detected = 0;

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
void colorbar_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void init_gui(PlotfileData *pf, int argc, char **argv);
void render_slice(PlotfileData *pf);
void update_info_label(PlotfileData *pf);
void var_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void axis_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void map_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void show_map_settings_dialog(PlotfileData *pf);
void show_map_unavailable_dialog(void);
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
void update_layer_label(PlotfileData *pf);
void canvas_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void canvas_motion_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void canvas_button_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
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
void sdm_settings_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_settings_apply_callback(Widget w, XtPointer client_data, XtPointer call_data);
void sdm_settings_close_callback(Widget w, XtPointer client_data, XtPointer call_data);

/* Global pointer for callbacks */
PlotfileData *global_pf = NULL;
Widget map_button_widget = NULL;

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

/* Switch to a different timestep */
void switch_timestep(PlotfileData *pf, int new_timestep) {
    if (new_timestep < 0 || new_timestep >= n_timesteps) return;

    current_timestep = new_timestep;

    /* Update plotfile directory */
    strncpy(pf->plotfile_dir, timestep_paths[current_timestep], MAX_PATH - 1);

    /* Always free old overlay data before reading new timestep */
    free_all_levels(pf);

    /* Save overlay_mode before read_header (which resets it) */
    int saved_overlay_mode = pf->overlay_mode;

    /* Re-read header for new timestep */
    read_header(pf);

    /* Restore overlay_mode */
    pf->overlay_mode = saved_overlay_mode;

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
            /* Parse FabOnDisk: Cell_D_XXXXX */
            char *p = strchr(line, ':');
            if (p) {
                p++;
                while (*p == ' ') p++;
                char *end = strchr(p, ' ');
                if (end) *end = '\0';
                end = strchr(p, '\n');
                if (end) *end = '\0';
                strncpy(pf->boxes[pf->n_boxes].filename, p, 63);
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
    /* Initialize any remaining dimensions to 0 for 2D cases */
    for (i = pf->ndim; i < 3; i++) {
        pf->level_lo[i] = 0;
        pf->level_hi[i] = 0;
    }

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
        
        /* Skip FAB header (read until newline) */
        int c;
        while ((c = fgetc(fp)) != EOF && c != '\n');
        
        /* Skip to variable data */
        fseek(fp, var_idx * box_size * sizeof(double), SEEK_CUR);
        
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
            /* Parse FabOnDisk: Cell_D_XXXXX */
            char *p = strchr(line, ':');
            if (p) {
                p++;
                while (*p == ' ') p++;
                char *end = strchr(p, ' ');
                if (end) *end = '\0';
                end = strchr(p, '\n');
                if (end) *end = '\0';
                strncpy(ld->boxes[ld->n_boxes].filename, p, 63);
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
    for (i = pf->ndim; i < 3; i++) {
        ld->level_lo[i] = 0;
        ld->level_hi[i] = 0;
        ld->grid_dims[i] = 1;
    }

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

        /* Skip FAB header */
        int c;
        while ((c = fgetc(fp)) != EOF && c != '\n');

        /* Skip to variable data */
        fseek(fp, var_idx * box_size * sizeof(double), SEEK_CUR);

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
        pf->levels[level].loaded = 0;
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

/* Apply colormap to data */
void apply_colormap(double *data, int width, int height, 
                   unsigned long *pixels, double vmin, double vmax, int cmap_type) {
    int i, j;
    double range = vmax - vmin;
    if (range < 1e-10) range = 1.0;
    
    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            double val = data[j * width + i];
            double t = (val - vmin) / range;
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
        double value = vmin + fraction * (vmax - vmin);

        /* Map to drawable area with margins */
        int y = top_margin + (canvas_height - top_margin - bottom_margin) -
                (int)(fraction * (canvas_height - top_margin - bottom_margin));

        /* Draw tick mark */
        XDrawLine(display, colorbar, text_gc, width, y, width + 5, y);

        /* Draw label with vertical centering adjustment */
        snprintf(text, sizeof(text), "%.2e", value);
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
    
    /* Variable buttons box */
    n = 0;
    XtSetArg(args[n], XtNfromVert, info_label); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientVertical); n++;
    XtSetArg(args[n], XtNtop, XawChainTop); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    var_box = XtCreateManagedWidget("varBox", boxWidgetClass, form, args, n);
    
    /* Add variable buttons (all available variables) */
    for (i = 0; i < pf->n_vars; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, pf->variables[i]); n++;
        button = XtCreateManagedWidget(pf->variables[i], commandWidgetClass, var_box, args, n);
        XtAddCallback(button, XtNcallback, var_button_callback, (XtPointer)(long)i);
    }
    
    /* Canvas drawing area */
    n = 0;
    XtSetArg(args[n], XtNfromVert, info_label); n++;
    XtSetArg(args[n], XtNfromHoriz, var_box); n++;
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
    XtSetArg(args[n], XtNfromHoriz, var_box); n++;
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
    
    /* Map button */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Map: OFF"); n++;
    map_button_widget = XtCreateManagedWidget("Map", commandWidgetClass, axis_box, args, n);
    XtAddCallback(map_button_widget, XtNcallback, map_button_callback, NULL);

    /* COLUMN 2, ROW 2: Tools box (Colormap, Range, Distrib) */
    Widget tools_box;
    n = 0;
    XtSetArg(args[n], XtNfromVert, axis_box); n++;
    XtSetArg(args[n], XtNfromHoriz, nav_box); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    tools_box = XtCreateManagedWidget("toolsBox", boxWidgetClass, form, args, n);

    /* Colormap button - opens popup with colormap options */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Colormap"); n++;
    button = XtCreateManagedWidget("colormap", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, colormap_button_callback, NULL);

    /* Range button for custom colorbar min/max */
    n = 0;
    XtSetArg(args[n], XtNlabel, "Range"); n++;
    button = XtCreateManagedWidget("range", commandWidgetClass, tools_box, args, n);
    XtAddCallback(button, XtNcallback, range_button_callback, NULL);

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
        XtSetArg(args[n], XtNfromHoriz, var_box); n++;
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
    pixel_data = (unsigned long *)malloc(canvas_width * canvas_height * sizeof(unsigned long));
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
}

/* Update info label */
void update_info_label(PlotfileData *pf) {
    char text[512];
    const char *axis_names[] = {"X", "Y", "Z"};
    int max_idx = pf->grid_dims[pf->slice_axis];
    
    if (hover_value_text[0] != '\0') {
        if (pf->n_levels > 1) {
            snprintf(text, sizeof(text), 
                     "%s | Level: %d | Axis: %s | Layer: %d/%d | Time: %.3f | %s",
                     pf->variables[pf->current_var],
                     pf->current_level,
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time,
                     hover_value_text);
        } else {
            snprintf(text, sizeof(text), 
                     "%s | Axis: %s | Layer: %d/%d | Time: %.3f | %s",
                     pf->variables[pf->current_var],
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time,
                     hover_value_text);
        }
    } else {
        if (pf->n_levels > 1) {
            snprintf(text, sizeof(text), 
                     "%s | Level: %d | Axis: %s | Layer: %d/%d | Time: %.3f",
                     pf->variables[pf->current_var],
                     pf->current_level,
                     axis_names[pf->slice_axis],
                     pf->slice_idx + 1, max_idx,
                     pf->time);
        } else {
            snprintf(text, sizeof(text), 
                     "%s | Axis: %s | Layer: %d/%d | Time: %.3f",
                     pf->variables[pf->current_var],
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

        update_info_label(global_pf);
        render_slice(global_pf);
    }
}

/* Axis button callback */
void axis_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int axis = (int)(long)client_data;
    if (global_pf) {
        global_pf->slice_axis = axis;
        global_pf->slice_idx = 0;  /* Start at first layer */

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
        update_info_label(global_pf);
        render_slice(global_pf);
    }
}

/* Map button callback */
void map_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        /* Check if lon_m and lat_m variables exist */
        int lon_idx = find_variable_index(global_pf, "lon_m");
        int lat_idx = find_variable_index(global_pf, "lat_m");
        
        if (lon_idx >= 0 && lat_idx >= 0) {
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
            show_map_unavailable_dialog();
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

static void scan_coastline_files(void) {
    if (n_coastlines > 0) return;

    DIR *dir = opendir("map_layers");
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
        snprintf(ce->filename, sizeof(ce->filename), "map_layers/%s", name);

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
    map_coastlines_enabled = 0;
    for (int i = 0; i < n_coastlines; i++) {
        coastlines[i].enabled = 0;
        update_coastline_button_label(&coastlines[i]);
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

void show_map_unavailable_dialog(void) {
    if (map_unavailable_shell) {
        XtPopup(map_unavailable_shell, XtGrabNone);
        return;
    }

    map_unavailable_shell = XtVaCreatePopupShell(
        "mapUnavailable", transientShellWidgetClass, toplevel,
        XtNtitle, "Map",
        NULL);

    Widget msg_form = XtVaCreateManagedWidget("mapUnavailableForm", formWidgetClass,
                                              map_unavailable_shell, NULL);
    Widget msg_label = XtVaCreateManagedWidget("mapUnavailableLabel", labelWidgetClass,
                                               msg_form,
                                               XtNlabel, "lat_m and lon_m are not available",
                                               XtNborderWidth, 0,
                                               NULL);
    Widget ok_btn = XtVaCreateManagedWidget("mapUnavailableOk", commandWidgetClass, msg_form,
                                            XtNlabel, "OK",
                                            XtNfromVert, msg_label,
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
    update_info_label(global_pf);
    render_slice(global_pf);
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

/* Structure for range dialog */
typedef struct {
    Widget min_text;
    Widget max_text;
    Widget dialog_shell;
} RangeDialogData;

/* Global pointer to range dialog data for keyboard input */
RangeDialogData *active_range_dialog = NULL;
int active_field = 0;  /* 0 = min, 1 = max */

/* Apply custom range callback */
void range_apply_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    RangeDialogData *data = (RangeDialogData *)client_data;

    if (data) {
        String min_str, max_str;
        Arg args[1];

        XtSetArg(args[0], XtNstring, &min_str);
        XtGetValues(data->min_text, args, 1);
        XtSetArg(args[0], XtNstring, &max_str);
        XtGetValues(data->max_text, args, 1);

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
        active_range_dialog = NULL;
    }
}

/* Auto range callback - reset to data-driven min/max */
void range_auto_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    RangeDialogData *data = (RangeDialogData *)client_data;

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
        active_range_dialog = NULL;
    }
}

/* Close range dialog callback */
void range_close_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    RangeDialogData *data = (RangeDialogData *)client_data;

    if (data) {
        XtPopdown(data->dialog_shell);
        XtDestroyWidget(data->dialog_shell);
        free(data);
        dialog_active = 0;
        active_text_widget = NULL;
        active_range_dialog = NULL;
    }
}

/* Focus on min field */
void range_min_focus_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    RangeDialogData *data = (RangeDialogData *)client_data;
    if (data) {
        active_text_widget = data->min_text;
        active_field = 0;
    }
}

/* Focus on max field */
void range_max_focus_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    RangeDialogData *data = (RangeDialogData *)client_data;
    if (data) {
        active_text_widget = data->max_text;
        active_field = 1;
    }
}

/* Ensure focus follows mouse clicks in range dialog text fields */
void range_text_click_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    RangeDialogData *data = (RangeDialogData *)client_data;
    if (!data || event->type != ButtonPress) return;

    XtSetKeyboardFocus(data->dialog_shell, w);
    Time time = CurrentTime;
    XtCallAcceptFocus(w, &time);

    active_text_widget = w;
    active_field = (w == data->max_text) ? 1 : 0;
}

/* Range button callback - dialog to set custom min/max */
void range_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        Arg args[10];
        int n;
        Widget dialog_shell, form, label, button, min_text, max_text;
        char min_str[64], max_str[64];

        /* Format current values */
        snprintf(min_str, sizeof(min_str), "%.6e", use_custom_range ? custom_vmin : current_vmin);
        snprintf(max_str, sizeof(max_str), "%.6e", use_custom_range ? custom_vmax : current_vmax);

        n = 0;
        XtSetArg(args[n], XtNtitle, "Set Colorbar Range"); n++;
        dialog_shell = XtCreatePopupShell("rangeDialog", transientShellWidgetClass, toplevel, args, n);

        n = 0;
        form = XtCreateManagedWidget("form", formWidgetClass, dialog_shell, args, n);

        /* Title label */
        n = 0;
        XtSetArg(args[n], XtNlabel, "Set colorbar min/max values:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        label = XtCreateManagedWidget("title", labelWidgetClass, form, args, n);

        /* Min label */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNlabel, "Min:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        Widget min_label = XtCreateManagedWidget("minLabel", labelWidgetClass, form, args, n);

        /* Min text input */
        n = 0;
        XtSetArg(args[n], XtNfromVert, label); n++;
        XtSetArg(args[n], XtNfromHoriz, min_label); n++;
        XtSetArg(args[n], XtNwidth, 150); n++;
        XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
        XtSetArg(args[n], XtNstring, min_str); n++;
        min_text = XtCreateManagedWidget("minInput", asciiTextWidgetClass, form, args, n);

        /* Max label */
        n = 0;
        XtSetArg(args[n], XtNfromVert, min_label); n++;
        XtSetArg(args[n], XtNlabel, "Max:"); n++;
        XtSetArg(args[n], XtNborderWidth, 0); n++;
        Widget max_label = XtCreateManagedWidget("maxLabel", labelWidgetClass, form, args, n);

        /* Max text input */
        n = 0;
        XtSetArg(args[n], XtNfromVert, min_label); n++;
        XtSetArg(args[n], XtNfromHoriz, max_label); n++;
        XtSetArg(args[n], XtNwidth, 150); n++;
        XtSetArg(args[n], XtNeditType, XawtextEdit); n++;
        XtSetArg(args[n], XtNstring, max_str); n++;
        max_text = XtCreateManagedWidget("maxInput", asciiTextWidgetClass, form, args, n);

        /* Create data structure */
        RangeDialogData *range_data = malloc(sizeof(RangeDialogData));
        range_data->min_text = min_text;
        range_data->max_text = max_text;
        range_data->dialog_shell = dialog_shell;

        /* Click-to-focus handlers for text fields */
        XtAddEventHandler(min_text, ButtonPressMask, False, range_text_click_handler, (XtPointer)range_data);
        XtAddEventHandler(max_text, ButtonPressMask, False, range_text_click_handler, (XtPointer)range_data);

        /* Apply button */
        n = 0;
        XtSetArg(args[n], XtNfromVert, max_label); n++;
        XtSetArg(args[n], XtNlabel, "Apply"); n++;
        button = XtCreateManagedWidget("apply", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, range_apply_callback, (XtPointer)range_data);

        /* Auto button */
        n = 0;
        XtSetArg(args[n], XtNfromVert, max_label); n++;
        XtSetArg(args[n], XtNfromHoriz, button); n++;
        XtSetArg(args[n], XtNlabel, "Auto"); n++;
        button = XtCreateManagedWidget("auto", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, range_auto_callback, (XtPointer)range_data);

        /* Close button */
        n = 0;
        XtSetArg(args[n], XtNfromVert, max_label); n++;
        XtSetArg(args[n], XtNfromHoriz, button); n++;
        XtSetArg(args[n], XtNlabel, "Close"); n++;
        button = XtCreateManagedWidget("close", commandWidgetClass, form, args, n);
        XtAddCallback(button, XtNcallback, range_close_callback, (XtPointer)range_data);

        XtRealizeWidget(dialog_shell);
        XtPopup(dialog_shell, XtGrabExclusive);

        /* Set keyboard focus to min text input */
        XtSetKeyboardFocus(dialog_shell, min_text);
        XSync(display, False);
        Time time = CurrentTime;
        XtCallAcceptFocus(min_text, &time);

        dialog_active = 1;
        active_text_widget = min_text;
        active_range_dialog = range_data;
        active_field = 0;
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
void render_slice(PlotfileData *pf) {
    int width, height;
    double *slice;
    double vmin = 1e30, vmax = -1e30;
    int i, j;
    char stats_text[128];

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
    extract_slice(pf, slice, pf->slice_axis, pf->slice_idx);

    /* Physical coordinate ranges for axes */
    double phys_xmin, phys_xmax, phys_ymin, phys_ymax;
    
    if (pf->map_mode) {
        /* Map mode: phys_xmin/max and phys_ymin/max will be set in the map rendering section */
    } else {
        /* Normal mode: use physical coordinates */
        phys_xmin = pf->prob_lo[x_axis];
        phys_xmax = pf->prob_hi[x_axis];
        phys_ymin = pf->prob_lo[y_axis];
        phys_ymax = pf->prob_hi[y_axis];
    }

    /* Store current slice for mouse interaction */
    if (current_slice_data) free(current_slice_data);
    current_slice_data = (double *)malloc(width * height * sizeof(double));
    memcpy(current_slice_data, slice, width * height * sizeof(double));
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

    /* Find data min/max, skipping gap cells when mask is active */
    for (i = 0; i < width * height; i++) {
        if (base_in_box && !base_in_box[i]) continue;
        if (slice[i] < vmin) vmin = slice[i];
        if (slice[i] > vmax) vmax = slice[i];
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
                
                /* Generate Z coordinates for this slice */
                for (j = 0; j < height; j++) {
                    for (i = 0; i < width; i++) {
                        int idx = j * width + i;
                        /* Map grid index to physical Z coordinate */
                        double z_coord = pf->prob_lo[2] + (j + 0.5) * (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                        y_coord_slice[idx] = z_coord;
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
                
                /* Generate Z coordinates for this slice */
                for (j = 0; j < height; j++) {
                    for (i = 0; i < width; i++) {
                        int idx = j * width + i;
                        /* Map grid index to physical Z coordinate */
                        double z_coord = pf->prob_lo[2] + (j + 0.5) * (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                        y_coord_slice[idx] = z_coord;
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
            
            /* Create pixel data for individual points */
            unsigned long *point_pixels = (unsigned long *)malloc(width * height * sizeof(unsigned long));
            apply_colormap(slice, width, height, point_pixels, display_vmin, display_vmax, pf->colormap);
            
            /* Render each data point at its coordinate */
            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    int idx = j * width + i;
                    double x_coord = x_geo_extent[idx];
                    double y_coord = y_coord_extent[idx];
                    
                    /* Map coordinates to screen coordinates */
                    if (x_coord >= phys_xmin && x_coord <= phys_xmax && y_coord >= phys_ymin && y_coord <= phys_ymax) {
                        int screen_x = offset_x + (int)((x_coord - phys_xmin) / (phys_xmax - phys_xmin) * local_render_width);
                        int screen_y = offset_y + (int)((phys_ymax - y_coord) / (phys_ymax - phys_ymin) * local_render_height);
                        
                        /* Draw a small rectangle for each data point */
                        XSetForeground(display, gc, point_pixels[idx]);
                        XFillRectangle(display, canvas, gc, screen_x - 1, screen_y - 1, 3, 3);
                    }
                }
            }
            
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

            double pixel_width = (double)local_render_width / width;
            double pixel_height = (double)local_render_height / height;

            for (j = 0; j < height; j++) {
                for (i = 0; i < width; i++) {
                    if (base_in_box && !base_in_box[j * width + i]) continue;

                    unsigned long pixel = pixel_data[j * width + i];
                    XSetForeground(display, gc, pixel);

                    int x = offset_x + (int)(i * pixel_width);
                    int flipped_j = height - 1 - j;
                    int y = offset_y + (int)(flipped_j * pixel_height);
                    int w = (int)((i + 1) * pixel_width) - (int)(i * pixel_width);
                    int h = (int)((flipped_j + 1) * pixel_height) - (int)(flipped_j * pixel_height);
                    if (w < 1) w = 1;
                    if (h < 1) h = 1;

                    XFillRectangle(display, canvas, gc, x, y, w, h);
                }
            }
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

        /* Draw pixels as filled rectangles with correct aspect ratio */
        double pixel_width = (double)local_render_width / width;
        double pixel_height = (double)local_render_height / height;

        for (j = 0; j < height; j++) {
            for (i = 0; i < width; i++) {
                if (base_in_box && !base_in_box[j * width + i]) continue;

                unsigned long pixel = pixel_data[j * width + i];
                XSetForeground(display, gc, pixel);

                int x = offset_x + (int)(i * pixel_width);
                /* Flip y-axis: higher j (higher y in data) should be at top of screen */
                int flipped_j = height - 1 - j;
                int y = offset_y + (int)(flipped_j * pixel_height);
                int w = (int)((i + 1) * pixel_width) - (int)(i * pixel_width);
                int h = (int)((flipped_j + 1) * pixel_height) - (int)(flipped_j * pixel_height);
                if (w < 1) w = 1;
                if (h < 1) h = 1;

                XFillRectangle(display, canvas, gc, x, y, w, h);
            }
        }
    }

    /* Store rendering parameters for mouse interaction */
    render_offset_x = offset_x;
    render_offset_y = offset_y;
    render_width = local_render_width;
    render_height = local_render_height;

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

            /* Extract slice from this level */
            double *level_slice = (double *)malloc(lwidth * lheight * sizeof(double));
            extract_slice_level(ld, level_slice, pf->slice_axis, level_slice_idx);

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

            /* Draw level pixels, skipping cells not inside any box */
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
            free(level_pixels);

            printf("Overlay level %d: slice %d, screen [%d,%d]-[%d,%d]\n",
                   level, level_slice_idx, screen_x0, screen_y0, screen_x1, screen_y1);
        }
    }

    /* Draw axis frame (border around data) */
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    XDrawRectangle(display, canvas, text_gc, offset_x, offset_y, local_render_width, local_render_height);

    /* Draw X-axis ticks and labels */
    int n_xticks = 5;
    char label[32];
    for (i = 0; i <= n_xticks; i++) {
        double frac = (double)i / n_xticks;
        int tick_x = offset_x + (int)(frac * local_render_width);
        double phys_val = phys_xmin + frac * (phys_xmax - phys_xmin);

        /* Draw tick mark */
        XDrawLine(display, canvas, text_gc, tick_x, offset_y + local_render_height,
              tick_x, offset_y + local_render_height + 5);

        /* Draw label */
        snprintf(label, sizeof(label), "%.3g", phys_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(display, canvas, text_gc, tick_x - label_width / 2,
                offset_y + local_render_height + 18, label, strlen(label));
    }

    /* Draw Y-axis ticks and labels */
    int n_yticks = 5;
    for (i = 0; i <= n_yticks; i++) {
        double frac = (double)i / n_yticks;
        int tick_y = offset_y + local_render_height - (int)(frac * local_render_height);
        double phys_val = phys_ymin + frac * (phys_ymax - phys_ymin);

        /* Draw tick mark */
        XDrawLine(display, canvas, text_gc, offset_x - 5, tick_y, offset_x, tick_y);

        /* Draw label */
        snprintf(label, sizeof(label), "%.3g", phys_val);
        int label_width = XTextWidth(font, label, strlen(label));
        XDrawString(display, canvas, text_gc, offset_x - label_width - 8,
                    tick_y + 4, label, strlen(label));
    }

    /* Draw axis labels with units */
    const char *axis_names[] = {"X", "Y", "Z"};
    char x_label[32], y_label[32];

    if (pf->map_mode) {
        /* Map mode: use longitude/latitude labels */
        strcpy(x_label, "Longitude (deg)");
        strcpy(y_label, "Latitude (deg)");
    } else {
        /* Normal mode: use physical coordinates with units */
        const char *unit_str = "(m)";
        snprintf(x_label, sizeof(x_label), "%s %s", axis_names[x_axis], unit_str);
        snprintf(y_label, sizeof(y_label), "%s %s", axis_names[y_axis], unit_str);
    }

    /* X-axis label (centered below ticks) */
    int xlabel_width = XTextWidth(font, x_label, strlen(x_label));
    XDrawString(display, canvas, text_gc,
                offset_x + local_render_width / 2 - xlabel_width / 2,
                offset_y + local_render_height + 35, x_label, strlen(x_label));

    /* Y-axis label (rotated text is hard in X11, so just draw at left) */
    XDrawString(display, canvas, text_gc, 5,
                offset_y + local_render_height / 2 + 4, y_label, strlen(y_label));

    /* Draw text overlay - show display range (custom if set) */
    if (use_custom_range) {
        snprintf(stats_text, sizeof(stats_text), "range: %.3e to %.3e (custom)", display_vmin, display_vmax);
    } else {
        snprintf(stats_text, sizeof(stats_text), "min: %.3e  max: %.3e", display_vmin, display_vmax);
    }
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    XSetBackground(display, text_gc, WhitePixel(display, screen));
    XDrawImageString(display, canvas, text_gc, left_margin, canvas_height - 5,
                    stats_text, strlen(stats_text));

    /* Draw colorbar */
    draw_colorbar(display_vmin, display_vmax, pf->colormap,
                  pf->variables[pf->current_var]);

    /* Draw quiver overlay if enabled */
    if (quiver_data.enabled) {
        render_quiver_overlay(pf);
    }
    
    /* Draw map overlay if in map mode */
    if (pf->map_mode) {
        render_map_overlay(pf, phys_xmin, phys_xmax, phys_ymin, phys_ymax);
    }

    XFlush(display);
    
    printf("Rendered: %s, slice %d/%d (%.3e to %.3e)\n", 
           pf->variables[pf->current_var], pf->slice_idx + 1,
           pf->grid_dims[pf->slice_axis], vmin, vmax);
    
    free(slice);
    if (base_in_box) free(base_in_box);
}

/* Mouse motion handler - show value at cursor */
void canvas_motion_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch) {
    if (!global_pf || !current_slice_data) return;
    
    int mouse_x = event->xmotion.x;
    int mouse_y = event->xmotion.y;
    
    /* Convert mouse coordinates to data coordinates */
    if (mouse_x < render_offset_x || mouse_x >= render_offset_x + render_width ||
        mouse_y < render_offset_y || mouse_y >= render_offset_y + render_height) {
        /* Outside data region - clear hover text */
        if (hover_value_text[0] != '\0') {
            hover_value_text[0] = '\0';
            update_info_label(global_pf);
        }
        return;
    }
    
    int data_x = (int)((mouse_x - render_offset_x) * slice_width / (double)render_width);
    /* Flip y-axis for mouse coordinates to match flipped rendering */
    int data_y = slice_height - 1 - (int)((mouse_y - render_offset_y) * slice_height / (double)render_height);
    
    if (data_x >= 0 && data_x < slice_width && data_y >= 0 && data_y < slice_height) {
        double value = current_slice_data[data_y * slice_width + data_x];
        
        /* Update hover value text and info label */
        snprintf(hover_value_text, sizeof(hover_value_text), "[%d,%d]: %.6e", data_x, data_y, value);
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

    if (event->xbutton.button != Button1) return;
    
    int mouse_x = event->xbutton.x;
    int mouse_y = event->xbutton.y;
    
    /* Convert mouse coordinates to data coordinates */
    if (mouse_x < render_offset_x || mouse_x >= render_offset_x + render_width ||
        mouse_y < render_offset_y || mouse_y >= render_offset_y + render_height) {
        return;
    }
    
    int data_x = (int)((mouse_x - render_offset_x) * slice_width / (double)render_width);
    /* Flip y-axis for mouse coordinates to match flipped rendering */
    int data_y = slice_height - 1 - (int)((mouse_y - render_offset_y) * slice_height / (double)render_height);
    
    if (data_x >= 0 && data_x < slice_width && data_y >= 0 && data_y < slice_height) {
        show_line_profiles(global_pf, data_x, data_y);
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

        /* Draw label */
        snprintf(label, sizeof(label), "%.0f", y_val);
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

/* Callback to destroy popup and free data */
void close_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    PopupData *popup_data = (PopupData *)client_data;
    
    if (popup_data) {
        /* Free plot data */
        for (int i = 0; i < 3; i++) {
            if (popup_data->plot_data_array[i]) {
                if (popup_data->plot_data_array[i]->data) 
                    free(popup_data->plot_data_array[i]->data);
                if (popup_data->plot_data_array[i]->x_values)
                    free(popup_data->plot_data_array[i]->x_values);
                free(popup_data->plot_data_array[i]);
            }
        }
        
        /* Destroy the popup shell */
        XtDestroyWidget(popup_data->shell);
        
        /* Free the popup data structure */
        free(popup_data);
    }
}

/* Show 1D line profiles through clicked point along x, y, z */
void show_line_profiles(PlotfileData *pf, int data_x, int data_y) {
    /* Get 3D coordinates based on current slice */
    int x_coord, y_coord, z_coord;
    
    if (pf->slice_axis == 2) {  /* Z slice */
        x_coord = data_x;
        y_coord = data_y;
        z_coord = pf->slice_idx;
    } else if (pf->slice_axis == 1) {  /* Y slice */
        x_coord = data_x;
        y_coord = pf->slice_idx;
        z_coord = data_y;
    } else {  /* X slice */
        x_coord = pf->slice_idx;
        y_coord = data_x;
        z_coord = data_y;
    }
    
    /* Create plot data structures */
    PlotData *x_plot_data = (PlotData *)malloc(sizeof(PlotData));
    PlotData *y_plot_data = (PlotData *)malloc(sizeof(PlotData));
    PlotData *z_plot_data = (PlotData *)malloc(sizeof(PlotData));
    
    /* Extract and store X profile data */
    x_plot_data->n_points = pf->grid_dims[0];
    x_plot_data->data = (double *)malloc(pf->grid_dims[0] * sizeof(double));
    x_plot_data->x_values = (double *)malloc(pf->grid_dims[0] * sizeof(double));
    x_plot_data->vmin = 1e30;
    x_plot_data->vmax = -1e30;
    for (int i = 0; i < pf->grid_dims[0]; i++) {
        x_plot_data->x_values[i] = i;
        int idx = z_coord * pf->grid_dims[0] * pf->grid_dims[1] + y_coord * pf->grid_dims[0] + i;
        x_plot_data->data[i] = pf->data[idx];
        if (x_plot_data->data[i] < x_plot_data->vmin) x_plot_data->vmin = x_plot_data->data[i];
        if (x_plot_data->data[i] > x_plot_data->vmax) x_plot_data->vmax = x_plot_data->data[i];
    }
    x_plot_data->xmin = 0;
    x_plot_data->xmax = pf->grid_dims[0] - 1;
    snprintf(x_plot_data->title, sizeof(x_plot_data->title), "%s along X (Y=%d, Z=%d)", 
             pf->variables[pf->current_var], y_coord, z_coord);
    snprintf(x_plot_data->xlabel, sizeof(x_plot_data->xlabel), "X");
    
    /* Extract and store Y profile data */
    y_plot_data->n_points = pf->grid_dims[1];
    y_plot_data->data = (double *)malloc(pf->grid_dims[1] * sizeof(double));
    y_plot_data->x_values = (double *)malloc(pf->grid_dims[1] * sizeof(double));
    y_plot_data->vmin = 1e30;
    y_plot_data->vmax = -1e30;
    for (int j = 0; j < pf->grid_dims[1]; j++) {
        y_plot_data->x_values[j] = j;
        int idx = z_coord * pf->grid_dims[0] * pf->grid_dims[1] + j * pf->grid_dims[0] + x_coord;
        y_plot_data->data[j] = pf->data[idx];
        if (y_plot_data->data[j] < y_plot_data->vmin) y_plot_data->vmin = y_plot_data->data[j];
        if (y_plot_data->data[j] > y_plot_data->vmax) y_plot_data->vmax = y_plot_data->data[j];
    }
    y_plot_data->xmin = 0;
    y_plot_data->xmax = pf->grid_dims[1] - 1;
    snprintf(y_plot_data->title, sizeof(y_plot_data->title), "%s along Y (X=%d, Z=%d)", 
             pf->variables[pf->current_var], x_coord, z_coord);
    snprintf(y_plot_data->xlabel, sizeof(y_plot_data->xlabel), "Y");
    
    /* Extract and store Z profile data */
    z_plot_data->n_points = pf->grid_dims[2];
    z_plot_data->data = (double *)malloc(pf->grid_dims[2] * sizeof(double));
    z_plot_data->x_values = (double *)malloc(pf->grid_dims[2] * sizeof(double));
    z_plot_data->vmin = 1e30;
    z_plot_data->vmax = -1e30;
    for (int k = 0; k < pf->grid_dims[2]; k++) {
        z_plot_data->x_values[k] = k;
        int idx = k * pf->grid_dims[0] * pf->grid_dims[1] + y_coord * pf->grid_dims[0] + x_coord;
        z_plot_data->data[k] = pf->data[idx];
        if (z_plot_data->data[k] < z_plot_data->vmin) z_plot_data->vmin = z_plot_data->data[k];
        if (z_plot_data->data[k] > z_plot_data->vmax) z_plot_data->vmax = z_plot_data->data[k];
    }
    z_plot_data->xmin = 0;
    z_plot_data->xmax = pf->grid_dims[2] - 1;
    snprintf(z_plot_data->title, sizeof(z_plot_data->title), "%s along Z (X=%d, Y=%d)", 
             pf->variables[pf->current_var], x_coord, y_coord);
    snprintf(z_plot_data->xlabel, sizeof(z_plot_data->xlabel), "Z");
    
    /* Create popup data structure */
    PopupData *popup_data = (PopupData *)malloc(sizeof(PopupData));
    popup_data->plot_data_array[0] = x_plot_data;
    popup_data->plot_data_array[1] = y_plot_data;
    popup_data->plot_data_array[2] = z_plot_data;
    
    /* Create popup shell */
    Widget popup_shell = XtVaCreatePopupShell("Line Profiles",
        transientShellWidgetClass, toplevel,
        XtNwidth, 900,
        XtNheight, 700,
        NULL);
    
    /* Store shell in popup data */
    popup_data->shell = popup_shell;
    
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
    
    /* Create three plot canvases with expose event handlers */
    Widget x_canvas = XtVaCreateManagedWidget("x_plot",
        simpleWidgetClass, popup_form,
        XtNfromVert, title_label,
        XtNwidth, 880,
        XtNheight, 180,
        XtNborderWidth, 1,
        NULL);
    
    Widget y_canvas = XtVaCreateManagedWidget("y_plot",
        simpleWidgetClass, popup_form,
        XtNfromVert, x_canvas,
        XtNwidth, 880,
        XtNheight, 180,
        XtNborderWidth, 1,
        NULL);
    
    Widget z_canvas = XtVaCreateManagedWidget("z_plot",
        simpleWidgetClass, popup_form,
        XtNfromVert, y_canvas,
        XtNwidth, 880,
        XtNheight, 180,
        XtNborderWidth, 1,
        NULL);
    
    /* Add expose event handlers to all three canvases */
    XtAddEventHandler(x_canvas, ExposureMask, False, plot_expose_handler, x_plot_data);
    XtAddEventHandler(y_canvas, ExposureMask, False, plot_expose_handler, y_plot_data);
    XtAddEventHandler(z_canvas, ExposureMask, False, plot_expose_handler, z_plot_data);
    
    /* Close button */
    Widget close_button = XtVaCreateManagedWidget("Close",
        commandWidgetClass, popup_form,
        XtNfromVert, z_canvas,
        NULL);
    
    XtAddCallback(close_button, XtNcallback, close_popup_callback, popup_data);
    
    /* Show popup */
    XtPopup(popup_shell, XtGrabNone);
}

/* Popup data for profile (3 plots) */
typedef struct {
    Widget shell;
    PlotData *mean_plot;
    PlotData *std_plot;
    PlotData *skewness_plot;
} ProfilePopupData;

/* Callback to destroy profile popup and free data */
void close_profile_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    ProfilePopupData *popup_data = (ProfilePopupData *)client_data;

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

    /* Create plot data for mean */
    PlotData *mean_plot = (PlotData *)malloc(sizeof(PlotData));
    mean_plot->n_points = n_slices;
    mean_plot->data = means;
    mean_plot->x_values = (double *)malloc(n_slices * sizeof(double));
    memcpy(mean_plot->x_values, layer_indices, n_slices * sizeof(double));
    mean_plot->vmin = 1e30;
    mean_plot->vmax = -1e30;
    for (int i = 0; i < n_slices; i++) {
        if (means[i] < mean_plot->vmin) mean_plot->vmin = means[i];
        if (means[i] > mean_plot->vmax) mean_plot->vmax = means[i];
    }
    mean_plot->xmin = 1;
    mean_plot->xmax = n_slices;
    snprintf(mean_plot->title, sizeof(mean_plot->title), "%s Mean along %s axis",
             pf->variables[pf->current_var], axis_names[axis]);
    snprintf(mean_plot->xlabel, sizeof(mean_plot->xlabel), "%s Layer", axis_names[axis]);
    snprintf(mean_plot->vlabel, sizeof(mean_plot->vlabel), "%s Mean", pf->variables[pf->current_var]);

    /* Create plot data for std */
    PlotData *std_plot = (PlotData *)malloc(sizeof(PlotData));
    std_plot->n_points = n_slices;
    std_plot->data = stds;
    std_plot->x_values = layer_indices;  /* Share with mean, will be freed once */
    std_plot->vmin = 1e30;
    std_plot->vmax = -1e30;
    for (int i = 0; i < n_slices; i++) {
        if (stds[i] < std_plot->vmin) std_plot->vmin = stds[i];
        if (stds[i] > std_plot->vmax) std_plot->vmax = stds[i];
    }
    std_plot->xmin = 1;
    std_plot->xmax = n_slices;
    snprintf(std_plot->title, sizeof(std_plot->title), "%s Std Dev along %s axis",
             pf->variables[pf->current_var], axis_names[axis]);
    snprintf(std_plot->xlabel, sizeof(std_plot->xlabel), "%s Layer", axis_names[axis]);
    snprintf(std_plot->vlabel, sizeof(std_plot->vlabel), "%s Std", pf->variables[pf->current_var]);

    /* Create plot data for skewness */
    PlotData *skewness_plot = (PlotData *)malloc(sizeof(PlotData));
    skewness_plot->n_points = n_slices;
    skewness_plot->data = skewness;
    skewness_plot->x_values = layer_indices;  /* Share with mean, will be freed once */
    skewness_plot->vmin = 1e30;
    skewness_plot->vmax = -1e30;
    for (int i = 0; i < n_slices; i++) {
        if (skewness[i] < skewness_plot->vmin) skewness_plot->vmin = skewness[i];
        if (skewness[i] > skewness_plot->vmax) skewness_plot->vmax = skewness[i];
    }
    skewness_plot->xmin = 1;
    skewness_plot->xmax = n_slices;
    snprintf(skewness_plot->title, sizeof(skewness_plot->title), "%s Skewness along %s axis",
             pf->variables[pf->current_var], axis_names[axis]);
    snprintf(skewness_plot->xlabel, sizeof(skewness_plot->xlabel), "%s Layer", axis_names[axis]);
    snprintf(skewness_plot->vlabel, sizeof(skewness_plot->vlabel), "%s Skewness", pf->variables[pf->current_var]);

    /* Create popup data structure */
    ProfilePopupData *popup_data = (ProfilePopupData *)malloc(sizeof(ProfilePopupData));
    popup_data->mean_plot = mean_plot;
    popup_data->std_plot = std_plot;
    popup_data->skewness_plot = skewness_plot;

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

    /* Show popup */
    XtPopup(popup_shell, XtGrabNone);
}

/* Profile button callback */
void profile_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data) {
        show_slice_statistics(global_pf);
    }
}

/* Popup data for distribution histogram */
typedef struct {
    Widget shell;
    double *bin_counts;
    double *bin_centers;
    int n_bins;
} DistributionPopupData;

/* Callback to destroy distribution popup and free data */
void close_distribution_popup_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    DistributionPopupData *popup_data = (DistributionPopupData *)client_data;

    if (popup_data) {
        if (popup_data->bin_counts) free(popup_data->bin_counts);
        if (popup_data->bin_centers) free(popup_data->bin_centers);
        XtDestroyWidget(popup_data->shell);
        free(popup_data);
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
            if (y_val >= 1e6 || (y_val != 0 && y_val < 0.01))
                snprintf(label, sizeof(label), "%.1e", y_val);
            else
                snprintf(label, sizeof(label), "%.0f", y_val);
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
        double bc = hist->bin_centers[i];
        double bcount = hist->bin_counts[i];
        if (bcount <= 0 && log_y) continue;

        int bar_x, bar_w, bar_h, bar_y;

        if (log_x) {
            double left_edge = bc - bin_width / 2;
            double right_edge = bc + bin_width / 2;
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
            bar_x = plot_left + (int)((bc - x_min - bin_width / 2) / (x_max - x_min) * plot_width);
            bar_w = (int)(bin_width / (x_max - x_min) * plot_width);
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

/* Show distribution histogram for current slice */
void show_distribution(PlotfileData *pf) {
    const char *axis_names[] = {"X", "Y", "Z"};
    int axis = pf->slice_axis;
    int slice_idx = pf->slice_idx;

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

    /* Extract slice data and calculate statistics */
    double *slice_data = (double *)malloc(slice_size * sizeof(double));
    double sum = 0.0, sum_sq = 0.0;
    double data_min = 1e30, data_max = -1e30;

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
            double val = pf->data[idx];
            slice_data[k++] = val;
            sum += val;
            sum_sq += val * val;
            if (val < data_min) data_min = val;
            if (val > data_max) data_max = val;
        }
    }

    /* Calculate mean and std */
    double mean = sum / slice_size;
    double variance = (sum_sq / slice_size) - (mean * mean);
    double std = (variance > 0) ? sqrt(variance) : 0.0;

    /* Calculate skewness (second pass) */
    double sum_third = 0.0;
    for (int i = 0; i < slice_size; i++) {
        double diff = slice_data[i] - mean;
        sum_third += diff * diff * diff;
    }
    double skewness = 0.0;
    if (std > 0) {
        double std3 = std * std * std;
        skewness = (sum_third / slice_size) / std3;
    }

    /* Determine number of bins using Sturges' rule */
    int n_bins = (int)(1 + 3.322 * log10((double)slice_size));
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
    for (int i = 0; i < slice_size; i++) {
        int bin = (int)((slice_data[i] - data_min) / bin_width);
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

    free(slice_data);

    /* Create histogram data structure */
    HistogramData *hist_data = (HistogramData *)malloc(sizeof(HistogramData));
    hist_data->bin_counts = bin_counts;
    hist_data->bin_centers = bin_centers;
    hist_data->n_bins = n_bins;
    hist_data->count_max = count_max;
    hist_data->bin_min = data_min;
    hist_data->bin_max = data_max;
    hist_data->mean = mean;
    hist_data->std = std;
    hist_data->skewness = skewness;
    snprintf(hist_data->title, sizeof(hist_data->title), "%s Distribution at %s Layer %d",
             pf->variables[pf->current_var], axis_names[axis], slice_idx + 1);
    snprintf(hist_data->xlabel, sizeof(hist_data->xlabel), "%s", pf->variables[pf->current_var]);

    /* Create popup shell */
    Widget popup_shell = XtVaCreatePopupShell("Distribution",
        transientShellWidgetClass, toplevel,
        XtNwidth, 600,
        XtNheight, 400,
        NULL);

    Widget popup_form = XtVaCreateManagedWidget("form",
        formWidgetClass, popup_shell,
        NULL);

    /* Histogram canvas */
    Widget hist_canvas = XtVaCreateManagedWidget("histogram",
        simpleWidgetClass, popup_form,
        XtNwidth, 580,
        XtNheight, 320,
        XtNborderWidth, 1,
        NULL);

    /* Add expose event handler */
    XtAddEventHandler(hist_canvas, ExposureMask, False, histogram_expose_handler, hist_data);

    /* Close button */
    Widget close_button = XtVaCreateManagedWidget("Close",
        commandWidgetClass, popup_form,
        XtNfromVert, hist_canvas,
        NULL);

    XtAddCallback(close_button, XtNcallback, close_histogram_popup_callback, hist_data);

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
                for (int jj = 0; jj < height; jj++) {
                    for (int ii = 0; ii < width; ii++) {
                        int idx = jj * width + ii;
                        double z_coord = pf->prob_lo[2] + (jj + 0.5) * (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                        y_coord_slice[idx] = z_coord;
                    }
                }
            } else {
                /* X-slice: lat vs Z */
                read_variable_data(pf, lat_idx);
                extract_slice_from_data(pf->data, pf, x_coord_slice, pf->slice_axis, pf->slice_idx);
                for (int jj = 0; jj < height; jj++) {
                    for (int ii = 0; ii < width; ii++) {
                        int idx = jj * width + ii;
                        double z_coord = pf->prob_lo[2] + (jj + 0.5) * (pf->prob_hi[2] - pf->prob_lo[2]) / pf->grid_dims[2];
                        y_coord_slice[idx] = z_coord;
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
                /* Flip Y to match image rendering (higher j = higher physical Y = screen top) */
                int flipped_j = height - 1 - j;
                screen_x = render_offset_x + (int)((double)i * render_width / width);
                screen_y = render_offset_y + (int)((double)flipped_j * render_height / height);

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

void cleanup(PlotfileData *pf) {
    if (pf->data) free(pf->data);
    if (pixel_data) free(pixel_data);
    if (current_slice_data) free(current_slice_data);
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

    /* Log-X binning: use log-spaced bins */
    int use_log_bins = pd->log_x && rmin > 0;
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
    double *bin_sd_counts = (double *)calloc(n_bins, sizeof(double));
    double *bin_mass = (double *)calloc(n_bins, sizeof(double));

    if (use_log_bins) {
        for (int i = 0; i < n_bins; i++) {
            double log_center = log_rmin + (i + 0.5) * log_bin_width;
            bin_centers[i] = pow(10.0, log_center);
        }
    } else {
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

    /* Find max for scaling */
    double count_max = 0;
    for (int i = 0; i < n_bins; i++) {
        if (display_values[i] > count_max) count_max = display_values[i];
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

    hist->bin_counts = display_values;
    hist->bin_centers = bin_centers;
    hist->n_bins = n_bins;
    hist->count_max = count_max;
    hist->bin_min = rmin;
    hist->bin_max = rmax;
    hist->mean = mean;
    hist->std = std;
    hist->skewness = skewness;
    hist->kurtosis = kurtosis;
    snprintf(hist->title, sizeof(hist->title), "%s", sdm_metric_titles[pd->current_metric]);

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

    const char *ylabel = sdm_metric_ylabels[pd->current_metric];

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
    const char *basename = strrchr(plotfile_dir, '/');
    basename = basename ? basename + 1 : plotfile_dir;

    if (n_timesteps > 1) {
        snprintf(text, sizeof(text), "SDM: %s  |  Particles: %d  |  Metric: %s  |  Step %d/%d",
                 basename, pd->n_particles, sdm_metric_labels[pd->current_metric],
                 current_timestep + 1, n_timesteps);
    } else {
        snprintf(text, sizeof(text), "SDM: %s  |  Particles: %d  |  Metric: %s",
                 basename, pd->n_particles, sdm_metric_labels[pd->current_metric]);
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

/* SDM Settings apply callback */
void sdm_settings_apply_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (!global_pd) return;

    Arg args[1];
    String cutoff_str, binwidth_str;

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

/* SDM Settings button callback - open popup dialog */
void sdm_settings_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Arg args[10];
    int n;
    Widget form_w, label, button_w;
    char cutoff_str[64], binwidth_str[64];

    /* Format current values */
    if (global_pd && global_pd->cutoff_radius > 0)
        snprintf(cutoff_str, sizeof(cutoff_str), "%.4f", global_pd->cutoff_radius);
    else
        cutoff_str[0] = '\0';

    if (global_pd && global_pd->custom_bin_width > 0)
        snprintf(binwidth_str, sizeof(binwidth_str), "%.4f", global_pd->custom_bin_width);
    else
        binwidth_str[0] = '\0';

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

    /* Apply button */
    n = 0;
    XtSetArg(args[n], XtNfromVert, bw_label); n++;
    XtSetArg(args[n], XtNlabel, "Apply"); n++;
    button_w = XtCreateManagedWidget("apply", commandWidgetClass, form_w, args, n);
    XtAddCallback(button_w, XtNcallback, sdm_settings_apply_callback, NULL);

    /* Close button */
    n = 0;
    XtSetArg(args[n], XtNfromVert, bw_label); n++;
    XtSetArg(args[n], XtNfromHoriz, button_w); n++;
    XtSetArg(args[n], XtNlabel, "Close"); n++;
    button_w = XtCreateManagedWidget("close", commandWidgetClass, form_w, args, n);
    XtAddCallback(button_w, XtNcallback, sdm_settings_close_callback, NULL);

    XtRealizeWidget(sdm_dialog_shell);
    XtPopup(sdm_dialog_shell, XtGrabExclusive);

    /* Set keyboard focus to cutoff text input */
    XtSetKeyboardFocus(sdm_dialog_shell, sdm_settings_text_cutoff);
    XSync(display, False);
    Time time_val = CurrentTime;
    XtCallAcceptFocus(sdm_settings_text_cutoff, &time_val);

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

int main(int argc, char **argv) {
    PlotfileData pf = {0};
    Arg args[2];
    char check_path[MAX_PATH];
    const char *prefix = "plt";  /* Default prefix */

    /* Check for --sdm flag */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--sdm") == 0) {
            sdm_mode = 1;
            /* Shift remaining args over this flag */
            for (int j = i; j < argc - 1; j++) {
                argv[j] = argv[j + 1];
            }
            argc--;
            i--;  /* Re-check this position */
        }
    }

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [--sdm] <plotfile_directory> [prefix]\n", argv[0]);
        fprintf(stderr, "  Single plotfile:    %s plt00100\n", argv[0]);
        fprintf(stderr, "  Multi-timestep:     %s /path/to/dir plt\n", argv[0]);
        fprintf(stderr, "  With prefix plt2d:  %s /path/to/dir plt2d\n", argv[0]);
        fprintf(stderr, "  SDM mode:           %s --sdm plt00100\n", argv[0]);
        fprintf(stderr, "  SDM multi-timestep: %s --sdm /path/to/dir plt\n", argv[0]);
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
    } else {
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
            free(sdm_hist_data);
        }
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
            }

            if (changed) {
                update_layer_label(global_pf);
                update_info_label(global_pf);
                render_slice(global_pf);
            }
        }
        
        XtDispatchEvent(&event);
    }
    
    cleanup(&pf);
    return 0;
}
