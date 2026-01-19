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

#define MAX_VARS 128
#define MAX_BOXES 1024
#define MAX_PATH 512
#define MAX_LINE 1024

/* Data structures */
typedef struct {
    int lo[3];
    int hi[3];
    char filename[64];
} Box;

typedef struct {
    char plotfile_dir[MAX_PATH];
    char variables[MAX_VARS][64];
    int n_vars;
    int ndim;
    double time;
    int grid_dims[3];
    Box boxes[MAX_BOXES];
    int n_boxes;
    double *data;  /* Current variable data */
    int current_var;
    int slice_axis;
    int slice_idx;
    int colormap;  /* 0=viridis, 1=jet, 2=turbo, 3=plasma */
    int current_level;  /* Current AMR level */
    int n_levels;       /* Number of AMR levels */
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
    char xlabel[64];    /* X-axis label */
} PlotData;

/* Popup window data */
typedef struct {
    Widget shell;
    PlotData *plot_data_array[3];
} PopupData;

/* X11 globals */
Display *display;
Widget toplevel, form, canvas_widget, var_box, info_label;
Widget axis_box, cmap_box, nav_box, colorbar_widget, layer_label;
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

/* Function prototypes */
int detect_levels(PlotfileData *pf);
int read_header(PlotfileData *pf);
int read_cell_h(PlotfileData *pf);
int read_variable_data(PlotfileData *pf, int var_idx);
void extract_slice(PlotfileData *pf, double *slice, int axis, int idx);
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
void draw_colorbar(double vmin, double vmax, int cmap_type);
void cmap_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void colorbar_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void init_gui(PlotfileData *pf, int argc, char **argv);
void render_slice(PlotfileData *pf);
void update_info_label(PlotfileData *pf);
void var_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void axis_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void level_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void nav_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void jump_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void update_layer_label(PlotfileData *pf);
void canvas_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void canvas_motion_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void canvas_button_handler(Widget w, XtPointer client_data, XEvent *event, Boolean *continue_dispatch);
void show_line_profiles(PlotfileData *pf, int data_x, int data_y);
void cleanup(PlotfileData *pf);

/* Global pointer for callbacks */
PlotfileData *global_pf = NULL;

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
    
    /* Skip to domain box (lines: low, high, refinement) */
    fgets(line, MAX_LINE, fp);  /* low */
    fgets(line, MAX_LINE, fp);  /* high */
    fgets(line, MAX_LINE, fp);  /* refinement */
    
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
    
    /* Update grid dimensions based on level domain */
    for (i = 0; i < pf->ndim; i++) {
        pf->grid_dims[i] = level_hi[i] - level_lo[i] + 1;
    }
    
    printf("Level %d: Found %d boxes, Grid: %d x %d x %d\n", 
           pf->current_level, pf->n_boxes, 
           pf->grid_dims[0], pf->grid_dims[1], pf->grid_dims[2]);
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
        size_t idx = 0;
        for (k = 0; k < box_dims[2]; k++) {
            for (j = 0; j < box_dims[1]; j++) {
                for (i = 0; i < box_dims[0]; i++) {
                    int gx = box->lo[0] + i;
                    int gy = box->lo[1] + j;
                    int gz = box->lo[2] + k;
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

/* Draw colorbar */
void draw_colorbar(double vmin, double vmax, int cmap_type) {
    int height = 256, width = 30;
    
    /* Clear colorbar with white background */
    XSetForeground(display, colorbar_gc, WhitePixel(display, screen));
    XFillRectangle(display, colorbar, colorbar_gc, 0, 0, 100, canvas_height);
    
    /* Draw colorbar as solid rectangles */
    for (int i = 0; i < height; i++) {
        double t = (double)(height - 1 - i) / (height - 1);
        RGB color = get_colormap_rgb(t, cmap_type);
        unsigned long pixel = (color.r << 16) | (color.g << 8) | color.b;
        
        XSetForeground(display, colorbar_gc, pixel);
        int y = (i * canvas_height) / height;
        int h = ((i + 1) * canvas_height) / height - y;
        if (h < 1) h = 1;
        XFillRectangle(display, colorbar, colorbar_gc, 0, y, width, h);
    }
    
    /* Draw labels */
    char text[32];
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    snprintf(text, sizeof(text), "%.2e", vmax);
    XDrawString(display, colorbar, text_gc, width + 5, 15, text, strlen(text));
    snprintf(text, sizeof(text), "%.2e", vmin);
    XDrawString(display, colorbar, text_gc, width + 5, canvas_height - 5, text, strlen(text));
    
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
    
    /* COLUMN 1: Axis and Navigation buttons */
    /* Axis buttons box */
    n = 0;
    XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
    XtSetArg(args[n], XtNfromHoriz, var_box); n++;
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
    
    /* Navigation buttons (+/-) in Column 1 */
    n = 0;
    XtSetArg(args[n], XtNfromVert, axis_box); n++;
    XtSetArg(args[n], XtNfromHoriz, var_box); n++;
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    nav_box = XtCreateManagedWidget("navBox", boxWidgetClass, form, args, n);
    
    const char *nav_labels[] = {"-", "+"};
    for (i = 0; i < 2; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, nav_labels[i]); n++;
        button = XtCreateManagedWidget(nav_labels[i], commandWidgetClass, nav_box, args, n);
        XtAddCallback(button, XtNcallback, nav_button_callback, (XtPointer)(long)i);
    }
    
    /* Layer display label */
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
    
    /* COLUMN 2: Level and Colormap buttons */
    /* Level buttons box (only if multiple levels exist) */
    Widget level_box = NULL;
    if (pf->n_levels > 1) {
        n = 0;
        XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
        XtSetArg(args[n], XtNfromHoriz, nav_box); n++;
        XtSetArg(args[n], XtNborderWidth, 1); n++;
        XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
        XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
        XtSetArg(args[n], XtNleft, XawChainLeft); n++;
        level_box = XtCreateManagedWidget("levelBox", boxWidgetClass, form, args, n);
        
        /* Add level buttons (limit to 10 levels) */
        int max_levels = pf->n_levels < 10 ? pf->n_levels : 10;
        for (i = 0; i < max_levels; i++) {
            n = 0;
            snprintf(label_text, sizeof(label_text), "Level %d", i);
            XtSetArg(args[n], XtNlabel, label_text); n++;
            button = XtCreateManagedWidget(label_text, commandWidgetClass, level_box, args, n);
            XtAddCallback(button, XtNcallback, level_button_callback, (XtPointer)(long)i);
        }
    }
    
    /* Colormap buttons in Column 2 */
    n = 0;
    if (pf->n_levels > 1) {
        XtSetArg(args[n], XtNfromVert, level_box); n++;
        XtSetArg(args[n], XtNfromHoriz, nav_box); n++;
    } else {
        XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
        XtSetArg(args[n], XtNfromHoriz, nav_box); n++;
    }
    XtSetArg(args[n], XtNborderWidth, 1); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    cmap_box = XtCreateManagedWidget("cmapBox", boxWidgetClass, form, args, n);
    
    const char *cmap_labels[] = {"viridis", "jet", "turbo", "plasma", "hot", "cool", "gray", "magma"};
    for (i = 0; i < 8; i++) {
        n = 0;
        XtSetArg(args[n], XtNlabel, cmap_labels[i]); n++;
        button = XtCreateManagedWidget(cmap_labels[i], commandWidgetClass, cmap_box, args, n);
        XtAddCallback(button, XtNcallback, cmap_button_callback, (XtPointer)(long)i);
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
        
        update_layer_label(global_pf);
        update_info_label(global_pf);
        render_slice(global_pf);
    }
}

/* Level button callback */
void level_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int level = (int)(long)client_data;
    if (global_pf && level < global_pf->n_levels) {
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

/* Dialog OK callback for layer jump */
void jump_dialog_ok_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Widget dialog = (Widget)client_data;
    Widget shell = XtParent(dialog);
    
    if (global_pf) {
        char *value = XawDialogGetValueString(dialog);
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
    
    XtPopdown(shell);
    XtDestroyWidget(shell);
}

/* Dialog Cancel callback */
void jump_dialog_cancel_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    Widget dialog = (Widget)client_data;
    Widget shell = XtParent(dialog);
    XtPopdown(shell);
    XtDestroyWidget(shell);
}

/* Jump button callback - open dialog to jump to layer */
void jump_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf) {
        Arg args[10];
        int n;
        Widget dialog, dialog_shell, text_widget;
        char prompt[128];
        int max_idx = global_pf->grid_dims[global_pf->slice_axis];
        
        snprintf(prompt, sizeof(prompt), "Jump to layer (1-%d):", max_idx);
        
        n = 0;
        XtSetArg(args[n], XtNtitle, "Jump to Layer"); n++;
        dialog_shell = XtCreatePopupShell("jumpDialog", transientShellWidgetClass, toplevel, args, n);
        
        n = 0;
        XtSetArg(args[n], XtNlabel, prompt); n++;
        XtSetArg(args[n], XtNvalue, ""); n++;
        dialog = XtCreateManagedWidget("dialog", dialogWidgetClass, dialog_shell, args, n);
        
        XawDialogAddButton(dialog, "OK", jump_dialog_ok_callback, (XtPointer)dialog);
        XawDialogAddButton(dialog, "Cancel", jump_dialog_cancel_callback, (XtPointer)dialog);
        
        XtRealizeWidget(dialog_shell);
        
        /* Get the text input widget and set focus to it for X11 forwarding */
        text_widget = XtNameToWidget(dialog, "value");
        if (text_widget) {
            XtSetKeyboardFocus(dialog_shell, text_widget);
        }
        
        XtPopup(dialog_shell, XtGrabNone);
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

/* Colormap button callback */
void cmap_button_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    int cmap = (int)(long)client_data;
    if (global_pf) {
        global_pf->colormap = cmap;
        render_slice(global_pf);
        draw_colorbar(current_vmin, current_vmax, cmap);
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
        draw_colorbar(current_vmin, current_vmax, global_pf->colormap);
    }
}
void render_slice(PlotfileData *pf) {
    int width, height;
    double *slice;
    double vmin = 1e30, vmax = -1e30;
    int i;
    char stats_text[128];
    
    /* Determine slice dimensions */
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
    
    slice = (double *)malloc(width * height * sizeof(double));
    extract_slice(pf, slice, pf->slice_axis, pf->slice_idx);
    
    /* Store current slice for mouse interaction */
    if (current_slice_data) free(current_slice_data);
    current_slice_data = (double *)malloc(width * height * sizeof(double));
    memcpy(current_slice_data, slice, width * height * sizeof(double));
    slice_width = width;
    slice_height = height;
    
    /* Find min/max */
    for (i = 0; i < width * height; i++) {
        if (slice[i] < vmin) vmin = slice[i];
        if (slice[i] > vmax) vmax = slice[i];
    }
    
    /* Store current vmin/vmax for colorbar */
    current_vmin = vmin;
    current_vmax = vmax;
    
    /* Apply colormap */
    apply_colormap(slice, width, height, pixel_data, vmin, vmax, pf->colormap);
    
    /* Clear canvas with white background */
    XSetForeground(display, gc, WhitePixel(display, screen));
    XFillRectangle(display, canvas, gc, 0, 0, canvas_width, canvas_height);
    
    /* Calculate scaling to maintain aspect ratio */
    double data_aspect = (double)width / height;
    double canvas_aspect = (double)canvas_width / canvas_height;
    
    int offset_x, offset_y;
    if (data_aspect > canvas_aspect) {
        /* Width-limited */
        render_width = canvas_width;
        render_height = (int)(canvas_width / data_aspect);
        offset_x = 0;
        offset_y = (canvas_height - render_height) / 2;
    } else {
        /* Height-limited */
        render_width = (int)(canvas_height * data_aspect);
        render_height = canvas_height;
        offset_x = (canvas_width - render_width) / 2;
        offset_y = 0;
    }
    
    /* Store rendering parameters for mouse interaction */
    render_offset_x = offset_x;
    render_offset_y = offset_y;
    
    /* Draw pixels as filled rectangles with correct aspect ratio */
    double pixel_width = (double)render_width / width;
    double pixel_height = (double)render_height / height;
    
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
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
    
    /* Draw text overlay */
    snprintf(stats_text, sizeof(stats_text), "min: %.3e  max: %.3e", vmin, vmax);
    XSetForeground(display, text_gc, BlackPixel(display, screen));
    XSetBackground(display, text_gc, WhitePixel(display, screen));
    XDrawImageString(display, canvas, text_gc, 10, canvas_height - 10, 
                    stats_text, strlen(stats_text));
    
    /* Draw colorbar */
    draw_colorbar(vmin, vmax, pf->colormap);
    
    XFlush(display);
    
    printf("Rendered: %s, slice %d/%d (%.3e to %.3e)\n", 
           pf->variables[pf->current_var], pf->slice_idx + 1,
           pf->grid_dims[pf->slice_axis], vmin, vmax);
    
    free(slice);
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

void cleanup(PlotfileData *pf) {
    if (pf->data) free(pf->data);
    if (pixel_data) free(pixel_data);
    if (current_slice_data) free(current_slice_data);
}

int main(int argc, char **argv) {
    PlotfileData pf = {0};
    Arg args[2];
    
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <plotfile_directory>\n", argv[0]);
        return 1;
    }
    
    strncpy(pf.plotfile_dir, argv[1], MAX_PATH - 1);
    
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
    update_info_label(&pf);
    render_slice(&pf);
    
    printf("\nGUI Controls:\n");
    printf("  Click variable buttons to change variable\n");
    printf("  Click X/Y/Z buttons to switch axis\n");
    if (pf.n_levels > 1) {
        printf("  Click Level 0/Level 1/... buttons to switch level\n");
    }
    printf("  Click colormap buttons (viridis/jet/turbo/plasma/hot/cool/gray/magma)\n");
    printf("  Click +/- buttons to navigate layers (or use keyboard arrows)\n\n");
    
    /* Main event loop with expose and keyboard handling */
    XtAppContext app_context = XtWidgetToApplicationContext(toplevel);
    while (1) {
        XEvent event;
        XtAppNextEvent(app_context, &event);
        
        /* Handle expose events */
        if (event.type == Expose) {
            if (event.xexpose.window == canvas && global_pf && global_pf->data) {
                render_slice(global_pf);
            } else if (event.xexpose.window == colorbar && global_pf) {
                draw_colorbar(current_vmin, current_vmax, global_pf->colormap);
            }
        }
        /* Handle keyboard events */
        else if (event.type == KeyPress && global_pf) {
            KeySym key = XLookupKeysym(&event.xkey, 0);
            int max_idx = global_pf->grid_dims[global_pf->slice_axis] - 1;
            int changed = 0;
            
            if (key == XK_plus || key == XK_equal || key == XK_Right) {
                if (global_pf->slice_idx < max_idx) {
                    global_pf->slice_idx++;
                    changed = 1;
                }
            } else if (key == XK_minus || key == XK_underscore || key == XK_Left) {
                if (global_pf->slice_idx > 0) {
                    global_pf->slice_idx--;
                    changed = 1;
                }
            } else if (key >= XK_0 && key <= XK_3) {
                /* Switch colormap with 0-3 keys */
                global_pf->colormap = key - XK_0;
                changed = 1;
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
