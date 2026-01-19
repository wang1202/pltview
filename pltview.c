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

#define MAX_VARS 64
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
} PlotfileData;

/* Colormap */
typedef struct {
    unsigned char r, g, b;
} RGB;

/* X11 globals */
Display *display;
Widget toplevel, form, canvas_widget, var_box, info_label, slice_scroll;
Widget axis_box;
Window canvas;
GC gc, text_gc;
XImage *ximage;
int screen;
unsigned long *pixel_data;
int canvas_width = 800;
int canvas_height = 600;
Pixmap pixmap;
XFontStruct *font;

/* Function prototypes */
int read_header(PlotfileData *pf);
int read_cell_h(PlotfileData *pf);
int read_variable_data(PlotfileData *pf, int var_idx);
void extract_slice(PlotfileData *pf, double *slice, int axis, int idx);
void apply_colormap(double *data, int width, int height, 
                   unsigned long *pixels, double vmin, double vmax);
RGB viridis_colormap(double t);
void init_gui(PlotfileData *pf, int argc, char **argv);
void render_slice(PlotfileData *pf);
void update_info_label(PlotfileData *pf);
void var_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void axis_button_callback(Widget w, XtPointer client_data, XtPointer call_data);
void scroll_callback(Widget w, XtPointer client_data, XtPointer call_data);
void canvas_expose_callback(Widget w, XtPointer client_data, XtPointer call_data);
void cleanup(PlotfileData *pf);

/* Global pointer for callbacks */
PlotfileData *global_pf = NULL;

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
    
    /* Skip to domain box (lines: num_levels, low, high, refinement) */
    fgets(line, MAX_LINE, fp);  /* num levels */
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
    
    printf("Loaded: %s\n", pf->plotfile_dir);
    printf("Variables: %d (", pf->n_vars);
    for (i = 0; i < pf->n_vars && i < 5; i++) {
        printf("%s%s", pf->variables[i], i < pf->n_vars-1 ? ", " : "");
    }
    if (pf->n_vars > 5) printf("...");
    printf(")\n");
    printf("Grid: %d x %d x %d\n", pf->grid_dims[0], pf->grid_dims[1], pf->grid_dims[2]);
    printf("Time: %.3f\n", pf->time);
    
    return 0;
}

/* Read Cell_H to get box layout and FabOnDisk mapping */
int read_cell_h(PlotfileData *pf) {
    char path[MAX_PATH];
    char line[MAX_LINE];
    FILE *fp;
    int i;
    
    snprintf(path, MAX_PATH, "%s/Level_0/Cell_H", pf->plotfile_dir);
    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open %s\n", path);
        return -1;
    }
    
    /* Skip first few lines until we find box definitions */
    int box_count = 0;
    while (fgets(line, MAX_LINE, fp)) {
        if (strncmp(line, "((", 2) == 0) {
            /* Parse box: ((lo_x,lo_y,lo_z) (hi_x,hi_y,hi_z) ...) */
            char *p = line + 2;
            for (i = 0; i < pf->ndim; i++) {
                while (*p && !isdigit(*p) && *p != '-') p++;
                pf->boxes[box_count].lo[i] = atoi(p);
                while (*p && (isdigit(*p) || *p == '-')) p++;
            }
            for (i = 0; i < pf->ndim; i++) {
                while (*p && !isdigit(*p) && *p != '-') p++;
                pf->boxes[box_count].hi[i] = atoi(p);
                while (*p && (isdigit(*p) || *p == '-')) p++;
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
    printf("Found %d boxes\n", pf->n_boxes);
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
        
        snprintf(path, MAX_PATH, "%s/Level_0/%s", pf->plotfile_dir, box->filename);
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

/* Apply colormap to data */
void apply_colormap(double *data, int width, int height, 
                   unsigned long *pixels, double vmin, double vmax) {
    int i, j;
    double range = vmax - vmin;
    if (range < 1e-10) range = 1.0;
    
    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++) {
            double val = data[j * width + i];
            double t = (val - vmin) / range;
            RGB color = viridis_colormap(t);
            pixels[j * width + i] = (color.r << 16) | (color.g << 8) | color.b;
        }
    }
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
    
    /* Add variable buttons (limit to first 10) */
    int max_vars = pf->n_vars < 10 ? pf->n_vars : 10;
    for (i = 0; i < max_vars; i++) {
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
    
    /* Slice scrollbar */
    n = 0;
    XtSetArg(args[n], XtNfromVert, canvas_widget); n++;
    XtSetArg(args[n], XtNfromHoriz, axis_box); n++;
    XtSetArg(args[n], XtNwidth, 400); n++;
    XtSetArg(args[n], XtNheight, 20); n++;
    XtSetArg(args[n], XtNorientation, XtorientHorizontal); n++;
    XtSetArg(args[n], XtNbottom, XawChainBottom); n++;
    XtSetArg(args[n], XtNleft, XawChainLeft); n++;
    XtSetArg(args[n], XtNright, XawChainRight); n++;
    slice_scroll = XtCreateManagedWidget("sliceScroll", scrollbarWidgetClass, form, args, n);
    XtAddCallback(slice_scroll, XtNscrollProc, scroll_callback, NULL);
    XtAddCallback(slice_scroll, XtNjumpProc, scroll_callback, NULL);
    
    XtRealizeWidget(toplevel);
    
    /* Get canvas window and create GC */
    canvas = XtWindow(canvas_widget);
    gc = XCreateGC(display, canvas, 0, NULL);
    
    /* Create text GC for overlay */
    text_gc = XCreateGC(display, canvas, 0, NULL);
    XSetForeground(display, text_gc, WhitePixel(display, screen));
    XSetBackground(display, text_gc, BlackPixel(display, screen));
    if (font) XSetFont(display, text_gc, font->fid);
    
    /* Allocate pixel buffer */
    pixel_data = (unsigned long *)malloc(canvas_width * canvas_height * sizeof(unsigned long));
    pixmap = XCreatePixmap(display, canvas, canvas_width, canvas_height, 
                          DefaultDepth(display, screen));
}

/* Update info label */
void update_info_label(PlotfileData *pf) {
    char text[256];
    const char *axis_names[] = {"X", "Y", "Z"};
    int max_idx = pf->grid_dims[pf->slice_axis] - 1;
    
    snprintf(text, sizeof(text), 
             "%s | Axis: %s | Slice: %d/%d | Time: %.3f",
             pf->variables[pf->current_var],
             axis_names[pf->slice_axis],
             pf->slice_idx, max_idx,
             pf->time);
    
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
        global_pf->slice_idx = global_pf->grid_dims[axis] / 2;
        
        /* Update scrollbar */
        float max_val = (float)global_pf->grid_dims[axis];
        float shown = 1.0 / max_val;
        Arg args[2];
        XtSetArg(args[0], XtNshown, shown);
        XtSetArg(args[1], XtNtopOfThumb, (float)global_pf->slice_idx / max_val);
        XtSetValues(slice_scroll, args, 2);
        
        update_info_label(global_pf);
        render_slice(global_pf);
    }
}

/* Scrollbar callback */
void scroll_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    float percent = *(float *)call_data;
    if (global_pf) {
        int max_idx = global_pf->grid_dims[global_pf->slice_axis] - 1;
        global_pf->slice_idx = (int)(percent * max_idx);
        if (global_pf->slice_idx > max_idx) global_pf->slice_idx = max_idx;
        if (global_pf->slice_idx < 0) global_pf->slice_idx = 0;
        
        update_info_label(global_pf);
        render_slice(global_pf);
    }
}

/* Canvas expose callback */
void canvas_expose_callback(Widget w, XtPointer client_data, XtPointer call_data) {
    if (global_pf && global_pf->data) {
        render_slice(global_pf);
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
    
    /* Find min/max */
    for (i = 0; i < width * height; i++) {
        if (slice[i] < vmin) vmin = slice[i];
        if (slice[i] > vmax) vmax = slice[i];
    }
    
    /* Apply colormap */
    apply_colormap(slice, width, height, pixel_data, vmin, vmax);
    
    /* Create XImage and draw to pixmap */
    Visual *visual = DefaultVisual(display, screen);
    ximage = XCreateImage(display, visual, 24, ZPixmap, 0,
                         (char *)pixel_data, width, height, 32, 0);
    
    XPutImage(display, pixmap, gc, ximage, 0, 0, 0, 0, width, height);
    
    /* Draw pixmap to canvas */
    XCopyArea(display, pixmap, canvas, gc, 0, 0, width, height, 0, 0);
    
    /* Draw text overlay */
    snprintf(stats_text, sizeof(stats_text), "min: %.3e  max: %.3e", vmin, vmax);
    XDrawImageString(display, canvas, text_gc, 10, height - 10, 
                    stats_text, strlen(stats_text));
    
    XFlush(display);
    
    printf("Rendered: %s, slice %d/%d (%.3e to %.3e)\n", 
           pf->variables[pf->current_var], pf->slice_idx,
           pf->grid_dims[pf->slice_axis]-1, vmin, vmax);
    
    free(slice);
    ximage->data = NULL;
    XDestroyImage(ximage);
}

void cleanup(PlotfileData *pf) {
    if (pf->data) free(pf->data);
    if (pixel_data) free(pixel_data);
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
    if (read_cell_h(&pf) < 0) return 1;
    
    /* Load first variable */
    pf.current_var = 0;
    pf.slice_axis = 2;  /* Z */
    pf.slice_idx = pf.grid_dims[2] / 2;
    
    read_variable_data(&pf, 0);
    
    /* Initialize GUI */
    init_gui(&pf, argc, argv);
    
    /* Set initial scrollbar position */
    float max_val = (float)pf.grid_dims[pf.slice_axis];
    float shown = 1.0 / max_val;
    XtSetArg(args[0], XtNshown, shown);
    XtSetArg(args[1], XtNtopOfThumb, (float)pf.slice_idx / max_val);
    XtSetValues(slice_scroll, args, 2);
    
    update_info_label(&pf);
    render_slice(&pf);
    
    printf("\nGUI Controls:\n");
    printf("  Click variable buttons to change variable\n");
    printf("  Click X/Y/Z buttons to switch axis\n");
    printf("  Drag scrollbar to navigate slices\n\n");
    
    XtAppMainLoop(XtWidgetToApplicationContext(toplevel));
    
    cleanup(&pf);
    return 0;
}
