#!/usr/bin/env python3
"""
PLTView - A simple viewer for AMReX plotfiles, inspired by ncview
"""

import os
# Set DISPLAY environment variable for XQuartz
if 'DISPLAY' not in os.environ:
    os.environ['DISPLAY'] = ':0'

import numpy as np
import matplotlib
# Use macOS backend which works with XQuartz when DISPLAY is set
matplotlib.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, Button
import struct


class AMReXPlotfile:
    """Class to read AMReX plotfile format"""
    
    def __init__(self, plotfile_dir):
        self.plotfile_dir = plotfile_dir
        self.variables = []
        self.grid_dims = None
        self.num_levels = None
        self.time = None
        self.data_cache = {}
        
        self.read_header()
    
    def read_header(self):
        """Parse the Header file"""
        header_path = os.path.join(self.plotfile_dir, 'Header')
        
        with open(header_path, 'r') as f:
            lines = f.readlines()
        
        # Parse header
        idx = 0
        version = lines[idx].strip()
        idx += 1
        
        # Number of variables
        n_vars = int(lines[idx].strip())
        idx += 1
        
        # Variable names
        self.variables = []
        for i in range(n_vars):
            self.variables.append(lines[idx].strip())
            idx += 1
        
        # Dimensionality
        self.ndim = int(lines[idx].strip())
        idx += 1
        
        # Time
        self.time = float(lines[idx].strip())
        idx += 1
        
        # Number of levels
        self.num_levels = int(lines[idx].strip())
        idx += 1
        
        # Skip to grid dimensions
        # Low corner
        idx += 1
        # High corner
        idx += 1
        
        # Refinement ratios
        idx += 1
        
        # Domain box
        domain_line = lines[idx].strip()
        idx += 1
        # Parse ((lo_x,lo_y,lo_z) (hi_x,hi_y,hi_z) (type_x,type_y,type_z))
        domain_line = domain_line.replace('(', '').replace(')', '').replace(',', ' ')
        parts = domain_line.split()
        lo = [int(parts[i]) for i in range(self.ndim)]
        hi = [int(parts[i]) for i in range(self.ndim, 2*self.ndim)]
        
        self.grid_dims = tuple([hi[i] - lo[i] + 1 for i in range(self.ndim)])
        
        # Number of steps
        self.timestep = int(lines[idx].strip())
        idx += 1
        
        print(f"Loaded plotfile: {self.plotfile_dir}")
        print(f"Variables: {', '.join(self.variables)}")
        print(f"Grid dimensions: {self.grid_dims}")
        print(f"Time: {self.time}")
    
    def read_variable_data(self, var_name, level=0):
        """Read data for a specific variable"""
        if var_name in self.data_cache:
            return self.data_cache[var_name]
        
        var_idx = self.variables.index(var_name)
        
        # Initialize data array
        data = np.zeros(self.grid_dims[::-1])  # AMReX uses Fortran ordering
        
        # Read from all Cell_D files
        level_dir = os.path.join(self.plotfile_dir, f'Level_{level}')
        
        # Read Cell_H header to get box layout
        cell_h_path = os.path.join(level_dir, 'Cell_H')
        with open(cell_h_path, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        
        # Parse header
        # Line 0: version
        # Line 1: how
        # Line 2: n_components
        n_components = int(lines[2])
        # Line 3: ngrow
        # Line 4: number of boxes info
        
        # Find the box definitions (between parentheses section)
        box_start_idx = None
        for i, line in enumerate(lines):
            if line.startswith('('):
                box_start_idx = i
                break
        
        if box_start_idx is None:
            raise ValueError("Could not find box definitions in Cell_H")
        
        # Count boxes - they start after the first '(' line
        boxes = []
        idx = box_start_idx + 1
        while idx < len(lines) and lines[idx].startswith('(('):
            box_line = lines[idx]
            # Parse ((lo_x,lo_y,lo_z) (hi_x,hi_y,hi_z) (type_x,type_y,type_z))
            box_line = box_line.replace('(', '').replace(')', '').replace(',', ' ')
            parts = box_line.split()
            
            lo = [int(parts[i]) for i in range(self.ndim)]
            hi = [int(parts[i]) for i in range(self.ndim, 2*self.ndim)]
            boxes.append((lo, hi))
            idx += 1
        
        n_boxes = len(boxes)
        
        # After the box definitions, there should be FabOnDisk mappings
        # Parse the FabOnDisk section to get the actual file names
        fab_files = []
        while idx < len(lines):
            if lines[idx].startswith('FabOnDisk:'):
                parts = lines[idx].split()
                fab_file = parts[1]  # e.g., "Cell_D_00000"
                fab_files.append(fab_file)
                idx += 1
            else:
                idx += 1
        
        if len(fab_files) != n_boxes:
            print(f"Warning: Found {len(fab_files)} FabOnDisk entries but {n_boxes} boxes")
            # Fallback to sequential numbering if no FabOnDisk section
            fab_files = [f'Cell_D_{i:05d}' for i in range(n_boxes)]
        
        # Read each box's data
        for box_idx in range(n_boxes):
            cell_d_path = os.path.join(level_dir, fab_files[box_idx])
            
            lo, hi = boxes[box_idx]
            box_shape = tuple([hi[i] - lo[i] + 1 for i in range(self.ndim)])
            box_size = np.prod(box_shape)
            
            # Read binary data from FAB file
            with open(cell_d_path, 'rb') as f:
                # Read FAB ASCII header
                header_line = b''
                while True:
                    char = f.read(1)
                    if char == b'\n':
                        break
                    header_line += char
                
                # Now we're at the binary data
                # Skip to the variable we want
                offset = var_idx * box_size * 8  # 8 bytes per double
                f.seek(f.tell() + offset)
                
                # Read data
                box_data = np.fromfile(f, dtype=np.float64, count=box_size)
                
                if len(box_data) != box_size:
                    print(f"Warning: Expected {box_size} values but got {len(box_data)}")
                
                # AMReX stores data in Fortran order (column-major): X varies fastest
                # For 3D: flat data is [x0y0z0, x1y0z0, ..., xN-1,y0,z0, x0y1z0, ..., x0y0z1, ...]
                # Reshape with order='F' to (nx, ny, nz), then transpose to get (nz, ny, nx)
                if self.ndim == 3:
                    box_data = box_data.reshape(box_shape, order='F')
                    # Swap axes: (nx, ny, nz) -> (nz, ny, nx)
                    box_data = np.moveaxis(box_data, [0, 1, 2], [2, 1, 0])
                elif self.ndim == 2:
                    box_data = box_data.reshape(box_shape, order='F')
                    box_data = np.moveaxis(box_data, [0, 1], [1, 0])
                
                # Debug: print box info and verify placement
                if box_idx in [0, 1, 37, 48, 63]:
                    print(f"  Box {box_idx}: lo={lo}, hi={hi}, box_shape={box_data.shape}")
                    print(f"    Slice: data[{lo[2]}:{hi[2]+1}, {lo[1]}:{hi[1]+1}, {lo[0]}:{hi[0]+1}]")
                    print(f"    Data range: {box_data.min():.6f} to {box_data.max():.6f}")
                    # Test: check what's in the target slice before assignment
                    if self.ndim == 3:
                        target_slice = data[lo[2]:hi[2]+1, lo[1]:hi[1]+1, lo[0]:hi[0]+1]
                        print(f"    Target slice shape: {target_slice.shape}, all zeros: {np.all(target_slice == 0)}")
                
                # Insert into global array
                if self.ndim == 3:
                    data[lo[2]:hi[2]+1, lo[1]:hi[1]+1, lo[0]:hi[0]+1] = box_data
                elif self.ndim == 2:
                    data[lo[1]:hi[1]+1, lo[0]:hi[0]+1] = box_data
        
        self.data_cache[var_name] = data
        return data


class PLTViewApp:
    """Main application window using matplotlib widgets"""
    
    def __init__(self, plotfile_dir):
        self.plotfile = AMReXPlotfile(plotfile_dir)
        self.current_var_idx = 0
        self.current_var = self.plotfile.variables[0]
        self.current_slice_idx = 0
        self.slice_axis = 2  # Default to Z-axis (0=X, 1=Y, 2=Z)
        self.current_cmap = 'viridis'  # Default colormap
        self.updating = False  # Flag to prevent re-entry
        self.current_data_3d = None  # Store current 3D data
        self.current_slice_data = None  # Store current 2D slice
        self.xlabel = None
        self.ylabel = None
        self.annotation = None  # For hover tooltip
        self.line_plot_window = None  # For line plot window
        
        self.setup_ui()
        self.update_plot()
    
    def setup_ui(self):
        """Setup the matplotlib-based user interface"""
        self.fig = plt.figure(figsize=(14, 8))
        
        # Main plot area
        self.ax_main = plt.subplot2grid((6, 5), (0, 1), colspan=3, rowspan=5)
        
        # Variable buttons area
        self.ax_vars = plt.subplot2grid((6, 5), (0, 0), rowspan=3)
        self.ax_vars.set_title('Variables', fontsize=10, loc='left')
        
        # Create variable selection buttons
        n_vars = len(self.plotfile.variables)
        var_labels = self.plotfile.variables
        
        # Limit visible variables if too many
        if n_vars > 15:
            var_labels = var_labels[:15] + [f'... ({n_vars - 15} more)']
        
        self.radio_vars = RadioButtons(self.ax_vars, var_labels[:15], active=0)
        self.radio_vars.on_clicked(self.on_variable_selected)
        
        # Colormap buttons area
        self.ax_cmap = plt.subplot2grid((6, 5), (3, 0), rowspan=2)
        self.ax_cmap.set_title('Colormap', fontsize=10, loc='left')
        cmap_options = ['viridis', 'turbo', 'jet', 'plasma', 'inferno', 'magma', 'seismic', 'RdBu_r']
        self.radio_cmap = RadioButtons(self.ax_cmap, cmap_options, active=0)
        self.radio_cmap.on_clicked(self.on_colormap_changed)
        
        # Slice axis buttons
        self.ax_axis = plt.subplot2grid((6, 5), (5, 0))
        self.radio_axis = RadioButtons(self.ax_axis, ['X', 'Y', 'Z'], active=2)
        self.radio_axis.on_clicked(self.on_axis_changed)
        
        # Slice slider
        self.ax_slider = plt.subplot2grid((6, 5), (5, 1), colspan=4)
        max_idx = self.plotfile.grid_dims[self.slice_axis] - 1
        self.slider = Slider(self.ax_slider, 'Slice', 0, max_idx, 
                            valinit=0, valstep=1)
        self.slider.on_changed(self.on_slice_changed)
        
        # Info text
        info_text = f"Time: {self.plotfile.time:.3f}\n"
        info_text += f"Grid: {' × '.join(map(str, self.plotfile.grid_dims))}"
        self.fig.text(0.02, 0.02, info_text, fontsize=9, 
                     verticalalignment='bottom', family='monospace')
        
        plt.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.12, 
                           hspace=0.3, wspace=0.3)
        
        # Connect mouse events
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.fig.canvas.mpl_connect('button_press_event', self.on_mouse_click)
    
    def on_variable_selected(self, label):
        """Handle variable selection"""
        try:
            idx = self.plotfile.variables.index(label)
            self.current_var_idx = idx
            self.current_var = label
            self.update_plot()
        except ValueError:
            pass
    
    def on_colormap_changed(self, label):
        """Handle colormap selection"""
        self.current_cmap = label
        self.update_plot()
    
    def on_axis_changed(self, label):
        """Handle slice axis change"""
        axis_map = {"X": 0, "Y": 1, "Z": 2}
        self.slice_axis = axis_map[label]
        
        # Update slider range
        max_idx = self.plotfile.grid_dims[self.slice_axis] - 1
        self.slider.valmin = 0
        self.slider.valmax = max_idx
        self.slider.set_val(min(self.current_slice_idx, max_idx))
        self.slider.ax.set_xlim(0, max_idx)
        
        self.current_slice_idx = int(self.slider.val)
        self.update_plot()
    
    def on_slice_changed(self, val):
        """Handle slice index change"""
        self.current_slice_idx = int(val)
        self.update_plot()
    
    def update_plot(self):
        """Update the visualization"""
        if self.updating:
            return
        
        self.updating = True
        try:
            print(f"Loading {self.current_var}...")
            
            # Load data
            data = self.plotfile.read_variable_data(self.current_var)
            self.current_data_3d = data  # Store for line plots
            
            # Extract 2D slice
            if self.slice_axis == 0:  # X
                slice_data = data[:, :, self.current_slice_idx]
                xlabel, ylabel = "Y", "Z"
            elif self.slice_axis == 1:  # Y
                slice_data = data[:, self.current_slice_idx, :]
                xlabel, ylabel = "X", "Z"
            else:  # Z
                slice_data = data[self.current_slice_idx, :, :]
                xlabel, ylabel = "X", "Y"
            
            self.current_slice_data = slice_data
            self.xlabel = xlabel
            self.ylabel = ylabel
            
            print(f"  Slice shape: {slice_data.shape}, min: {slice_data.min():.3e}, max: {slice_data.max():.3e}")
            
            # Check for invalid data
            if np.all(slice_data == 0):
                print("Warning: All data is zero")
            if np.any(np.isnan(slice_data)):
                print("Warning: Data contains NaN values")
                slice_data = np.nan_to_num(slice_data, nan=0.0)
            if np.any(np.isinf(slice_data)):
                print("Warning: Data contains Inf values")
                slice_data = np.nan_to_num(slice_data, posinf=0.0, neginf=0.0)
            
            # Clear and plot
            self.ax_main.clear()
            
            # Use actual min/max for colormap
            vmin, vmax = slice_data.min(), slice_data.max()
            if vmin == vmax:
                vmin, vmax = vmin - 1, vmax + 1
            
            # Don't use extent - let matplotlib handle pixel coordinates directly
            # Use aspect='equal' to preserve the physical aspect ratio of the grid cells
            im = self.ax_main.imshow(slice_data, origin='lower', aspect='equal', 
                                     interpolation='nearest', cmap=self.current_cmap,
                                     vmin=vmin, vmax=vmax)
            self.ax_main.set_xlabel(xlabel)
            self.ax_main.set_ylabel(ylabel)
            
            axis_names = ["X", "Y", "Z"]
            max_idx = self.plotfile.grid_dims[self.slice_axis] - 1
            title = f"{self.current_var} ({axis_names[self.slice_axis]}={self.current_slice_idx}/{max_idx})"
            self.ax_main.set_title(title)
            
            # Colorbar - update existing or create new
            if hasattr(self, 'colorbar') and self.colorbar is not None:
                # Update existing colorbar
                self.colorbar.update_normal(im)
            else:
                # Create new colorbar
                self.colorbar = self.fig.colorbar(im, ax=self.ax_main)
            
            # Add statistics
            stats = f"min: {slice_data.min():.3e}  max: {slice_data.max():.3e}  mean: {slice_data.mean():.3e}"
            self.ax_main.text(0.5, -0.15, stats, transform=self.ax_main.transAxes,
                             ha='center', fontsize=9, family='monospace')
            
            plt.draw()
        finally:
            self.updating = False
    
    def on_mouse_move(self, event):
        """Handle mouse movement to show value tooltip"""
        if event.inaxes != self.ax_main or self.current_slice_data is None:
            if self.annotation:
                self.annotation.set_visible(False)
                self.fig.canvas.draw_idle()
            return
        
        # Get pixel coordinates
        x, y = int(event.xdata + 0.5), int(event.ydata + 0.5)
        
        # Check bounds
        if 0 <= y < self.current_slice_data.shape[0] and 0 <= x < self.current_slice_data.shape[1]:
            value = self.current_slice_data[y, x]
            
            # Create or update annotation
            if self.annotation is None:
                self.annotation = self.ax_main.annotate(
                    '', xy=(0, 0), xytext=(10, 10),
                    textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8),
                    fontsize=9, family='monospace'
                )
            
            self.annotation.set_text(f'{self.xlabel}={x}, {self.ylabel}={y}\nValue: {value:.6e}')
            self.annotation.xy = (event.xdata, event.ydata)
            self.annotation.set_visible(True)
            self.fig.canvas.draw_idle()
    
    def on_mouse_click(self, event):
        """Handle mouse click to show line plot along direction"""
        if event.inaxes != self.ax_main or self.current_data_3d is None:
            return
        
        # Get pixel coordinates
        x, y = int(event.xdata + 0.5), int(event.ydata + 0.5)
        
        # Check bounds
        if not (0 <= y < self.current_slice_data.shape[0] and 0 <= x < self.current_slice_data.shape[1]):
            return
        
        # Create line plot window
        if self.line_plot_window is not None:
            plt.close(self.line_plot_window)
        
        self.line_plot_window = plt.figure(figsize=(10, 6))
        
        # Create 3 subplots for X, Y, Z directions
        axis_names = ['X', 'Y', 'Z']
        
        # Map pixel coordinates to 3D indices based on current view
        if self.slice_axis == 0:  # X-axis slice (viewing Y-Z plane)
            # x -> Y, y -> Z, fixed X
            indices = [self.current_slice_idx, x, y]  # [X, Y, Z]
        elif self.slice_axis == 1:  # Y-axis slice (viewing X-Z plane)
            # x -> X, y -> Z, fixed Y
            indices = [x, self.current_slice_idx, y]  # [X, Y, Z]
        else:  # Z-axis slice (viewing X-Y plane)
            # x -> X, y -> Y, fixed Z
            indices = [x, y, self.current_slice_idx]  # [X, Y, Z]
        
        for i, axis_name in enumerate(axis_names):
            ax = self.line_plot_window.add_subplot(1, 3, i + 1)
            
            # Extract 1D line along this axis
            if i == 0:  # Along X
                line_data = self.current_data_3d[indices[2], indices[1], :]
                x_coords = np.arange(self.plotfile.grid_dims[0])
                fixed_info = f'Y={indices[1]}, Z={indices[2]}'
            elif i == 1:  # Along Y
                line_data = self.current_data_3d[indices[2], :, indices[0]]
                x_coords = np.arange(self.plotfile.grid_dims[1])
                fixed_info = f'X={indices[0]}, Z={indices[2]}'
            else:  # Along Z
                line_data = self.current_data_3d[:, indices[1], indices[0]]
                x_coords = np.arange(self.plotfile.grid_dims[2])
                fixed_info = f'X={indices[0]}, Y={indices[1]}'
            
            # Plot
            ax.plot(x_coords, line_data, 'b-', linewidth=1.5)
            ax.grid(True, alpha=0.3)
            ax.set_xlabel(f'{axis_name} index')
            ax.set_ylabel(self.current_var)
            ax.set_title(f'Along {axis_name} ({fixed_info})')
            
            # Mark the current position
            current_pos = indices[i]
            ax.axvline(current_pos, color='r', linestyle='--', alpha=0.7, label=f'Current: {axis_name}={current_pos}')
            ax.legend(fontsize=8)
        
        self.line_plot_window.suptitle(f'{self.current_var} - Line plots through clicked point', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.show()
    
    def show(self):
        """Display the viewer"""
        plt.show()


def main():
    import sys
    
    # Check if plotfile directory provided as argument
    if len(sys.argv) < 2:
        print("Usage: python pltview.py <plotfile_directory>")
        sys.exit(1)
    
    plotfile_dir = sys.argv[1]
    
    if not os.path.isdir(plotfile_dir):
        print(f"Error: '{plotfile_dir}' is not a valid directory")
        sys.exit(1)
    
    try:
        app = PLTViewApp(plotfile_dir)
        app.show()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()