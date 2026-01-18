# PLTView

A simple viewer for AMReX plotfiles, inspired by ncview.

## Features

- Browse variables in AMReX plotfile directories
- View 2D slices of 3D data
- Navigate through different slice axes (X, Y, Z)
- Interactive matplotlib-based GUI with XQuartz support
- Display data statistics (min, max, mean)
- Hover tooltip showing values at cursor position
- Click to view 1D line plots along X, Y, Z directions
- Multiple colormap options

## Installation

### From Git Repository

```bash
# Clone the repository
git clone https://github.com/yourusername/pltview.git
cd pltview

# Install in editable mode
pip install -e .
```

### Development Installation

```bash
# With virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On macOS/Linux
pip install -e .
```

## Requirements

- Python >= 3.8
- numpy >= 1.20.0
- matplotlib >= 3.0.0

## Usage

After installation, you can run pltview from anywhere:

```bash
# Run with plotfile directory
pltview plt00100

# Or use the Python module directly
python -m pltview plt00100
```

## Controls

- **Variable Buttons** (left panel, top): Click to select which variable to visualize
- **Colormap Buttons** (left panel, middle): Select colormap (viridis, turbo, jet, etc.)
- **Slice Axis** (bottom left): Choose X, Y, or Z axis for slicing  
- **Slice Slider** (bottom): Drag to navigate through slices along the chosen axis
- **Mouse Hover**: Hover over the plot to see coordinates and values
- **Mouse Click**: Click on any point to open line plots along X, Y, Z directions
- **Matplotlib Toolbar**: Use built-in tools to zoom, pan, and save images

## File Format

This tool reads AMReX plotfile format (used by ERF, AMReX-Hydro, etc.):
- `Header`: Metadata about variables and grid structure
- `Level_X/`: Data directories for each AMR level
- `Cell_D_XXXXX`: FAB binary data files for each MPI domain
- `Cell_H`: Cell data header with box layout information

Each Cell_D file contains a FAB (Fortran Array Box) header followed by binary double-precision floating-point data in Fortran (column-major) order.

## Example Data

Your plotfile contains 13 variables on a 140×80×100 grid:
- density, rhoadv_0
- x_velocity, y_velocity, z_velocity  
- temp, theta
- pres_hse, dens_hse, pressure, pert_pres
- z_phys, mapfac
