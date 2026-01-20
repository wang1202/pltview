# pltview

A lightweight X11 viewer for AMReX plotfiles, inspired by ncview and built with the assistance of Claude Sonnet 4.5.

![Example Screenshot](Example.png)

## Features

- **Direct X11 rendering** with minimal dependencies (10-100x faster than Python alternatives)
- **Multi-level AMR support**: Automatically detects and visualizes multiple refinement levels
- **Interactive 3D slicing**: View 2D slices of 3D data along X, Y, Z axes with wrap-around navigation
- **Mouse interaction**:
  - Hover to see values at cursor position
  - Click to view 1D line profiles along X, Y, Z directions in popup window
- **Multiple colormap options**: viridis, jet, turbo, plasma, hot, cool, gray, magma
- **Level handling**: Preserves slice position when switching between AMR levels
- **Dynamic grid adaptation**: Automatically adjusts to different grid dimensions per level
- **Variables supported**: Displays all available variables (up to 128) in the sidebar

## Installation

### Via pip (Recommended)

```bash
pip install git+https://github.com/wang1202/pltview.git
```

### Editable Install (For Development)

If you want to modify the code and have changes take effect immediately:

```bash
git clone https://github.com/wang1202/pltview.git
cd pltview
pip install -e .
```

After editable install, you can modify `pltview.c` and rebuild:

```bash
make
# Changes take effect immediately - pltview command uses the updated binary
```

### From Source (Manual Build)

```bash
git clone https://github.com/wang1202/pltview.git
cd pltview
make
./pltview plt00100
```

**Prerequisites:**

- **macOS**: Install XQuartz from https://www.xquartz.org/
- **Linux**: Install X11 development libraries:
  - Debian/Ubuntu: `sudo apt-get install libx11-dev libxt-dev libxaw7-dev libxmu-dev`
  - RHEL/CentOS: `sudo yum install libX11-devel libXt-devel libXaw-devel libXmu-devel`

## Usage

After installation:

```bash
pltview plt00100
```

## Controls

**GUI Layout:**

- **Left sidebar**: Variable selection buttons (all available variables, up to 128 supported)
- **Main canvas**: Data visualization with white background and aspect ratio preservation
- **Right colorbar**: Data range and colormap scale
- **Bottom controls** (organized in 2 columns):
  - **Column 1**: Axis buttons (X/Y/Z) and navigation (+/-)
  - **Column 2**: Level selection (Level 0/Level 1/...) and colormap buttons

**Mouse Interaction:**

- **Hover**: Shows value at cursor position in info label at top
- **Click**: Opens popup window with line profiles along X, Y, Z directions

**Buttons:**

- **Variable Buttons**: Select which variable to visualize
- **X/Y/Z Buttons**: Switch viewing axis (perpendicular to slice)
- **Level Buttons**: Switch between AMR refinement levels (appears when multiple levels detected)
- **Colormap Buttons**: Choose from 8 colormaps (viridis/jet/turbo/plasma/hot/cool/gray/magma)
- **+/- Buttons**: Navigate through slices with wrap-around (layer 1 → - → last layer)
- **Keyboard**: Arrow keys also navigate slices

**Line Profile Popup:**
The popup window displays three graphs showing how the variable value changes along each spatial dimension (X, Y, Z) through the clicked point, with proper axis labels and tick marks.

## Requirements

- **C Compiler**: gcc or clang
- **X11 Libraries**:
  - macOS: XQuartz (https://www.xquartz.org/)
  - Linux: libX11, libXt, libXaw, libXmu development packages
- **Python**: >= 3.6 (for pip installation wrapper)

**Runtime dependencies**: None beyond X11 (no Python runtime dependencies for the viewer itself)

## File Format

This tool reads AMReX plotfile format (used by ERF, AMReX-Hydro, etc.):

- `Header`: Metadata about variables and grid structure
- `Level_0/`, `Level_1/`, ...: Data directories for each AMR refinement level
- `Cell_D_XXXXX`: FAB binary data files for each MPI domain
- `Cell_H`: Cell data header with box layout and FabOnDisk mapping

**Multi-level Support:**

- Automatically detects available AMR levels by scanning for `Level_X` directories
- Handles varying grid dimensions across different refinement levels
- Preserves slice position when switching levels (clamped to valid range if needed)

Each Cell_D file contains a FAB (Fortran Array Box) header followed by binary double-precision floating-point data in Fortran (column-major) order.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please consider citing:

```bibtex
@software{pltview,
  author = {Wang, Aaron},
  title = {pltview: A simple viewer for AMReX plotfiles},
  year = {2025},
  url = {https://github.com/wang1202/pltview}
}
```
