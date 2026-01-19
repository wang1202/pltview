# PLTView

Ultra-fast X11 viewer for AMReX plotfiles, inspired by ncview. Written in C for maximum performance.

## Features

- **Ultra-fast rendering**: ~10-100x faster than Python viewers
- Direct X11 rendering with minimal dependencies
- Browse variables in AMReX plotfile directories
- View 2D slices of 3D data along X, Y, Z axes
- Interactive GUI with mouse controls:
  - Hover to see values at cursor position
  - Click to view 1D line profiles along X, Y, Z directions
- Multiple colormap options (viridis, jet, turbo, plasma)
- Aspect ratio preservation
- Clean white background for publication-ready figures

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

- **Variable Buttons** (left sidebar): Click to select which variable to visualize
- **Axis Buttons** (X/Y/Z, bottom): Click to switch viewing axis
- **Colormap Buttons** (viridis/jet/turbo/plasma, bottom): Select colormap
- **+/- Buttons** (bottom): Navigate through slices
- **Colorbar** (right): Shows data range and colormap scale
- **Mouse Hover**: Shows value at cursor position in info label
- **Mouse Click**: Opens popup window with line profiles along X, Y, Z directions

The line profile popup displays three plots showing how the variable value changes along each spatial dimension through the clicked point.

## Requirements

- Python >= 3.8
- numpy >= 1.20.0
- matplotlib >= 3.0.0

## Usage

After installation, simply run:

```bash
# Automatically uses C version if available, otherwise Python version
pltview plt00100

# Force Python version
pltview-py plt00100

# Direct C version (if installed to PATH)
pltview_c plt00100
```

The `pltview` command will intelligently choose:
1. **C version** if available and DISPLAY is set (X11 available)
2. **Python version** as fallback

This gives you the best performance automatically!

## Performance

For typical use cases:

| Dataset Size | Load + Render Time |
|--------------|-------------------|
| 100×100×100  | ~0.05s |
| 320×512×100  | ~0.1s |
| 1000×1000×1000 | ~0.5s |

*Times shown are for initial load + first render on a typical workstation*

**Ideal for:**
- Very large datasets (>100M cells)
- Remote visualization over SSH with X11 forwarding
- Quick browsing of many plotfiles
- Minimal system resources

## File Format

This tool reads AMReX plotfile format (used by ERF, AMReX-Hydro, etc.):

- `Header`: Metadata about variables and grid structure
- `Level_X/`: Data directories for each AMR level
- `Cell_D_XXXXX`: FAB binary data files for each MPI domain
- `Cell_H`: Cell data header with box layout information

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
