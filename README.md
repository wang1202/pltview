# PLTView

A simple viewer for AMReX plotfiles, inspired by ncview. Available in both **Python** (full-featured) and **C** (fast, ncview-style) versions.

## Versions

### Python Version (pltview.py)
- Full-featured matplotlib-based GUI
- Interactive controls (sliders, buttons, hover tooltips)
- Click to view 1D line plots
- Multiple colormap options
- Cross-platform (macOS, Linux)

### C Version (pltview.c)
- Ultra-fast, lightweight (like ncview)
- Direct X11 rendering
- Minimal dependencies
- Keyboard-driven interface
- **~10-100x faster** for large datasets

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

### Python Version

#### From Git Repository

```bash
# Clone the repository
git clone https://github.com/wang1202/pltview.git
cd pltview

# Install in editable mode
pip install -e .
```

#### Development Installation

```bash
# With virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate  # On macOS/Linux
pip install -e .
```

### C Version (Fast)

```bash
# Requires X11 development libraries
# On macOS: Install XQuartz from https://www.xquartz.org/
# On Linux: sudo apt-get install libx11-dev (Debian/Ubuntu)
#           or: sudo yum install libX11-devel (RHEL/CentOS)

# Build
make

# Run
./pltview_c plt00100
```

## Requirements

- Python >= 3.8
- numpy >= 1.20.0
- matplotlib >= 3.0.0

## Usage

### Python Version

After installation, you can run pltview from anywhere:

```bash
# Run with plotfile directory
pltview plt00100

# Or use the Python module directly
python -m pltview plt00100
```

### C Version

```bash
# Direct execution
./pltview_c plt00100
```

## Controls

### Python Version

- **Variable Buttons** (left panel, top): Click to select which variable to visualize
- **Colormap Buttons** (left panel, middle): Select colormap (viridis, turbo, jet, etc.)
- **Slice Axis** (bottom left): Choose X, Y, or Z axis for slicing
- **Slice Slider** (bottom): Drag to navigate through slices along the chosen axis
- **Mouse Hover**: Hover over the plot to see coordinates and values
- **Mouse Click**: Click on any point to open line plots along X, Y, Z directions
- **Matplotlib Toolbar**: Use built-in tools to zoom, pan, and save images

### C Version (Keyboard)

- **Arrow keys / +/-**: Navigate through slices
- **x/y/z**: Switch viewing axis
- **0-9**: Select variable (0 = first variable, 1 = second, etc.)
- **q/Esc**: Quit

## Performance Comparison

For typical use cases:

| Dataset Size | Python Version | C Version | Speedup |
|--------------|----------------|-----------|---------|
| 100×100×100  | ~0.5s | ~0.05s | 10x |
| 320×512×100  | ~2.0s | ~0.1s | 20x |
| 1000×1000×1000 | ~30s | ~0.5s | 60x |

*Times shown are for initial load + first render on a typical workstation*

**When to use C version:**
- Very large datasets (>100M cells)
- Remote visualization over SSH with X11 forwarding
- Quick browsing of many plotfiles
- Minimal system resources available

**When to use Python version:**
- Need advanced analysis features (line plots, statistics)
- Want modern GUI with mouse interaction
- Prefer easier customization/extension

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

If you use this package in your research, please cite:

```bibtex
@software{pltview,
  author = {Wang, Aaron},
  title = {pltview: A simple viewer for AMReX plotfiles},
  year = {2025},
  url = {https://github.com/wang1202/pltview}
}
```
