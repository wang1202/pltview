# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pltview is a lightweight X11 viewer for AMReX plotfiles written primarily in C with a Python wrapper for installation. It's designed for ultra-fast visualization of large scientific simulation datasets with minimal dependencies.

## Architecture

### Core Components
- **pltview.c**: Single-file C application (~6400 lines, no header files) containing all visualization logic, X11/Athena widget GUI, AMReX plotfile parsing, and data rendering
- **pltview_entry.py**: Python wrapper script that locates and executes the compiled C binary
- **setup.py**: Custom setuptools build system that compiles the C binary during installation

### Key Features
- Multi-level AMR (Adaptive Mesh Refinement) support with automatic level detection
- Multi-timestep navigation for time series data  
- Interactive 3D slicing along X/Y/Z axes
- **Enhanced Quiver (vector field) visualization** with comprehensive controls:
  - Component selection dialog with smart defaults
  - Adjustable arrow density (1-5 levels)
  - Variable arrow scale (0.2x-3.0x)
  - Color options (black, white, red, blue)
- SDM (Super Droplet Method) particle data visualization mode
- Statistical analysis (histograms, profiles, time series)
- Multiple colormap support (viridis, jet, turbo, plasma, hot, cool, gray, magma)

## Development Commands

### Building
```bash
# Quick compile during development (produces ./pltview_c)
make

# Editable development install (recommended — produces ./pltview, wired to `pltview` CLI command)
pip install -e .
```

### Testing
There is no test suite. Verify changes manually by running `pltview` against AMReX plotfile directories.

### Running the Application
```bash
# Single plotfile
pltview plt00100

# Multi-timestep mode
pltview /path/to/output plt

# SDM particle mode
pltview --sdm plt00100
```

## Code Structure

### C Application (pltview.c)
- **Data Structures**: PlotfileData, Box, RGB, PlotData, PopupData, QuiverData for managing AMReX data and visualization state
- **AMReX Parser**: Functions for reading plotfile headers, FAB data, and multi-level AMR structures
- **X11 GUI**: Athena widget-based interface with buttons, canvas, and popup windows
- **Rendering Engine**: Direct pixel manipulation for fast data visualization
- **Colormap System**: Multiple built-in colormaps with dynamic range adjustment
- **Quiver System**: Vector field overlay with automatic component defaults (x_velocity/y_velocity for Z plane, etc.)
- **SDM Mode**: Specialized particle data histogram analysis

### Python Components
- **pltview_entry.py**: Entry point wrapper that finds and executes the C binary
- **setup.py**: Custom build commands (BuildC, InstallC, DevelopC) that handle X11 dependency checking and C compilation

## File Format Support

Reads AMReX plotfile format:
- `Header`: Variable metadata and grid structure
- `Level_*/`: AMR refinement level directories
- `Cell_D_*`: FAB binary data files (double-precision, Fortran column-major order)
- `Cell_H`: Box layout and domain mapping
- `super_droplets_moisture/`: Particle data for SDM mode

## Development Notes

### Prerequisites
- X11 development libraries (libX11-dev, libXt-dev, libXaw7-dev, libXmu-dev on Linux)
- XQuartz on macOS
- gcc compiler

### Key Considerations
- No runtime dependencies beyond X11 (no Python libraries needed at runtime)
- Direct pixel manipulation for performance (10-100x faster than Python alternatives)
- Handles large datasets (320×512×100+ grid points) interactively
- Preserves slice positions when switching AMR levels
- Memory management for large multi-timestep datasets

### Build System Details
- `make` produces `pltview_c`; `pip install -e .` produces `pltview` — both from pltview.c but with different output names
- Custom setuptools commands (BuildC, InstallC, DevelopC, EditableWheel in setup.py) detect X11 paths automatically
- After `pip install -e .`, recompile with `make` or re-run `pip install -e .` to pick up C code changes
- Compiler flags: `-O3 -Wall -march=native`; linked against `-lX11 -lXt -lXaw -lXmu -lm`
- macOS requires XQuartz; include/lib paths at `/opt/X11/`