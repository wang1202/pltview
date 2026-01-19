# Installation Guide

## Quick Start

On your Linux server or macOS system with X11:

```bash
pip install git+https://github.com/wang1202/pltview.git
```

The installation will automatically compile the C version.

## Prerequisites

### Linux
**Debian/Ubuntu:**
```bash
sudo apt-get install libx11-dev libxt-dev libxaw7-dev libxmu-dev gcc
```

**RHEL/CentOS/Fedora:**
```bash
sudo yum install libX11-devel libXt-devel libXaw-devel libXmu-devel gcc
```

### macOS
Install XQuartz from https://www.xquartz.org/

## Usage

```bash
# Run the viewer
pltview plt00100
```

## What Happens During Installation

1. **pip install** downloads the package
2. **setup.py** checks for X11 libraries (required)
3. Compiles `pltview.c` → `pltview` binary
4. Installs the binary to your PATH

### Example Installation Output

```
Building pltview (C version)...
============================================================
Using X11 from: /usr/include
Running: gcc -O3 -Wall -march=native -I/usr/include -o pltview pltview.c ...
✓ pltview built successfully!
============================================================
Installing pltview to /home/user/.local/bin/pltview
✓ Installation complete!
```

## Verification

After installation:

```bash
# Check if pltview is in PATH
which pltview

# Run it
pltview plt00100
```

## Performance

On a typical Linux server with a 320×512×100 dataset:
- **Load + render**: ~0.1 seconds
- **Ultra-fast** for interactive browsing of large AMReX datasets
