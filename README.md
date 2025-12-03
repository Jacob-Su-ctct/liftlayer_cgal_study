# LiftLayer CGAL Study

This repository contains CGAL examples, visualization tools, and build configurations for our team's CGAL study and development work.

## Overview

- **CGAL Integration**: Uses CGAL as a Git submodule for clean dependency management
- **Custom Examples**: Sample code demonstrating CGAL functionality
- **Visualization Tools**: Python scripts for viewing OFF files and mesh data
- **Docker Support**: Build configurations for embedded device deployment

## Prerequisites

- **CMake** 3.12 or higher
- **C++ Compiler** with C++17 support (GCC 9+, Clang 10+, or MSVC 2019+)
- **Python 3** (for visualization scripts)
  - `trimesh`
  - `matplotlib`
  - `numpy`

### Install Python Dependencies

```bash
pip3 install trimesh matplotlib numpy
```

## Quick Start

### Clone the Repository

```bash
# Clone with submodules (recommended)
git clone --recursive https://github.com/YOUR_USERNAME/liftlayer_cgal_study.git
cd liftlayer_cgal_study

# Or if already cloned without --recursive
git clone https://github.com/YOUR_USERNAME/liftlayer_cgal_study.git
cd liftlayer_cgal_study
git submodule update --init --recursive
```

### Build Examples

```bash
# Create build directory
mkdir build && cd build

# Configure with CMake
cmake -DCGAL_WITH_examples=ON -DCGAL_DIR=../cgal ..

# Build
make -j$(nproc)

# Run examples
./examples/example_name
```

## Repository Structure

```
liftlayer_cgal_study/
├── cgal/                    # CGAL library (Git submodule)
├── examples/                # Custom CGAL examples
├── scripts/                # Utility scripts
├── docker/                 # Docker configurations
├── CMakeLists.txt          # Root CMake configuration
└── README.md
```

## Visualization Tools

### View Multiple OFF Files

```bash
# View one or more OFF mesh files
python3 scripts/view_multiple_off.py file1.off file2.off

# List available OFF files in current directory
python3 scripts/view_multiple_off.py
```

Features:
- Display multiple meshes with different colors
- Automatic equal aspect ratio (1:1:1)
- Support for both mesh faces and point cloud curves

## Docker Build

Build CGAL for embedded devices:

```bash
# Build Docker image
docker build --platform linux/arm/v7 -t cgal-armv7-builder-buster-glibc2.28 .

# Run container
docker run -it cgal-embedded
```

## Development Workflow

### Update CGAL Version

```bash
cd cgal
git fetch
git checkout v6.1  # or desired version/branch
cd ..
git add cgal
git commit -m "Update CGAL to v6.1"
git push
```

### Add New Examples

1. Create your example in `examples/` directory
2. Update `examples/CMakeLists.txt` if needed
3. Build and test
4. Commit and push

```bash
cd examples
# Create your example
nano my_new_example.cpp

cd ../build
make
./examples/my_new_example
```

## Troubleshooting

### CGAL Directory is Empty

If the `cgal/` directory is empty after cloning:

```bash
git submodule update --init --recursive
```

### CMake Can't Find CGAL

Ensure you're pointing to the submodule:

```bash
cmake -DCGAL_DIR=../cgal ..
```

Or set it as an environment variable:

```bash
export CGAL_DIR=/absolute/path/to/liftlayer_cgal_study/cgal
cmake ..
```

### Build Errors

Make sure you have all dependencies:

```bash
# Ubuntu/Debian
sudo apt-get install build-essential cmake libgmp-dev libmpfr-dev libboost-all-dev

# Fedora/RHEL
sudo dnf install gcc-c++ cmake gmp-devel mpfr-devel boost-devel

# macOS
brew install cmake gmp mpfr boost
```

## Useful CGAL Resources

- [CGAL Documentation](https://doc.cgal.org/)
- [CGAL User Manual](https://doc.cgal.org/latest/Manual/index.html)
- [CGAL GitHub Repository](https://github.com/CGAL/cgal)

## Contributing

1. Create a feature branch
2. Make your changes
3. Test thoroughly
4. Submit a pull request

## License

This project uses CGAL which is licensed under GPL/LGPL. Please refer to the CGAL license for details.

## Team Notes

- Main CGAL repository: https://github.com/CGAL/cgal
- Current CGAL version tracked: Check `cgal/` submodule commit
- Build issues? Check the troubleshooting section above

---

**Note**: The `cgal/` directory is a Git submodule. Always clone with `--recursive` or run `git submodule update --init --recursive` after cloning.
