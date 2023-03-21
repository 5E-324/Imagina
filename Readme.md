# Imagina

Imagina is a fast fractal renderer & explorer.

_The project is being rewritten. This repository may be renamed and replaced by the rewritten version when it's available._

## Controls

* Left mouse button
    * Drag - Pan view
    * Double click - Find nearby feature and open Feature Finder dialog
* Mouse scroll wheel - Zoom
* A - Increase color density
* S - Decrease color density
* D - Increase iteration limit
* F - Decrease iteration limit
* E - Cycle color forward
* R - Cycle color backward

## Build
### Dependencies
* mpir
* libpng

### CMake variables
* INCLUDE_PATH Path to include files.
* LIBRARY_PATH Path to libraries.
* MPIR_USE_DLL Use dll version of MPIR. Default: ON (To use statically linked version, set this to OFF)