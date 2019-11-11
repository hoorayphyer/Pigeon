# BUICK: Build Usher In CmaKe
## `cd build && cmake ..` no more! Let's `buick .`!

`buick` is a snippet to automate running `cmake` in a user-defined build directory to prevent the notorious cmake build cache files contamination.

## Install
Download and `chmod` the script `buick`. To make it system-wide available, just create an alias in your shell rc file.

## Usage
1. `buick` with no arguments will `cd` to the build directory and `make`.
2. `buick <path_to_CMakeLists> <cmake_commandline_arguments>` will create `./build` and run `cmake <cmake_commandline_arguments>` in it. `<path_to_CMakeLists>` can be relative.
3. The above command will create `./.buick` in which the build directory and last successful cmake command are cached.
4. `buick --last-cmake` runs the cached cmake command.
5. `buick --build-dir=<your_build_dir>` uses a different directory than the default `./build`.
6. `buick --clean` will `rm -rf ` the cached build directory if it exists. Passing `--clean` will make `buick` ignore all other arguments.

## Caveats
* I only included gnu make. For other build tools, please modify the script as needed.
* Relative paths used in cmake commandline options may not be properly supported.
