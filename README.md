# Fast Editing of Singularities in Field-Aligned Stripe Patterns [SIGGRRAPH Asia 2022]

This is the implementation of the SIGGRAPH Asia 2022 paper "[Fast Editing of Singularities in Field-Aligned Stripe Patterns"](https://yutanoma.com/projects/singularity-editing), by Yuta Noma, Nobuyuki Umetani, and Yoshihiro Kawahara.

## Installation

This project is built upon [libigl](http://libigl.github.io/) and [directional](https://avaxman.github.io/Directional/).

The only dependencies are STL, Eigen, [libigl](http://libigl.github.io/libigl/), [directional](https://avaxman.github.io/Directional/), and the dependencies
of the `igl::opengl::glfw::Viewer` (OpenGL, glad and GLFW).
The CMake build system will automatically download libigl and its dependencies using
[CMake FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html),
thus requiring no setup on your part.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `singularity-editing` binary.

## Run

From within the `build` directory just issue:

    ./singularity-editing

