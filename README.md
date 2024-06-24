# Fast Editing of Singularities in Field-Aligned Stripe Patterns [SIGGRRAPH Asia 2022]

Public code release for [Fast Editing of Singularities in Field-Aligned Stripe Patterns](https://yutanoma.com/projects/singularity-editing).

**Fast Editing of Singularities in Field-Aligned Stripe Patterns**<br>
Yuta Noma, Nobuyuki Umetani, Yoshihiro Kawahara. <br>
Proceedings of SIGGRAPH Asia 2022 (Conference track) <br>
[[Project page](https://yutanoma.com/projects/singularity-editing)]

## Build the code

**Unix-like machines**: configure (with cmake) and compile
```
cd /path/to/directory
mkdir build
cd build
cmake ..
make -j6
```

**Windows / Visual Studio**

Install CMake, and use either the CMake GUI or the command line interface (as on unix) to generate a Visual Studio solution.  Build the solution with Visual Studio.

## Run the code

```
./singularity-editing ../data/models/bunny.obj
```

## UI options

### Move singularities

By selecting "move singularities", you can click on a singularity and move it to another position.

### Add singularities

By selecting "add singularity pairs" and clicking on two different faces, you can add a pair of singularities.

### Remove singularities

By selecting "remove singularities" and clicking a pair of positive (red) and negative (blue) singularities, you can annihilate them.

## Contact

If any issues or questions, please contact yn.devilstick@gmail.com.


