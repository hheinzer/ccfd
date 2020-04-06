![](doc/ccfd.svg)

This code is intended to become a drop-in replacement for `cfdfv`, a CFD code written in Fortran by the [Institute of Aerodynamics and Gas Dynamics](https://www.iag.uni-stuttgart.de/en/) at the University of Stuttgart for a CFD programming course. This code itself is not available online, as far as I know, but it features a similar code structure to [FLEXI](https://www.flexi-project.org/).

The program uses the [CGNS](https://cgns.github.io/) library version 4.1.1, for storing the calculation results.

The logo is heavily inspired by the software over at [suckless.org](https://suckless.org/).

![](doc/comparison.svg)

# Dependencies

- `git`
- `make`, and `cmake`
- `gcc`
- `gnuplot` (optional, for displaying calculation residuals, available [here](http://www.gnuplot.info/))
- `gmsh` (optional, for mesh generation, available [here](http://gmsh.info/))
- `paraview` (optional, for post-processing the results, available [here](https://www.paraview.org/))

# Installation

The installation process is easiest on Linux, but possible on MacOS and Windows. For best performance you should build an run it on your native OS.

## Linux

First make sure that all necessary dependencies are all installed. These can usually be obtained through your distributions package manager, on Arch based systems the following command should suffice
```
# pacman -S git base-devel cmake
```
For an Ubuntu like system the following command should do the trick
```
# apt-get install git build-essential cmake
```

Next, navigate to the directory where you want to keep `ccfd` files and clone the git repository
```
$ cd path/to/directory
$ git clone https://github.com/hhh95/ccfd.git
```
Then, you can enter change to the `ccfd` directory and build the program using `make`
```
$ cd ccfd
$ make
```
This should compile the code and create the two folders `obj` and `bin`, the last one containing the `ccfd` executable.

Continue with [Usage](#usage-under-linux)

## MacOS

It should work basically the same as with Linux, but detailed instructions will follow soon...

## Windows

It should work basically the same as with Linux, but detailed instructions will follow soon...

# Usage under Linux

As a first step you should add the `ccfd` executable to your path. You can do so by running

```
$ source ccfdrc
```
Next, navigate to the `calc` folder. Here you will find example case files, contained in folder. First, let us try the Riemann problems. Navigate to the `riemann` folder with
```
$ cd riemann
```
Here you will find the SOD test case, as well as different versions of the case. Start a calculation with
```
$ ccfd sod.ini
```
and observe the output. There should be four new files. The initial condition of the case and the calculation results at t = 0.25 s. They should all be `.csv` files. You can examine them with any spreadsheet program you like. Alternatively you can use `paraview`, a free post-processing program, that can visualize 1D, 2D, and 3D data. For most distributions you can install it from the package manager. On Arch Linux, run
```
# pacman -S paraview
```

On Ubuntu, the Paraview program in the repositories does not read CGNS files for some reason. You will need to compile [Paraview from source](https://github.com/Kitware/ParaView/blob/master/Documentation/dev/build.md), with CGNS enabled:
```
# apt-get install git cmake build-essential libgl1-mesa-dev libxt-dev qt5-default libqt5x11extras5-dev libqt5help5 qttools5-dev qtxmlpatterns5-dev-tools libqt5svg5-dev python3-dev python3-numpy libopenmpi-dev libtbb-dev ninja-build
$ git clone https://gitlab.kitware.com/paraview/paraview.git
$ mkdir paraview_build
$ cd paraview
$ git checkout v5.8.0
$ git submodule update --init --recursive
$ cd ../paraview_build
$ cmake -GNinja -DPARAVIEW_ENABLE_PYTHON=ON -DPARAVIEW_USE_MPI=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB -DCMAKE_BUILD_TYPE=Release -DPARAVIEW_ENABLE_CGNS=ON ../paraview
$ ninja
```


Next open the Paraview program
```
$ paraview &
```
Now, click on *File*->*Open* and select both sets of `.csv` files. Because we will be looking at 1D data, we need to change from *Render View* to *Line Chart View*. In the top right of the viewing area, click on the `X` button. Now select *Line Chart View* from the list. You should now see an empty grid with an x-, and a y-axis. Now in the *Pipeline Browser* to the left, click on the eye icons of the `.csv` files you have loaded in. A plot should appear in on the axis grid. It will probably show the initial state. In the top bar, click on the play button. Now, the final state should be shown. You will see the analytical, or exact, solution, as well as the numerical solution.

For more information on the theory, maybe have a look at [Wikipedia](https://en.wikipedia.org/wiki/Sod_shock_tube).

The procedure for running the other cases is the same. However, if the solution data is 2D, then you do not need to switch to *Line Chart View*. You can stay in *Render View* and click *Apply*, once you have loaded the solution file. The 2D CGNS output files will usually have more than just the solution. You can can load everything at once by selecting the file that hast `_Master` in its name.
