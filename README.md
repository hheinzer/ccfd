![](doc/ccfd.svg)

This code is intended to become a drop-in replacement for `cfdfv`, a CFD code written in Fortran by the [Institute of Aerodynamics and Gas Dynamics](https://www.iag.uni-stuttgart.de/en/) at the University of Stuttgart for a CFD programming course. This code itself is not available online, as far as I know, but it features a similar code structure to [FLEXI](https://www.flexi-project.org/).

The program uses the [CGNS](https://cgns.github.io/) library version 3.1.4, for storing the calculation results.

The logo is heavily inspired by the software over at [suckless.org](https://suckless.org/).

![](doc/comparison.svg)

# YOU can help this project

Due to the Corona situation, I was separated from my more powerful tower PC at home and now only have my ThinkPad X220 with me. Some of the example cases just take for ever to complete, so I am looking for people that can test `ccfd` on these cases:
- [ ] `blasius`
- [ ] `cylinder`
- [ ] `doubleMachReflection`
- [ ] `forwardFacingStep`
- [ ] `naca0012`
- [ ] `naca2312`
- [ ] `pulse`
- [ ] `richtmyerMeshkov`
- [ ] `vortexStreet`

# Dependencies

- `git`
- `make`
- `cmake`
- `gcc`
- `gnuplot` (optional, for displaying calculation residuals, available [here](http://www.gnuplot.info/))
- `gmsh` (optional, for mesh generation, available [here](http://gmsh.info/))
- `paraview` (optional, for post-processing the results, available [here](https://www.paraview.org/))

# Installation

The installation process is easiest on Linux, but possible on MacOS and Windows.

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

Continue with [Usage under Linux](#usage-under-linux)

## MacOS

It should work basically the same as with Linux, but detailed instructions will follow soon...

## Windows

I was not really able to install it on windows directly, so you will have to use a [Linux Bash shell](https://docs.microsoft.com/en-us/windows/wsl/install-win10) for compiling and running the calculations. Get the latest Ubuntu shell. After the installation process is done, start the Ubuntu shell and install the necessary utilities
```
# apt update && apt upgrade
# apt install git make cmake gcc
```
Now go to your desired working directory. I would suggest selecting something in your actual Windows drive, as this will make accessing the calculation results a lot easier
```
$ cd /mnt/c/Users/<YourUserName>/Desktop
$ git clone https://github.com/hhh95/ccfd.git
$ cd ccfd
$ make
```
This should compile the code and create the two folders `obj` and `bin`, the last one containing the `ccfd` executable.

After everything is set up, continue with [Usage under Linux](#usage-under-linux). However, when installing ParaView, do not install it in the Ubuntu shell, but rather install it normally for [Windows](https://www.paraview.org/download/).

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
and observe the output. There should be four new files. The initial condition of the case and the calculation results at t = 0.25 s. They should all be `.csv` files. You can examine them with any spreadsheet program you like. Alternatively you can use `paraview`, a free post-processing program, that can visualize 1D, 2D, and 3D data. For Arch-based distributions you can install it from the package manager
```
# pacman -S paraview
```

On Ubuntu, the ParaView program in the repositories does not read CGNS files correctly for some reason. You will need to download ParaView 5.8 from the [ParaView website](https://www.paraview.org/download/). Next do the following
```
$ cd folder/where/you/downloaded/paraview
$ tar -xvf ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz
# mv ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz /opt
$ echo "export PATH:$PATH:/opt/ParaView-5.8.0-MPI-Linux-Python3.7-64bit/bin" >> ~/.bashrc
```
Next open the ParaView program in the directory where you performed the calculations
```
$ cd path/to/ccfd/calc/riemann
$ paraview &
```
Now, click on *File*->*Open* and select both sets of `.csv` files. Because we will be looking at 1D data, we need to change from *Render View* to *Line Chart View*. In the top right of the viewing area, click on the `X` button. Now select *Line Chart View* from the list. You should now see an empty grid with an x-, and a y-axis. Now in the *Pipeline Browser* to the left, click on the eye icons of the `.csv` files you have loaded in. A plot should appear in on the axis grid. It will probably show the initial state. In the top bar, click on the play button. Now, the final state should be shown. You will see the analytical, or exact, solution, as well as the numerical solution.

For more information on the theory, maybe have a look at [Wikipedia](https://en.wikipedia.org/wiki/Sod_shock_tube).

The procedure for running the other cases is the same. However, if the solution data is 2D, then you do not need to switch to *Line Chart View*. The 2D CGNS output files will usually have more than just the solution file. You can load everything at once by selecting the file that has `_Master` in its name. After loading the file, select all *Cell Arrays* in the *Pipeline Browser* and click on *Apply*. Then you can look at the different fields of the solution, by selecting them in the top bar (where it first says *Solid Color*).
