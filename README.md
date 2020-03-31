![](doc/ccfd.svg)

A two-dimensional finite volume computational fluid dynamics code, written in C

This code is intended to become a drop-in replacement for `cfdfv`, a CFD code written in Fortran by the [Institute of Aerodynamics and Gas Dynamics](http://www.iag.uni-stuttgart.de) at the University of Stuttgart for a CFD programming course. This code itself is not available online, as far as I know, but it is a modified version of [FLEXI](https://www.flexi-project.org/).

The program uses the [CGNS](https://cgns.github.io/) library version 4.1.1, for storing the calculation results.

# Installation

## Dependencies

- `git`
- `make`
- `gcc`
- `gmsh` (optional, for mesh generation)
- `paraview` (optional, for post-processing the results)
- `gnuplot` (optional, for displaying calculation residuals)

## Linux

### Build

First make sure that all necessary dependencies are all installed. These can usually be obtained through your distributions package manager, on Arch based systems the following command should suffice
```
# pacman -S git base-devel
```
For an Ubuntu like system the following command should do the trick
```
# apt-get install git build-essential
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

## MacOS

It should work basically the same as with Linux, but detailed instructions will follow soon...

## Windows

It should work basically the same as with Linux, but detailed instructions will follow soon...

# Usage

Navigate to the `calc` folder. Here you will find example input files, such as the SOD test case. The calculation can be started by simply executing `ccfd` followed by the `.ini` file, containing the case setup
```
$ ../bin/ccfd sod.ini
```
The program will display information about the initialization and time stepping process until the calculation is finished. Depending on what output format was set in the `.ini` file, the results at the specified output times/intervals are written to disk.
