# Iltis,

a light-weight line emission transfer code

## Requirements

You will need at least
- a C++ compiler with C++11 features
- OpenMP support is required, although you can just remove the -fopenmp option in the `Makefile` to get around this

Additionally, it is useful to have
- MPI libraries/headers to run in MPI parallel mode
- for the Python tools, you need an installation of Python 2
- for the tools that convert the output to HDF5 (in `./text2hdf5/`), you need the HDF5 library (also in the C++ flavor!), and probably you need to adapt the pathes in `text2hdf5/Makefile` 

## Installation

Set CC to your preferred C++ compiler in `Makefile` and just do `make`

## Usage

Run Iltis using

`./Iltis.exe \<inputs file\>`
  
For testing purposes, you can use the inputs file in `RegressionTests/Shell` (should also work on a Desktop)

## Known Issues

--


## A quick documentation

### Introduction

Iltis is a Monte-Carlo radiative transfer code specialized on the Lyman-$\alpha$ line. However, it is relatively simple to change the code to run with a different emission line. At startup, all necessary parameters for a run a read from parameter file that is given as the only argument to the executable. Besides the parameters in those files, there are a few preprocessor directives in use, mainly for the acceleration schemes implemented to speed up the RT.

The code follows Zheng et al. 2002, Laursen et al. 2009, Dijkstra et al. 2006 in terms of the RT scheme. Detailed references are given in the code.

Iltis is relatively agnostic about the way the geometry of the problem, that is, density, velocity, temperature, .., fields are configured. It is therefore not very complicated to add own datasources. In this release, only the SphericalShell, the InfiniteSlab, and the Unigrid datasource/problem type are included. An implementation for reading and using Ramses files is available on request.

### Problem types

#### The InfiniteSlab

The infinite slab is a frequently used test problem. It considers a slab of gas, infinite in two dimensions, and finite in one, with a length 2 times `boxsize`. Together with the column_density (or density), the column_density_dust (or dust_density), and the temperature, this defines the problem.  

