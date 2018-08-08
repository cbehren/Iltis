# Iltis,

a light-weight line emission transfer code

## Requirements

You will need
- a C++ compiler ready for C++11 features
- MPI libraries/headers
- The BaseOctet template class (resides in a private repo on Bitbucket), cloned into ../BaseOctet if . is the Iltis directory
- for the python tools, you need an installation of Python 2
- for the text2hdf5 files, you need the hdf5 library (also in the C++ flavor!), and probably you need to adapt the pathes in text2hdf5/Makefile 

## Installation

Set CC to your preferred c++ compiler in Makefile and just do make!

## Usage

Run Iltis using

./Iltis.exe \<inputs file\>
  
For testing purposes, you can use the inputs file in RegressionTests/Shell (should also work on a Desktop

## Known Issues

The strict MPI requirement comes from the fact that the BaseOctet template used is not able to build without MPI. I need to fix that at some point.
