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


## A quick documentation

### Introduction

Iltis is a Monte-Carlo radiative transfer code specialized on the Lyman-$\alpha$ line. However, it is relatively simple to change the code to run with a different emission line. At startup, all necessary parameters for a run a read from parameter file that is given as the only argument to the executable. Besides the parameters in those files, there are a few preprocessor directives in use, mainly for the acceleration schemes implemented to speed up the RT.

The code follows Zheng et al. 2002, Laursen et al. 2009, Dijkstra et al. 2006 in terms of the RT scheme. Detailed references are given in the code.

Iltis is relatively agnostic about the way the geometry of the problem, that is, density, velocity, temperature, .., fields are configured. It is therefore not very complicated to add own datasources. In this release, only the SphericalShell, the InfiniteSlab, and the Unigrid datasource/problem type are included. An implementation for reading and using Ramses files is available on request.

### Problem/Dataset types

#### InfiniteSlab

The infinite slab is a frequently used test problem. It considers a slab of gas, infinite in two dimensions, and finite in one, with a length 2 times `boxsize`. The problem is completely specified by (colum) density, dust (column) density, and temperature. 

#### SphericalShell

The spherical shell is a well-known, simple model of LAEs pioneered by Verhamme et al. 2006. It considers a spherical shell with an inner and outer radius, a gas/dust (column) density, a temperature, and an constant outflow velocity.

#### Unigrid   

A single, fixed grid read from disk, specifying the density, dust density, temperature, and velocity field in each cell. The format of the grid file is very simple and can be understood from the example in `RegressionTests/Unigrid/testgrid.py`: In short, the first line of the file should contain a hash, followed by the linear grid size. All other lines contain the data for one cell, in this order: density, temperature, velocity, dust density, each separated by a whitespace.




### Parameters
#### general parameters
- `boxsize`: the size of the domain in cm. 
- `dataset_type`: type of dataset used
- `seed`: the initial random seed of the simulation 
- `use_peeling_off`: set to true to enable peeling off
- `verbosity`: set to value > 1 to get more information about the progress of a simulation
- `cosmology.H0`: Hubble constant at current time. Used for calculation of Hubble rate.
- `cosmology.Omega_L`: Omega Lambda
- `cosmology.Omega_M`: Omega Matter
- `redshift`: redshift of the simulation, used to calculate the Hubble flow rate if not overriden by hubble_flow
- `redshift_shifted`: deprecated
- `hubble_flow`: if set, determines the hubble rate in km/s/Mpc. 
- `no_hubble_flow`: set to true to disable the Hubble flow
- `max_num_peeling_off_photons`: maximum number of peeling off photons launched in one go (to cope with memory limitations)
- `max_step`: the maximum allowed spatial size of a step a photon makes in code units. Main purpose is to limit the space a photon travels between two evaluations of the Hubble flow the photon sees
- `use_biasing`: experimental feature to use biasing in dust absorption. 
- `minimum_bias`: sets the bias factor below which photons are disarded
- `number_of_instruments`: the total number of virtual instruments used for the peeling off algorithm
- `number_of_photons`: the total number of photons launched. Might be overwritten by `emission.minimum_number_of_photons_per_source`
- `tau_max`: the optical depth above which we discard peeling off photons


#### dust modelling
- `dust.albedo`: the albedo of dust in the simulation
- `dust.dahlia.type`: use MW/LMC/SMC type dust in the Dahlia model
- `dust.dust_metal_ratio`: ratio of dust to metals in the simulation
- `dust.model`: the model used for dust in the simulation

#### emission
- `emission.bounding_box_hi`: if set, defines the upper bound of region in which photons can be emitted (applicable only for the EmissionList module)
- `emission.bounding_box_lo`: if set, defines the lower bound of region in which photons can be emitted (applicable only for the EmissionList module)
- `emission.center`: the point in code coordinates where emission of photons takes place (applicable only for the BaseEmission module)
- `emission.emitter_file`: if set, the EmissionList module is used to read a list of point emitters
- `emission.fixed_width`: sets the width of the intrinsic spectrum of photons (if the spectrum is Gaussian)
- `emission.minimum_luminosity`: minimum luminosity of a source to be considered, in code units (1e42 erg/s) (applicable only for the EmissionList module)
- `emission.minimum_number_of_photons_per_source`: the minimum number of photons launched per source. Photon weights are scaled in order to ensure energy conservation
- `emission.spectrum`: type of the intrinsic spectrum (delta or Gaussian)

#### neutral fraction
- `neutral_fraction.model`: the type of model used for calculating the neutral fraction of gas



#### Output
- `output.binary`: set to true to generate output directly as binary
- `output.do_slices`: set to true to write a centered slice of the simulation domain
- `output.peeling_off_prefix`: prefix of the outputfiles for the peeling off
- `output.prefix`: prefix of the output files for the MC photons
- `output.slices.resolution`: resolution of slice

#### SphericalShell
- `shell.column_density`: column density of neutral hydrogen. Either this or shell.density must be set. Unit is 1/cm^2
- `shell.column_density_dust`: column density of dust. Either this or shell.dust_density must be set. Unit is g/cm^2
- `shell.density`: see above
- `shell.density_dust`: see above
- `shell.inner_radius`: inner radius of the shell in code units
- `shell.outer_radius`: outer radius of the shell in code units
- `shell.outflow_velocity`: the (constant) outflow velocity of the shell in km/s
- `shell.temperature`: temperature of the gas in the shell


#### InfiniteSlab
- `slab.column_density`: column density of neutral hydrogen. Either this or slab.density must be set. Unit is 1/cm^2
- `slab.column_density_dust`: column density of dust. Either this or slab.dust_density must be set. Unit is g/cm^2
- `slab.density`: see above
- `slab.density_dust`: see above
- `slab.temperature`: temperature of the gas in the slab


#### Unigrid
- `unigrid.filename`: filename for reading in the grid
- `unigrid.output_tau_stats`: set to true to generate statistics on the optical depths in the grid


### preprocessor parameters
#### defined in `LymanAlphaLine.H`
`USE_ACCELERATION`: use an acceleration scheme based on truncating the scattering atoms' velocitiy to avoid core scatterings
`FULL_ACC`: truncate at a fixed value of the dimensionless frequency, x, in the atom's frame. If not defined, use a scheme by Laursen et al. 2008 to dynamically determine the cutoff value
`ÀCC_XCRIT`: the fixed value at which we truncate in case ``FULL_ACC`` is defined
`DIPOLE_SCATTERING`: instead of an isotropic phase function, use a dipole
`X_CRIT_DIPOLE`: use dipole scattering above this value (dimensionless frequency)
#### defined in `BaseSimulation.H`
`NEUFELD_ACCELERATION_SCHEME`: use a semi-analytic scheme to exit cells with extremely high optical depths, see Laursen et al. 2008 for details

