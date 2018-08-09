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
- `redshiftShifted`: TODO rename
- `hubbleFlow`: if set, determines the hubble rate in km/s/Mpc. TODO rename
- `no_hubble_flow`: set to true to disable the Hubble flow
- `max_num_peeling_off_photons`: maximum number of peeling off photons launched in one go (to cope with memory limitations)
- `max_step`: the maximum allowed spatial size of a step a photon makes in code units. Main purpose is to limit the space a photon travels between two evaluations of the Hubble flow the photon sees
- `use_biasing`: experimental feature to use biasing in dust absorption. 
- `minimum_bias`: sets the bias factor below which photons are disarded
- `number_of_instruments`: the total number of virtual instruments used for the peeling off algorithm
- `number_of_photons`: the total number of photons launched. Might be overwritten by `emission.minimum_number_of_photons_per_source`
- `split_domain`: TODO dont remember this one
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


#### TPG 
- `tpg.max_refinement_level`: maximum refinement level in TPG 
- `tpg.refine_on_density`: density threshold we need to reach in order to refine a cell
- `tpg.rootlevel`: the root level of the initial grid, with 2^rootlevel being the linear size of grid

#### Unigrid
- `unigrid.filename`: filename for reading in the grid
- `unigrid.output_tau_stats`: set to true to generate statistics on the optical depths in the grid


#### Octet
- `octet.distribution_method`: 
- `octet.distribution_method.row_dir`: 
- `octet.name`: 
- `octet.root`: 
- `octet.rootlevel`: 

#### Ramses
- `ramses.distribution_method`: 
- `ramses.distribution_method.row_dir`: 
- `ramses.ionize_cells`: 
- `ramses.ionize_cells.remove_dust`: 
- `ramses.output_tau_stats`: 
- `ramses.root`: 
- `ramses.scale_density`: 
- `ramses.scale_dust_density`: 
- `ramses.scale_velocity`: 
- `ramses.star_file`: 

#### Plotter

- `plotter.depth`: 
- `plotter.npixels`: 
- `plotter.opticalDepth.wavelength`: 
- `plotter.opticalDepthGrid.divide_by_density`: 
- `plotter.opticalDepthGrid.lambda_left`: 
- `plotter.opticalDepthGrid.lambda_right`: 
- `plotter.opticalDepthGrid.nlambda`: 
- `plotter.output_prefix`: 
- `plotter.oversampling`: 
- `plotter.width`: 
