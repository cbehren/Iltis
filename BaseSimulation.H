#include "BaseDataset.H"
#include "BaseEmissionLine.H"
#include "BaseParticleVector.H"
#include "BaseEmissionModel.H"
#include "BaseOutput.H"
#include "Utilities.H"
#include <string>
#include <vector> 

//here is where the magic happens. Defines the main simulation object, running the main loop(s).


class BaseSimulation{
public:
    //setup the simulation
    virtual int setup();
    //to a single simulation step
    int do_step();
    //return true if simulation is done
    bool is_done();
    //cleanup after finishing simulation
    virtual int cleanup();
    //read the relevant parameters from the parameter file
    int read_params();

protected:
    //propagate peeling off photons
    int do_peeling_off();
    //propagate photons
    int do_rays();
    //exchange photons if necessary, check for escaped/absorbed photons, write output.
    int do_communication();
    //true if we use peeling off
    bool use_peeling_off=false;
    //maximum optical depth allowed for a peeling off photon before it is discarded
    double tau_max = 20.0;
    //number of directions in which we keep track of peeling off photons
    int number_of_instruments = 0;
    //holds these directions
    std::vector<line_of_sight> instruments;
    //initialize a new peeling off photon given the original particle, the cell, the atom velocity from the relevant scattering.
    void launch_peeling(const BaseParticle &p, const BaseCell& cell, double atom_velocity[3],bool on_dust=false);
    //maximum number of peeling off photons to be launched before we stop the current step and cycle through the peeling off photons. necessary to reduce memory usage.
    long max_num_peeling_off_photons=1e9;
    //true if we use a file specifying our emitters
    bool use_emitter_file = false;
    //true if use biasing
    bool use_biasing = false;
    //bias below which we discard a photon
    double minimum_bias = 0.1;
    //string indicating which dataset type we use
    std::string dataset_type;
    //vector of particles
    BaseParticleVector particles;
    //array of vectors of peeling off particles, one for each instrument.
    BaseParticleVector  *vparticles = NULL;
    //emission line object
    BaseEmissionLine* line = NULL;
    //dataset object
    BaseDataset* ds = NULL;
    //output module
    BaseOutput* o = NULL;
    //emission model object
    BaseEmissionModel* e = NULL;
    
    int verbosity = 1;
#ifdef REMOVE_CGM              
    //TODO remove
    const double center[3] = {0.53033, 0.51041, 0.50396};
    const double radius = 5e-2;
#endif
    
    
};
