#include "BaseDataset.H"
#include "BaseEmissionLine.H"
#include "BaseParticleVector.H"
#include "BaseEmissionModel.H"
#include "PlottingOperators.H"
#include "BaseOutput.H"
#include "Utilities.H"

#include <string>
#include <vector> 

//here is where the magic happens. Defines the main simulation object, running the main loop(s).



class PlottingInterface{
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
    //initialize photons
    int generate_rays();
    //propagate photons
    int do_rays();
    //exchange photons if necessary, check for escaped/absorbed photons, write output.
    int do_communication();
    //number of directions in which we keep track of peeling off photons
    int number_of_instruments = 0;
    //holds these directions
    std::vector<line_of_sight> instruments;
    //initialize a new peeling off photon given the original particle, the cell, the atom velocity from the relevant scattering.
    std::string dataset_type;
    //vector of particles
    BaseParticleVector particles;
    //emission line object
    BaseEmissionLine* line = NULL;
    //dataset object
    BaseDataset* ds = NULL;
    //output module
    BaseOutput* o = NULL;
    //emission model object
    BaseEmissionModel* e = NULL;
    //plotting parameters
    line_of_sight los;
    line_of_sight center;
    line_of_sight up;
    double depth;
    double width;
    int npixels;
    int oversampling;
    std::string output_prefix;
    
    std::vector<BasePlotOperator*> plotops;

    int verbosity = 1;
};
