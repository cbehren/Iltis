#ifndef __BASE_EMISSION_MODEL_
#define __BASE_EMISSION_MODEL_
#include "BaseParticleVector.H"
#include "EmissionList.H"
#include "BaseDataset.H"
#include "BaseEmissionLine.H"

//defines the base class for how and where photons are initially emitted. The work is primiarly done in the virtual setup routine.


class BaseEmissionModel{
public:
    //setup the emitter(s)
    virtual void setup(BaseParticleVector& particles,const BaseEmissionLine& line, const BaseDataset& ds);
    //helper function to init a single point source with nphotons photons. CAUTION: You need to shift the photons into the fluid frame afterwards!
    void init_point_source(BaseParticleVector& particles,const long nphotons,const Emitter& emitter,const double weight,const BaseEmissionLine& line);
    //helper function to shift the photons into their fluid frame initially. 
    void shift_into_observed_frame(BaseParticle &p,const BaseEmissionLine& line, const BaseDataset& ds);
private:
    
    
};
#endif
