#ifndef __SLIM_OUTPUT_
#define __SLIM_OUTPUT_

#include "BaseParticle.H"
#include "BaseParticleVector.H"
#include "BaseEmissionLine.H"
#include "BaseDataset.H"
#include "BaseOutput.H"
#include <fstream>
#include "Histogram.H"

//defines the default output module.
class SlimOutput : public BaseOutput{
public:
    //TODO BaseParticles can be const arguments.
    //write out the initial properties of emitted photons
    virtual void writeInputPhotons(BaseParticleVector &vp);
    //write out properties of photons that escaped
    virtual void writePhoton(BaseParticle &p);
    //write out properties of photons that have been absorbed by dust
    virtual void writeDeadPhoton(BaseParticle &p);
    //write out properties of peeling off photons
    virtual void writePeelingOff(BaseParticle &p,int index);
    //close and combine all relevant output files.
    virtual void cleanup();
    virtual ~SlimOutput();
    
protected:
    //open output files
    virtual void open_files();
    virtual void open_peeling_files();
    //close output files
    virtual void close_files();
    
    Histogram* escaped_spectrum=NULL;
    Histogram* destroyed_spectrum=NULL;
    Histogram** peeling_spectra=NULL;
    
    
};
#endif
