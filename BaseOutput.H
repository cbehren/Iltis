#ifndef __BASE_OUTPUT_
#define __BASE_OUTPUT_

#include "BaseParticle.H"
#include "BaseParticleVector.H"
#include "BaseEmissionLine.H"
#include "BaseDataset.H"
#include <fstream>

//defines the default output module.
class BaseOutput{
public:
    //setup the module
    virtual void setup(BaseEmissionLine* l);
    //write out the initial properties of emitted photons
    virtual void writeInputPhotons(BaseParticleVector &vp);
    //write out properties of photons that escaped
    virtual void writePhoton(BaseParticle &p);
    //write out properties of photons that have been absorbed by dust
    virtual void writeDeadPhoton(BaseParticle &p);
    //write out properties of peeling off photons
    virtual void writePeelingOff(BaseParticle &p,int index);
    //write out a slice through the hydro data (for debugging purposes)
    virtual void writeSlices(const BaseDataset& ds, std::string snap_time);
    //same, but for distributed data. TODO should be private
    virtual void writeSlicesParallel(const BaseDataset& ds);
    //close and combine all relevant output files.
    virtual void cleanup();
    
private:
    //open output files
    void open_files();
    //close output files
    void close_files();
    //holds the files for the peeling off writes
    std::ofstream* pfiles=NULL;
    //holds the file for the escaped photons
    std::ofstream file;
    //same, but for absorbed photons
    std::ofstream file_dead;
    //and here we have the filenames.
    std::string fname;
    std::string fname_dead;
    std::string fname_peeling;
    //true if we want to do output in a binary fashion. saves some spaces.
    bool binary = false;
    //the number of files used.
    int nfiles = 0;
    //the emission line object
    BaseEmissionLine* line = NULL;
    //total number of scatterings
    long totalScatterings = 0;
    //total number of escaped photons
    long totalDone = 0;
    //total number of absorbed photons
    long totalDeadDone = 0;
    //total number of peeling off contributions.
    long contributions = 0;
    
};
#endif
