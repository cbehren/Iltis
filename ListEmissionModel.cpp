#include "ListEmissionModel.H"
#include "LightParmParse/LParmParse.H"
#include "RandomNumbers.H"
#include "Parallel.H"
#include "EmissionList.H"
#include <string>
#include <vector>

void ListEmissionModel::setup(BaseParticleVector& particles,const BaseEmissionLine& line, const BaseDataset& ds)
{
    const int       MyProc   = Parallel::MyProc();
    const int       NProcs   = Parallel::NProcs();
    LParmParse pp;
    std::string filename;
    int minimum_number_of_photons_per_source=0;
    double minimum_luminosity=0.0;
    long number_of_photons;
    pp.get("number_of_photons",number_of_photons);
    pp.get("emission.emitter_file",filename);
    pp.query("emission.minimum_number_of_photons_per_source",minimum_number_of_photons_per_source);
    pp.query("emission.minimum_luminosity",minimum_luminosity);
    if(Parallel::IOProcessor())
    {
        std::cout << "Generating emitters from point sources in file " << filename << std::endl;
        std::cout << "Each source will emit at least " << minimum_number_of_photons_per_source << " photons " << std::endl;
    }
    EmissionList emitters;
    std::vector<double> bounding_box_lo;
    std::vector<double> bounding_box_hi;
    double applyBoundingBox=false;
    pp.queryarr("emission.bounding_box_lo",bounding_box_lo);
    pp.queryarr("emission.bounding_box_hi",bounding_box_hi);
    if(bounding_box_lo.size()==3 && bounding_box_hi.size()==3) 
    {
      applyBoundingBox=true;
      emitters.setUpBoundingBox(bounding_box_lo[0],bounding_box_lo[1],bounding_box_lo[2],bounding_box_hi[0],bounding_box_hi[1],bounding_box_hi[2]);
      if(Parallel::IOProcessor())
        std::cout << "Using Bounding Box for Emission!" << std::endl;
    }
    EmissionModelEnum em = EM_DEFAULT;
    SpectrumModelEnum sp = SP_DELTA;
    std::string modelstring("delta");
    pp.query("emission.spectrum",modelstring);
    if(modelstring.compare(std::string("gaussian"))==0)
    {
      sp = SP_GAUSSIAN_FIXED;
      double fixed_width = -1.0;
      pp.get("emission.fixed_width",fixed_width);
      emitters.setFixedWidth(fixed_width);
      if(Parallel::IOProcessor())
        std::cout << "Using Fixed Gaussian Width of Emission Line of " << fixed_width << " km/s\n";
      
    }
    else if(modelstring.compare(std::string("ZZ10"))==0)
    {
        sp = SP_ZZ10;
        if(Parallel::IOProcessor())
            std::cout << "Using ZZ10 model for spectrum\n";
        
    }
    else if(!(modelstring.compare(std::string("delta"))==0))
    {
        std::cout << "unknown model for spectrum" << std::endl;
        std::abort();
        
    }
    emitters.readFromFile(filename,em,sp,minimum_luminosity,applyBoundingBox);
    if(Parallel::IOProcessor())
        std::cout << "Found " << emitters.size() << " emitters in file\n";
    std::vector<double> relative_emissivity = emitters.relativeEmissivity();
    std::vector<double> absolute_emissivity = emitters.absoluteEmissivity();
    for(unsigned int i=0;i<emitters.size();i++)
    {
        if ((i%NProcs) == (unsigned int)MyProc)
        {
            Emitter &emitter = emitters[i];
            int count = round(relative_emissivity[i]*number_of_photons);
            if(count<minimum_number_of_photons_per_source)
            {
                count = minimum_number_of_photons_per_source;
                
            }
//             std::cout << "Process " << MyProc << " initializing source " << i << " with " << count << " photons" << std::endl;
            double weight = absolute_emissivity[i]/count;
            
            init_point_source(particles,count,emitter,weight,line);
        }
    }
    ds.set_domain(&particles);
    particles.exchange();
    for(unsigned long i=0;i<particles.size();i++)
        shift_into_observed_frame(particles[i],line,ds);
    long size = particles.size();
    Parallel::ReduceLongSum(size);
    if(Parallel::IOProcessor())
        std::cout << "Launched "<< size << " photons" << std::endl;
    
   
}