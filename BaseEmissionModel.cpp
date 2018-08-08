#include "BaseEmissionModel.H"
#include "LightParmParse/LParmParse.H"
#include "RandomNumbers.H"
#include "Parallel.H"
#include "LymanAlphaLine.H"
#include "Utilities.H"
#include <string>

void BaseEmissionModel::setup()
{
    
    
}


void BaseEmissionModel::launch_bunch(BaseParticleVector& particles, const BaseEmissionLine& line, const BaseDataset& ds, double dt)
{
    //note that this is the implementation for the base class; it is basically a stub routine putting photons into the center of the box.
    if(all_done)
        return;
    LParmParse pp;
    long number_of_photons;
    pp.get("number_of_photons",number_of_photons);
    double center[3] = {0.5, 0.5, 0.5};
    std::vector<double> temp_center;
    temp_center.resize(3,-1.0);
    pp.queryarr("emission.center",temp_center,LParmParse::FIRST,3);
    if(temp_center[0]>=0.0)
    {
        if(Parallel::IOProcessor())
            std::cout << "BaseEmissionModel::setup(): Set center of emission to " << temp_center[0] << " " << temp_center[1] << " " << temp_center[2] << std::endl;
        for(int i=0;i<3;i++)
            center[i] = temp_center[i];
    }
    
    int nprocs = Parallel::NProcs();
    //get the amount of particles we have to launch.
    int photons_per_process = number_of_photons/nprocs;
    if(Parallel::IOProcessor())
    {
        if(number_of_photons<nprocs)
            photons_per_process = number_of_photons;
        else
            photons_per_process += number_of_photons % nprocs;
    }
    //get the emission spectrum
    std::string str("delta");
    pp.query("emission.spectrum",str);
    
    //initialize every photon
    for(long i=0;i<photons_per_process;i++)
    {
        BaseParticle p;
        for(int j=0;j<3;j++)
        {
            p.x[j] = center[j];
             
        }
        RNG::point_on_sphere(p.k);
        p.optical_depth = RNG::exponential();
        //we use -1 as a flag here, see BaseSimulation.cpp for details.
        p.number_of_scatterings = -1;
        if(str.compare(std::string("delta"))==0)
        {            
            p.frequency = 0.0;
        }
        else if(str.compare(std::string("gaussian"))==0)
        {
            //in case we want a gaussian, figure out the width of the gaussian.
            double width;
            pp.get("emission.fixed_width",width);
            if(Parallel::IOProcessor() && i==0)
                std::cout << "Using gaussian input spectrum of width " << width;
            double lineshift_in_kms = RNG::gaussian(0,width);
            p.frequency = line.shift_to_frequency(lineshift_in_kms);            
        }
        else
        {
            std::cout << " BaseEmissionModel::setup(): Unsupported model for spectrum!" << std::endl;
            std::abort();
            
        }
        particles.push_back(p);
    }
    
    long size = particles.size();
    Parallel::ReduceLongSum(size);
    if(Parallel::IOProcessor())
        std::cout << "Launched "<< size << " photons" << std::endl;
    //in case the data is distributed, we need to exchange particles so that each process has the right photons. We need to do that before we call shift_into_observed_frame, because that routine requires access to the hydro cell in which each photons lives.
    ds.set_domain(&particles);
    particles.exchange();
    for(unsigned long i=0;i<particles.size();i++)
        shift_into_observed_frame(particles[i],line,ds);
    all_done = true;
}

void BaseEmissionModel::init_point_source(BaseParticleVector& particles,const long nphotons,const Emitter& emitter,const double weight,const BaseEmissionLine& line)
{
    
    double center[3];
    for (int i = 0; i < 3; i++)
       center[i] = emitter.position[i];
    for (long icnt = 0; icnt < nphotons; icnt++)
    {
        BaseParticle p;
        for(int i=0;i<3;i++)
        {
            p.x[i] = center[i];
            p.lsp[i] = center[i];
        }
        p.weight = weight;
        RNG::point_on_sphere(p.k);       
        p.optical_depth = RNG::exponential();
        p.number_of_scatterings = -1;
        p.emitter = emitter.index;
        if(emitter.width>0.0)
        {
            double lineshift_in_kms = RNG::gaussian(0,emitter.width);
            p.frequency = line.shift_to_frequency(lineshift_in_kms);
        }
        particles.push_back(p);
       
        
    }
    //CAUTION remember that after calling this routine, the photons are not yet shifted into the fluid cells' frame! 
 
    
}

void BaseEmissionModel::shift_into_observed_frame(BaseParticle &p,const BaseEmissionLine& line, const BaseDataset& ds)
{
    double pathlength=0.0;
    DatasetStatus status;
    const BaseCell* cell = ds.data_at(p.x,p.k,pathlength,status);
    if(status!=DS_ALL_GOOD)
    {
        std::cout << "BaseEmissionModel::shift_into_observed_frame(): Could not access data for particle" << std::endl;
        std::abort();
    }
    const double* velocity = cell->velocity;
    p.frequency = p.frequency + scalar(p.k,velocity)*line.Nu0()/BaseEmissionLine::c_val;
    p.status = BP_OKAY;
}
