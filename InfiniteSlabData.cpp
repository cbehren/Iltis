#include "InfiniteSlabData.H"
#include <cmath>
#include "Utilities.H"
#include "LightParmParse/LParmParse.H"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Parallel.H"

const BaseCell* InfiniteSlabData::data_at(const double x[3], const double k[3], double &pathlength,DatasetStatus &status) const
{
    int tid=0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if(x[0]<0.0 || x[0]>=1.0)
    {
        status = DS_NOT_IN_DOMAIN;
        return NULL;
    }
    BaseCell *mycell = &cells[tid];
    
    mycell->density = params->density;
    mycell->temperature = params->temperature;
    mycell->dust_density = params->density_dust;
    pathlength = intersection(x,k);
    
    pathlength = std::min(max_step,pathlength);
            
    if(pathlength!=pathlength)
    {
        std::cout << "Something went wrong; pathlength is NaN!" << std::endl;
        std::abort();
        
    }
    pathlength = std::max(1e-8,pathlength);
    
    return mycell;
}

int InfiniteSlabData::setup()
{
  int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    if(Parallel::IOProcessor())
    {
        std::cout << "Do setup for InfiniteSlabData with " << nthreads << " threads." << std::endl;
    }
   
 cells = new BaseCell[nthreads];
 
 //read the geometry from the inputs file.
 params = new SlabParameters;
 LParmParse pp;
 double dummy;
 if(pp.query("slab.column_density",dummy) && pp.query("slab.density",dummy))
 {
     std::cout << "Please provide either density or column_density" << std::endl;
     std::abort();
     
 }
 if(pp.query("slab.column_density_dust",dummy) && pp.query("slab.density_dust",dummy))
 {
     std::cout << "Please provide either density_dust or column_density_dust" << std::endl;
     std::abort();
     
 }

 if(pp.query("slab.column_density",dummy))
    pp.get("slab.column_density",params->column_density);
 else
    pp.get("slab.density",params->density);
 if(pp.query("slab.column_density_dust",dummy))
    pp.get("slab.column_density_dust",params->column_density_dust);
 else
    pp.get("slab.density_dust",params->density_dust);

 pp.get("slab.temperature",params->temperature);
 max_step = 1e-4;
 pp.query("max_step",max_step);
 pp.get("boxsize",params->boxsize); 
 

  if(params->column_density < 0.0 || params->column_density_dust < 0.0)
 {
     std::cout << "InfiniteSlabData::setup(): densities can not be negative!\n";
     std::abort();
 } 
   if(params->boxsize < 0.0)
 {
     std::cout << "InfiniteSlabData::setup(): boxsize can not be negative!\n";
     std::abort();
 } 
 //calculate the density in the slab.
  if(pp.query("slab.column_density",dummy))
    params->density = params->column_density/params->boxsize;
  else
    params->column_density = params->density*params->boxsize;
      
  if(pp.query("slab.column_density_dust",dummy))
    params->density_dust = params->column_density_dust/params->boxsize;
  else
    params->column_density_dust = params->density_dust*params->boxsize;
if(Parallel::IOProcessor())
{
    std::cout << "Configured spherical slab with density " << params->density << " dust density " << params->density_dust <<  std::endl;
    
}
 
 
 
 
 
 return 0;   
}

int InfiniteSlabData::set_parameters(SlabParameters sh)
{
    params = new SlabParameters;
    *params = sh; 
    return 0;
    
}



double InfiniteSlabData::intersection(const double x[3], const double k[3]) const
{
    //figure out the distance to the boundary of the slab
    if(k[0]==0.0)
        return 10000.0;
    
    double n[3] = {1,0,0};
    double p0[3] = {0,0,0};
    if(k[0]>0.0)
        p0[0]=1.0;
    double p0minusx[3];
    for(int i=0;i<3;i++)
        p0minusx[i] = p0[i]-x[i];
    double d = scalar(p0minusx,n)/scalar(k,n);
    return d;
    
    
}
