#include "SphericalShellData.H"
#include <cmath>
#include "Utilities.H"
#include "LightParmParse/LParmParse.H"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Parallel.H"

const BaseCell* SphericalShellData::data_at(const double x[3], const double k[3], double &pathlength,DatasetStatus &status) const
{
    int tid=0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    if(x[0]<0.0 || x[1]<0.0 || x[2]<0.0 || x[0]>=1.0 || x[1]>=1.0 || x[2]>=1.0)
    {
        status = DS_NOT_IN_DOMAIN;
        return NULL;
    }
    double radius=0.0;
    double r[3] = {0,0,0};
    double p = -1;
    for(int i=0;i<3;i++)
    {
        r[i] = x[i]-0.5;
        radius += r[i]*r[i];
    }
    status = DS_ALL_GOOD;
    radius = sqrt(radius);
    BaseCell *mycell = &cells[tid];
    if(radius >= params->inner_radius && radius < params->outer_radius)
    { 
        mycell->density = params->density;
        mycell->temperature = params->temperature;
        for(int i=0;i<3;i++)
            if(radius>0.0)
                mycell->velocity[i] = params->outflow_velocity*r[i]/radius;
            else
                mycell->velocity[i] = 0.0;
        mycell->dust_density = params->density_dust;
        double p1 = intersection(x,k,params->inner_radius);
        double p2 = intersection(x,k,params->outer_radius);
        if(params->inner_radius == 0.0)
            p1 = -1.0;
        if(p1 >= 0.0 && p2 >= 0.0)
            pathlength = std::min(p1,p2);
        else
            pathlength = std::max(p1,p2);
//         std::cout << "Comparing " << p1 << " and " << p2 << std::endl;
        
    }
    else
    {
        
        mycell->density = 0.0;
        mycell->dust_density = 0.0;
        for(int i=0;i<3;i++)
            mycell->velocity[i] = 0.0;
        mycell->temperature = params->temperature;
        if(radius > params->outer_radius)
            pathlength = 1.0;
        else
        {
            if(radius < params->inner_radius)
            {
                p = intersection(x,k,params->inner_radius);
//                 std::cout << "Distance to inner radius: " << p << std::endl;
            }
            else
            {
                p = intersection(x,k,params->outer_radius);
//                std::cout << "Distance to outer radius: " << p << std::endl;
            }
        }
        if(p<=0.0)
            pathlength = 1.0;
        else
            pathlength = p;
        
    }
    
            
    if(pathlength!=pathlength)
    {
        std::cout << "Something went wrong; pathlength is NaN!" << std::endl;
        std::abort();
        
    }
    pathlength = std::min(max_step,pathlength);
    pathlength = std::max(1e-8,pathlength);
    
    return mycell;
}

int SphericalShellData::setup()
{
  int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    if(Parallel::IOProcessor())
    {
        std::cout << "Do setup for SphericalShellData with " << nthreads << " threads." << std::endl;
    }
   
 cells = new BaseCell[nthreads];
 
 //read the geometry from the inputs file.
 params = new ShellParameters;
 LParmParse pp;
 double dummy;
 if(pp.query("shell.column_density",dummy) && pp.query("shell.density",dummy))
 {
     std::cout << "Please provide either density or column_density" << std::endl;
     std::abort();
     
 }
 if(pp.query("shell.column_density_dust",dummy) && pp.query("shell.density_dust",dummy))
 {
     std::cout << "Please provide either density_dust or column_density_dust" << std::endl;
     std::abort();
     
 }
 pp.get("shell.inner_radius",params->inner_radius);
 pp.get("shell.outer_radius",params->outer_radius);
 if(pp.query("shell.column_density",dummy))
    pp.get("shell.column_density",params->column_density);
 else
    pp.get("shell.density",params->density);
 if(pp.query("shell.column_density_dust",dummy))
    pp.get("shell.column_density_dust",params->column_density_dust);
 else
    pp.get("shell.density_dust",params->density_dust);

 pp.get("shell.temperature",params->temperature);
 pp.get("shell.outflow_velocity",params->outflow_velocity);
 max_step = 1e-4;
 pp.query("max_step",max_step);
 pp.get("boxsize",params->boxsize); 
 
 //some checks
 if(params->inner_radius >= params->outer_radius)
 {
     std::cout << "SphericalShellData::setup(): inner radius can not be smaller than outer!\n";
     std::abort();
 }  
  if(params->column_density < 0.0 || params->column_density_dust < 0.0)
 {
     std::cout << "SphericalShellData::setup(): densities can not be negative!\n";
     std::abort();
 } 
   if(params->boxsize < 0.0)
 {
     std::cout << "SphericalShellData::setup(): boxsize can not be negative!\n";
     std::abort();
 } 
 //calculate the density in the shell.
  if(pp.query("shell.column_density",dummy))
    params->density = params->column_density/((params->outer_radius-params->inner_radius)*params->boxsize);
  else
    params->column_density = params->density*((params->outer_radius-params->inner_radius)*params->boxsize);
      
  if(pp.query("shell.column_density_dust",dummy))
    params->density_dust = params->column_density_dust/((params->outer_radius-params->inner_radius)*params->boxsize);
  else
    params->column_density_dust = params->density_dust*((params->outer_radius-params->inner_radius)*params->boxsize);
if(Parallel::IOProcessor())
{
    std::cout << "Configured spherical shell with density " << params->density << " dust density " << params->density_dust <<  std::endl;
    
}
 
 
 
 
 
 return 0;   
}

int SphericalShellData::set_parameters(ShellParameters sh)
{
    params = new ShellParameters;
    *params = sh; 
    return 0;
    
}



double SphericalShellData::intersection(const double x[3], const double k[3],double radius) const
{
    //figure out the distance to the sphere to speed up.
    
    //for reference, see wikipedia: https://en.wikipedia.org/wiki/Line–sphere_intersection
    
    double ominusc[3];
    for(int i=0;i<3;i++)
        ominusc[i] = x[i] - 0.5;
    double c1 = pow(scalar(k,ominusc),2.0);
    double c2 = norm(ominusc)*norm(ominusc);
    double c3 = radius*radius;
    
    double csquared = c1 - c2 + c3;
    if(csquared<=0.0)//we do exclude the tangent case here.
    {
        return -1.0;
    }
    else
    {
        double solution1 = - scalar(k,ominusc) + sqrt(csquared);
        double solution2 = - scalar(k,ominusc) - sqrt(csquared);
        if(solution1<0.0 && solution2 < 0.0)
            return -1;
        if(solution1>0.0 && solution2>0.0)
            return std::min(solution1,solution2);
        else
            return std::max(solution1,solution2);
    }
    
    
}
