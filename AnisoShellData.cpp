#include "AnisoShellData.H"
#include <cmath>
#include "Utilities.H"
#include "LightParmParse/LParmParse.H"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Parallel.H"

const BaseCell* AnisoShellData::data_at(const double x[3], const double k[3], double &pathlength,DatasetStatus &status) const
{
    int tid = 0;
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #endif
    
    BaseCell *mycell = &cells[tid];
    
    double zdir[3] = {0.,0.,1.};
    const BaseCell* cell = SphericalShellData::data_at(x,k,pathlength,status);
    if(cell && (cell->density > 0. || cell->dust_density > 0.))
    {
        double r[3];
        for(int i=0;i<3;i++)
            r[i] = x[i]-0.5;
        double angle = acos(scalar(r,zdir)/sqrt(scalar(r,r)))*180./M_PI;
        if(angle >= 90.)
            angle -= 180.;
        
        if(fabs(angle)<opening_angle/2.)
        {
//            std::cout << angle << std::endl;
            mycell->density = 0.0;
            mycell->dust_density = 0.0;
        }
//         else
//         {
//             std::cout << "OUTSIDE angle is " << angle << std::endl;
//             
//         }
    }
    return cell;
        
}


int AnisoShellData::setup()
{
    int nthreads = 1;
    #ifdef _OPENMP
    nthreads = omp_get_max_threads();
    #endif
    if(Parallel::IOProcessor())
    {
        std::cout << "Do setup for AnisoShellData with " << nthreads << " threads." << std::endl;
    }

    cells = new BaseCell[nthreads];

    //read the geometry from the inputs file.
    params = new ShellParameters;
    LParmParse pp;
    double dummy;
    if(pp.query("anisoshell.column_density",dummy) && pp.query("anisoshell.density",dummy))
    {
        std::cout << "Please provide either density or column_density" << std::endl;
        std::abort();
        
    }
    if(pp.query("anisoshell.column_density_dust",dummy) && pp.query("anisoshell.density_dust",dummy))
    {
        std::cout << "Please provide either density_dust or column_density_dust" << std::endl;
        std::abort();
    }
    pp.get("anisoshell.inner_radius",params->inner_radius);
    pp.get("anisoshell.outer_radius",params->outer_radius);
    if(pp.query("anisoshell.column_density",dummy))
    pp.get("anisoshell.column_density",params->column_density);
    else
    pp.get("anisoshell.density",params->density);
    if(pp.query("anisoshell.column_density_dust",dummy))
    pp.get("anisoshell.column_density_dust",params->column_density_dust);
    else
    pp.get("anisoshell.density_dust",params->density_dust);

    pp.get("anisoshell.temperature",params->temperature);
    pp.get("anisoshell.outflow_velocity",params->outflow_velocity);
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
    //calculate the density in the anisoshell.
    if(pp.query("anisoshell.column_density",dummy))
    params->density = params->column_density/((params->outer_radius-params->inner_radius)*params->boxsize);
    else
    params->column_density = params->density*((params->outer_radius-params->inner_radius)*params->boxsize);
        
    if(pp.query("anisoshell.column_density_dust",dummy))
    params->density_dust = params->column_density_dust/((params->outer_radius-params->inner_radius)*params->boxsize);
    else
    params->column_density_dust = params->density_dust*((params->outer_radius-params->inner_radius)*params->boxsize);
    
    pp.get("anisoshell.opening_angle",opening_angle); 
    
    if(Parallel::IOProcessor())
    {
        std::cout << "Configured anisotropic shell with density " << params->density << " dust density " << params->density_dust <<  std::endl;
    }
    return 0;   
    
    


}
