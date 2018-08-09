#include "BaseEmissionLine.H"
#include "BaseCell.H"
#include "BaseParticle.H"
#include <cmath>
#include <iostream>
#include "LightParmParse/LParmParse.H"
#include "Parallel.H"
#include "Selectors.H"
#include <string>

double BaseEmissionLine::boxsize = 0;
double BaseEmissionLine::hubble_flow_value = 0;
const double BaseEmissionLine::c_val = 299792458.0*100;//cm/s
DustModule  *BaseEmissionLine::dm = NULL;


BaseEmissionLine::BaseEmissionLine()
{
    if(dm==NULL)
    {
        //initialize the dust module.
        LParmParse pp;
        std::string dustModelString("Null");
        pp.query("dust.model",dustModelString);
        dm=get_dust_module(dustModelString);
        
    }

}


double BaseEmissionLine::get_local_thermal_velocity(const BaseCell* cell)
{
    double T = cell->temperature;
    double vturb = cell->turbulent_velocity;
    double vthermal = 12.85*sqrt(T/1e4)*1000.0*100.0;//cm/s
    return sqrt(vturb*vturb + vthermal*vthermal);
    
}

double getHofz(double z)
{
  LParmParse pp;
  //this is the default cosmology, taken from Althaea.
  double H0 = 0.677900009155273E+02;
  double Omega_M = 0.307999998331070E+00;
  double Omega_L = 0.691999971866608E+00;
  
  pp.query("cosmology.H0",H0);
  pp.query("cosmology.Omega_M",Omega_M);
  pp.query("cosmology.Omega_L",Omega_L);
  return H0*sqrt(Omega_M*(z+1)*(z+1)*(z+1)+Omega_L);
}

void BaseEmissionLine::getLocalHubbleFlow(const BaseParticle& p, double v[3])
{
    double r[3];
    for(int i=0;i<3;i++)
        r[i] = p.x[i]-p.lsp[i];

    for(int i=0;i<3;i++)
        v[i]=hubble_flow_value*r[i];
}
void BaseEmissionLine::getLocalHubbleFlowAtMidpoint(const BaseParticle& p, double v[3], double pathlength)
{
    double r[3];
    for(int i=0;i<3;i++)
        r[i] = (p.x[i]+0.5*p.k[i]*pathlength)-p.lsp[i];

    for(int i=0;i<3;i++)
        v[i]=hubble_flow_value*r[i];
}
void BaseEmissionLine::setup()
{
    LParmParse pp;
    //figure out whether we have the boxsize parameter given implicitly or explicitly given.
    if(!pp.query("boxsize",boxsize))
    {
        std::string dataset_type("SphericalShell");
        pp.query("dataset_type",dataset_type);
        if(dataset_type.compare("Ramses")!=0 && dataset_type.compare("RamsesMPI")!=0)
        {
            std::cout << "BaseEmissionLine::setup(): Need boxsize parameter" << std::endl;
            std::abort();
            
        }
        
        
    }
    //figure out the hubble flow parameter.
    double hubble_flow = 0;
    double redshift = -1;
    if(pp.query("redshift",redshift))
    {
        
        if(pp.query("hubble_flow",hubble_flow))
        {
            std::cout << "Cant use both redshift and hubble flow value!" << std::endl;
            std::abort();
        }
        hubble_flow = getHofz(redshift);//km/s/Mpc.
        
    }
    else
    {
        pp.query("hubble_flow",hubble_flow);
        
    }
    //convert to cm/s/boxsize.
    bool no_hubble_flow;
    if(pp.query("no_hubble_flow",no_hubble_flow) && no_hubble_flow)
        hubble_flow = 0.0;
    if(Parallel::IOProcessor() && hubble_flow != 0.0)
        std::cout << "Using Hubble flow of " << hubble_flow << " km/s/Mpc." << std::endl;
    hubble_flow *= 1e5;//cm/s/Mpc
    hubble_flow /= (3.086e+24/boxsize);
    
    hubble_flow_value = hubble_flow;
//     std::cout << hubble_flow_value << std::endl;
        
        
        
}

double BaseEmissionLine::get_boxsize()
{
 return boxsize;   
}
void BaseEmissionLine::set_boxsize(double size)
{
    if(Parallel::IOProcessor())
        std::cout << "Setting boxsize to " << size << " cm." << std::endl;
  boxsize = size;   
}
