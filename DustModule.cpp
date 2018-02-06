#include "DustModule.H"
#define SCATTERDEBUG 0
#include "Utilities.H"
#include "LymanAlphaLine.H"
#include "RandomNumbers.H"
#include "LightParmParse/LParmParse.H"
#include <string>
#include "Parallel.H"

double DustModule::G_Greenstein = 0.0;
double DustModule::albedo = 0.0;
DUST_MODULE_TYPE DustModule::type = DM_NULL;

double DustModule::GreensteinPhase(const double R)
{
  //return cos theta = mu for a the Heyney-Greenstein phase function. The line below comes from direct inversion of the distribution function.
  return 1.0/(2.0*G_Greenstein)*(1+G_Greenstein*G_Greenstein-((1-G_Greenstein*G_Greenstein)/(1+G_Greenstein*(2*R-1)))*((1-G_Greenstein*G_Greenstein)/(1+G_Greenstein*(2*R-1))));  
}

double DustModule::GreensteinProbability(const double mu)
{
  //returns the probability of obtaining cos theta = mu at a scattering
  double g = G_Greenstein;
  return 0.5*(1.-g*g)/pow(1.+g*g-2.*g*mu,1.5);
  //this will be 0.5 if the G_Greenstein -> 0.
}


void DustModule::AnisotropicReemission(const double oldk[],double newk[])
{
 
  double perp_vec[3];
  RNG::point_on_sphere(newk);
  crossproduct(newk,oldk,perp_vec);  
  normalize(perp_vec);
  
  double mu=GreensteinPhase(RNG::uniform());
  for(int i=0;i<3;i++)
     newk[i]=mu*oldk[i]+sqrt(1-(mu*mu))*perp_vec[i];
  normalize(newk);
}

//
//NULL IMPLEMENTATION 
//

void DustModuleNull::init()
{
  G_Greenstein = 0.0;
  type = DM_NULL;
  
}
double DustModuleNull::scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given)
{
  std::cout << "DustModuleNull::scatter(): This should never be called." << std::endl;
  std::abort();
  return 0;
}

double DustModuleNull::DustConversionFactor()
{
  return 0.0;
  
}

// int DustModuleNull::ConvertDust(MultiFab& fab,const double dx)
// {
//  return 0; 
// }
double DustModuleNull::tauConstFactor()
{
  return 0;
}
double DustModuleNull::getAlbedo()
{
  return 0.;
}
//
//VERHAMME12 IMPLEMENTATION 
//
void DustModuleVerhamme12::init()
{
  G_Greenstein = 0.7;
  albedo = 0.5;
  LParmParse pp;
  if(pp.query("dust.albedo",albedo))
   {
       if(Parallel::IOProcessor())
           std::cout << "Setting albedo to " << albedo << std::endl;
   }
  type = DM_VERHAMME12;
  
}
double DustModuleVerhamme12::scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given)
{
  double local_v_thermal = BaseEmissionLine::get_local_thermal_velocity(cell);
  const double *bulk = cell->velocity;
  const double *k = p->k;
  const double x =  LymanAlphaLine::frequency_to_x(p->frequency,local_v_thermal);
  
  double velocity[3];
  
  double hflow[3];
  BaseEmissionLine::getLocalHubbleFlow(*p,hflow);
  
  double k_new[3];
  double k_new_atomframe[3];
  double k_atomframe[3];
  
  double k_difference[3];
  if(SCATTERDEBUG)
    std::cout << "Old data: k " << k[0] << " "  << k[1] << " "  << k[2] << " freq " << x << '\n';
  
  for(int i=0;i<3;i++)
    velocity[i] = (hflow[i]+bulk[i])/BaseEmissionLine::c_val;
  LymanAlphaLine::Lorentz_Transform_in(velocity, k, k_atomframe);
  
  //get the new k-vector. 
  if(!k_and_atom_velocity_given)
  {
      AnisotropicReemission(k_atomframe,k_new_atomframe);
      LymanAlphaLine::Lorentz_Transform_out(velocity, k_new_atomframe, k_new);
  }
  else
  {
      for(int i=0;i<3;i++)
          k_new[i] = newk[i]; 
  }

  
  //get the difference of the k-vectors, save new one to the particle
  for(int i=0;i<3;i++)
    {
      k_difference[i] = (k_new[i]-k[i]);
      if(!k_and_atom_velocity_given)
        newk[i] = k_new[i];
    }
    
    //calculate new freq in observer frame
  
   if(SCATTERDEBUG)
    std::cout << "New data: k " << newk[0] << " "  << newk[1] << " "  << newk[2] << " freq " << x + scalar(k_difference,velocity)*BaseEmissionLine::c_val/local_v_thermal << '\n';
  double new_freq_obs = (x + scalar(k_difference,velocity)*BaseEmissionLine::c_val/local_v_thermal);

  return LymanAlphaLine::x_to_frequency(new_freq_obs,local_v_thermal);

}

double DustModuleVerhamme12::DustConversionFactor()
{
  //this is the conversion factor we need to go from the metal density in Nyx units to the dust density in 1/m^3 as done in Verhamme 12.
  double ratiometaldust=0.3;
  double mass_dustgrain=3.0e-17;//in g
    
  return ratiometaldust/mass_dustgrain;
  
}
double DustModuleVerhamme12::ConvertDust(double density, double temperature, double metallicity)
{
    const double mh = 1.6726219e-24;//g
    if(temperature>1e5)
        return 0.0;
    double metal_density = density*metallicity;
    double dust_conversion_factor = DustConversionFactor()*mh;
    return dust_conversion_factor*metal_density;
    
}

// int DustModuleVerhamme12::ConvertDust(MultiFab& fab,const double dx)
// {
//     //currently set up for Ramses simulations. 
//   //expected input units; hydrogen number density in 1/m^3, metallicity z.
//   //needed; dust grain number density in m^3.
//   const double mh = 1.6726219e-24;
//  
//   
//   for (MFIter mfi(fab); mfi.isValid(); ++mfi)
//   {
//     FArrayBox& b = fab[mfi];
// 
//     const Box& box = mfi.validbox();
//     const int* lo = box.loVect();
//     const int* hi = box.hiVect();
//     //next line removes dust where T>100000 K
//     BL_FORT_PROC_CALL(DUST_AND_TEMPERATURE,dust_and_temperature)
//       (BL_TO_FORTRAN(b),&Sedona::NUM_STATE, lo, hi, &dx, &Sedona::Temp,&Sedona::rhoZ);
//   }
//    //next line converts yields metal mass density in units of hydrogen atoms per m^3
//   fab.Multiply(fab,fab,Sedona::Density,Sedona::rhoZ,1,fab.nGrow());
//   //next line converts first into g/m^3, then into dust grains per 1/m^3
//   double dust_conversion_factor=DustConversionFactor()*mh;
//   fab.mult(dust_conversion_factor,Sedona::rhoZ,1,fab.nGrow());
//   return 0;
// }

double DustModuleVerhamme12::tauConstFactor()
{
  return 1.e-12*M_PI*2.0;//following verhamme 06. these units are cgs.  
}

double DustModuleVerhamme12::getAlbedo()
{
  return albedo;
}


//
//DAHLIA IMPLEMENTATION 
//
void DustModuleDahlia::init()
{
   G_Greenstein = 0.7;
   type = DM_DAHLIA;
   double a = 0.2548;
   albedo =  a/(a+1.);//from Weingartner & Draine (2001) and Li & Draine (2001) for the LMC at the Lyman alpha line center.
   LParmParse pp;
   if(pp.query("dust.albedo",albedo))
   {
       if(Parallel::IOProcessor())
           std::cout << "Setting albedo to " << albedo << std::endl;
   }
  
}
double DustModuleDahlia::scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given)
{
  double local_v_thermal = BaseEmissionLine::get_local_thermal_velocity(cell);
  const double *bulk = cell->velocity;
  const double *k = p->k;
  const double x =  LymanAlphaLine::frequency_to_x(p->frequency,local_v_thermal);
  
  double velocity[3];
  
  double hflow[3];
  BaseEmissionLine::getLocalHubbleFlow(*p,hflow);
  
  double k_new[3];
  double k_new_atomframe[3];
  double k_atomframe[3];
  
  double k_difference[3];
  if(SCATTERDEBUG)
    std::cout << "Old data: k " << k[0] << " "  << k[1] << " "  << k[2] << " freq " << x << '\n';
  
  for(int i=0;i<3;i++)
    velocity[i] = (hflow[i]+bulk[i])/BaseEmissionLine::c_val;
  LymanAlphaLine::Lorentz_Transform_in(velocity, k, k_atomframe);
  
  //get the new k-vector. 
  if(!k_and_atom_velocity_given)
  {
      AnisotropicReemission(k_atomframe,k_new_atomframe);
      LymanAlphaLine::Lorentz_Transform_out(velocity, k_new_atomframe, k_new);
  }
  else
  {
      for(int i=0;i<3;i++)
          k_new[i] = newk[i]; 
  }

  
  //get the difference of the k-vectors, save new one to the particle
  for(int i=0;i<3;i++)
    {
      k_difference[i] = (k_new[i]-k[i]);
      if(!k_and_atom_velocity_given)
        newk[i] = k_new[i];
    }
    
    //calculate new freq in observer frame
  
   if(SCATTERDEBUG)
    std::cout << "New data: k " << newk[0] << " "  << newk[1] << " "  << newk[2] << " freq " << x + scalar(k_difference,velocity)*BaseEmissionLine::c_val/local_v_thermal << '\n';
  double new_freq_obs = (x + scalar(k_difference,velocity)*BaseEmissionLine::c_val/local_v_thermal);

  return LymanAlphaLine::x_to_frequency(new_freq_obs,local_v_thermal);
}


double DustModuleDahlia::DustConversionFactor()
{
  double ratiometaldust=0.3;     
  LParmParse pp;
  pp.query("dust.dust_metal_ratio",ratiometaldust);
  //std::cout << "using " << ratiometaldust << " as the dust to metal ratio." << std::endl;
  return ratiometaldust;
  
}
/*
int DustModuleDahlia::ConvertDust(MultiFab& fab,const double dx)
{
  //expected input units; hydrogen number density in 1/m^3, metallicity z.
  //needed; dust density in g/m^3. 
  const double mh = 1.6726219e-24;
 
  
  for (MFIter mfi(fab); mfi.isValid(); ++mfi)
  {
    FArrayBox& b = fab[mfi];

    const Box& box = mfi.validbox();
    const int* lo = box.loVect();
    const int* hi = box.hiVect();
    //next line removes dust where T>100000 K
//     BL_FORT_PROC_CALL(DUST_AND_TEMPERATURE,dust_and_temperature)
//       (BL_TO_FORTRAN(b),&Sedona::NUM_STATE, lo, hi, &dx, &Sedona::Temp,&Sedona::rhoZ);
  }
  //next line converts yields metal mass density in units of hydrogen atoms per m^3
  fab.Multiply(fab,fab,Sedona::Density,Sedona::rhoZ,1,fab.nGrow());
  //next line converts into g/m^3
  double dust_conversion_factor=DustConversionFactor()*mh;
  fab.mult(dust_conversion_factor,Sedona::rhoZ,1,fab.nGrow());
  return 0;
}*/

DustModuleDahlia::DustModuleDahlia()
{
  init();
  LParmParse pp;
  //setup the absorption cross section
  std::string dust_type("MW");
  pp.query("dust.dahlia.type",dust_type);
  if(!dust_type.compare("MW"))
  {
    if(Parallel::IOProcessor())
        std::cout << "DustModuleDahlia(): Using MW dust." << std::endl;
    //from ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_MW_3.1_60_D03.all for the MW, R_V=3.1
    tauConst = 5.802E+04;
  }
  if(!dust_type.compare("LMC"))
  {
    if(Parallel::IOProcessor())
        std::cout << "DustModuleDahlia(): Using LMC dust." << std::endl;
    //from ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_LMCavg_20
    tauConst = 7.252E+04;
  }
  if(!dust_type.compare("SMC"))
  {
    if(Parallel::IOProcessor())
        std::cout << "DustModuleDahlia(): Using SMC dust." << std::endl;
    //from ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_SMCbar_0
    tauConst = 7.169E+04;
  }
  if(tauConst<0.0)
  {
    std::cout << "DustModuleDahlia(): Unknown dust.dahlia.type!" << std::endl;
    std::abort();

  }
}
double DustModuleDahlia::tauConstFactor()
{
  return 72520.0; //from Weingartner & Draine (2001) and Li & Draine (2001) for the LMC at the Lyman alpha line center. is in units of cm^2/g
  //in integrate tau, CodeUnits_Density will eat up the m^3 in the dust density.
  //to be more precise, this is the value of K_abs at 1216 Angstrom in the table on Draines website,
  //ftp://ftp.astro.princeton.edu/draine/dust/mix/kext_albedo_WD_LMCavg_20
}

double DustModuleDahlia::getAlbedo()
{
    return albedo;
 
}

double DustModuleDahlia::ConvertDust(double density, double temperature, double metallicity)
{
    const double mh = 1.6726219e-24;//g
    //TODO add temperature floor
//     if(temperature>1e5)
//         return 0.0;
    double metal_density = density*metallicity;
    double dust_conversion_factor = DustConversionFactor()*mh;
    return dust_conversion_factor*metal_density;
    
}



