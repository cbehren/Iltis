#include "LymanAlphaLine.H"
#include "RejectionMethod.H"
#include "Utilities.H"
#include "RandomNumbers.H"
#include <cmath>

#ifdef _OPENMP
#include "omp.h"
#endif

#include "BaseCell.H"
#include "BaseParticle.H"

#include "LightParmParse/LParmParse.H"
#include "Parallel.H"

const double LymanAlphaLine::lambda_0 = 1215.668;//in angstrom
const double LymanAlphaLine::nu_0 = c_val/(lambda_0/1.0e8);//1/s

void LymanAlphaLine::setup()
{
    BaseEmissionLine::setup();
#ifndef USE_ACCELERATION
    if(Parallel::IOProcessor())
        std::cout << "LymanAlphaLine::setup(): NOT using acceleration scheme" << std::endl; 
#else
    #ifdef FULL_ACC
    if(Parallel::IOProcessor())
        std::cout << "LymanAlphaLine::setup(): Using static acceleration scheme with xcrit="<<  ACC_XCRIT << std::endl; 
    #else
    if(Parallel::IOProcessor())
        std::cout << "LymanAlphaLine::setup(): Using dynamic acceleration scheme" << std::endl; 
    #endif
#endif
    
}

double LymanAlphaLine::get_optical_depth_line_center(const BaseCell* cell, const BaseParticle* p, double pathlength, double *fraction_in_dust) const
{
    const double const_tau0 = boxsize*8.3e6/2.0e20; //following faucher-giguere 09. these units are cgs.
    const double density = cell->density;
    const double temperature = cell->temperature;
    const double dust_density = cell->dust_density;
   
    double tau0 = density *const_tau0 * pathlength * 1.0/(sqrt(temperature/20000.0));//following Laursen 06
    double tau0_dust = dust_density * dm->tauConstFactor() * pathlength * boxsize;
    if(tau0>0.0)
        *fraction_in_dust = tau0_dust/(tau0+tau0_dust);
    else
        *fraction_in_dust = 1.0;
    return tau0+tau0_dust;
    
}

double LymanAlphaLine::get_optical_depth(const BaseCell* cell, const BaseParticle* p, double pathlength, double *fraction_in_dust) const
{
    const double const_tau0 = boxsize*8.3e6/2.0e20; //following faucher-giguere 09. these units are cgs.
    double total_velocity[3];
    const double density = cell->density;
    const double *velocity = cell->velocity;
    const double temperature = cell->temperature;
    const double local_v_thermal = get_local_thermal_velocity(cell);
    const double dust_density = cell->dust_density;
    const double a_param = 4.7e-4*(12.85/(local_v_thermal/1.e5));//following Verhamme 2008   
    const double *k = p->k;
    const double x = frequency_to_x(p->frequency,local_v_thermal);
    if(hubble_flow_value!=0)
    {
        getLocalHubbleFlowAtMidpoint(*p,total_velocity,pathlength);
        for(int i=0;i<3;i++)
            total_velocity[i]+=velocity[i];
    }
    else
        for(int i=0;i<3;i++)
            total_velocity[i]=cell->velocity[i];
        
    const double x_fluid = x - scalar(k,total_velocity)/local_v_thermal;
    //TODO substepping for hubble flow.
    double tau0 = density *const_tau0 * pathlength * 1.0/(sqrt(temperature/20000.0));//following Laursen 06
    double tau = tau0 * get_voigt(x_fluid,a_param);
    double tau_dust = dust_density * (dm->tauConstFactor() * pathlength * boxsize);
    if(tau+tau_dust>0.0)
        *fraction_in_dust = tau_dust/(tau+tau_dust);
    else
        *fraction_in_dust = 0.0;
    return tau+tau_dust;
    
    
}
double LymanAlphaLine::get_voigt(double x,double a) //following Laursen et al 2009, 08050.3153
{

  double q=0;
  double z=(x*x-0.855)/(x*x+3.42);
  

  if(z>0)
    q=(1.0+21.0/(x*x))*(a/(M_PI*(x*x+1)))*(5.674*(z*z*z*z)-9.207*(z*z*z)+4.421*(z*z)+0.1117*z);
  return (q*sqrt(M_PI)+exp(-(x*x)));
  
  
}

int LymanAlphaLine::get_atom_velocity(double vel[], double local_v_thermal, const double k[],double freq, const double atimestau0)
{
  
 int tid = 0;
#ifdef _OPENMP
 tid = omp_get_thread_num();
#endif
  double a=4.7e-4*(12.85/(local_v_thermal/1.e5));//following Verhamme 2008   
  
  //thermal velocity is decomposed into components orthogonal and parallel to the direction of the photon.
  
  // generate two vectors orthogonal two each other and orthogonal to k.
  double perp_vec1[3];
  double perp_vec2[3];
  if(k[0]==0 && k[1]==0)
    {
      perp_vec1[0]=0;
      perp_vec1[1]=1;
      perp_vec1[2]=0;
    }
  else if(k[0]!=0)
      {
	perp_vec1[0]=-1*k[1]/k[0];
	perp_vec1[1]=1;  
	perp_vec1[2]=0;
       
      }
      
    else
      {
        perp_vec1[0]=1;
        perp_vec1[1]=-1*k[0]/k[1];
        perp_vec1[2]=0;
      }
      
  normalize(perp_vec1);
      
  crossproduct(perp_vec1,k,perp_vec2);  
  
  // get the parallel velocity.
//   std::cout << "input to rejection_method: " << freq << " " << a << std::endl;
  double parallel_velocity=rejection_method(freq,a)*local_v_thermal;
#ifdef USE_ACCELERATION
  #ifndef FULL_ACC
  double x_crit = 0.0;
  double x_cw = 1.59 -0.60*log10(a)-0.03*log10(a)*log10(a);
   if(atimestau0<=1)
	x_crit=0;
      else
      {
	if(atimestau0<=60)
	  x_crit=0.02*exp(0.6*pow(log(atimestau0),1.2));
	else
	  x_crit=0.02*exp(1.4*pow(log(atimestau0),0.6));
      }
  if(fabs(freq)>=x_cw)
    x_crit=0.0;
  #else    
  double x_crit = ACC_XCRIT;
  if(fabs(freq)>=x_crit)
    x_crit=0.0;
  #endif  
#else
  double x_crit = 0.0; //stub for upcoming acceleration scheme
#endif
  //get the orthogonal components
  double perp_velocity=sqrt(x_crit*x_crit-log(RNG::uniform(tid)));
  double r = RNG::uniform(tid);
  for(int i=0;i<3;i++)
  {
    perp_vec1[i]*=perp_velocity*local_v_thermal*cos(2*M_PI*r);
    perp_vec2[i]*=perp_velocity*local_v_thermal*sin(2*M_PI*r);
  }

  //write them into vel.
  for(int i=0;i<3;i++)
  {
    vel[i] = (perp_vec1[i]+perp_vec2[i]+parallel_velocity*k[i]);
//    std::cout << "atom velocity " << vel[i] << std::endl; 
  }
  //This is the THERMAL velocity of the scattering atom
  return 0;
}

double LymanAlphaLine::frequency_to_x(double freq,double vthermal)
{
    return freq*(c_val/nu_0)/vthermal;
    
}
double LymanAlphaLine::x_to_frequency(double x, double vthermal)
{
    return x*vthermal*nu_0/c_val;
}

#define SCATTERDEBUG 0
double LymanAlphaLine::scatter(const BaseCell* cell,const BaseParticle* p,double newk[], const double tau0,double atom_velocity[],bool k_and_atom_velocity_given=false)
{
    //TODO ADD hubble flow here!
  double local_v_thermal = get_local_thermal_velocity(cell);
  const double *bulk = cell->velocity;
  const double *k = p->k;
  const double x = frequency_to_x(p->frequency,local_v_thermal);
  const double a_param=4.7e-4*(12.85/(local_v_thermal/1.e5));
  const double atimestau0 = tau0 * a_param;
  double velocity[3];
  double bulk_and_hflow[3] = {0,0,0} ;
  double k_new_atomframe[3];
  double k_new[3];
  double k_atomframe[3];
  double k_difference[3];
//   double perp_vec1[3];
  
  if(hubble_flow_value!=0)
  {
      getLocalHubbleFlow(*p,bulk_and_hflow);
      for(int i=0;i<3;i++)
          bulk_and_hflow[i]+=bulk[i];
  }
  else
      for(int i=0;i<3;i++)
          bulk_and_hflow[i]=bulk[i];
  const double x_fluid = x - scalar(k,bulk_and_hflow)/local_v_thermal;

  
  if(SCATTERDEBUG)
    std::cout << "Old data: k " << k[0] << " "  << k[1] << " "  << k[2] << " freq " << x << '\n';
//    int tid = 0;
// #ifdef _OPENMP
//  tid = omp_get_thread_num();
// #endif
  if(local_v_thermal<0)
    std::abort();
  if(atimestau0<0)
    std::abort();

   if(!k_and_atom_velocity_given)
   {
    get_atom_velocity(velocity, local_v_thermal, k, x_fluid,atimestau0);
    
    //go into the atom's restframe, adding the bulk velocity first
    for(int i=0;i<3;i++)
        velocity[i] = (velocity[i] + bulk_and_hflow[i])/c_val;
   }
   else
   {
    for(int i=0;i<3;i++)
        velocity[i] = atom_velocity[i];
   }
  Lorentz_Transform_in(velocity, k, k_atomframe);
  
  //now get the new k-vector. We only use isotropic scattering here.
  if(k_and_atom_velocity_given)
  {
      Lorentz_Transform_in(velocity,newk,k_new_atomframe);
  }
  else
  {
      RNG::point_on_sphere(k_new_atomframe);      
  }
          
//TODO reimplement dipole  
// #ifdef DIPOLE
//   if(fabs(x_fluid)>=3 && !k_given)
//   {
//       double mu=0.0;
//      normalize(k_new_atomframe);
//       crossproduct(k_new_atomframe,k_atomframe,perp_vec1);  
//       normalize(perp_vec1);
// 
//       mu=solve_cubic(1.0,3,4.0-8.0*rn[tid].d2_value());
//       for(int i=0;i<3;i++)
//         k_new_atomframe[i]=mu*k[i]+sqrt(1-(mu*mu))*perp_vec1[i];
//   }
// #endif

  normalize(k_new_atomframe);
  
  //transform out
  Lorentz_Transform_out(velocity, k_new_atomframe, k_new);
  
  normalize(k_new);

  //get the difference of the k-vectors, save new one to the particle
  for(int i=0;i<3;i++)
    {
      k_difference[i] = (k_new[i]-k[i]);
      //TODO parantheses missing below
      if(!k_and_atom_velocity_given)
        newk[i] = k_new[i];
      atom_velocity[i] = velocity[i];
    }
    
    //calculate new freq in observer frame
  
   if(SCATTERDEBUG)
    std::cout << "New data: k " << newk[0] << " "  << newk[1] << " "  << newk[2] << " freq " << x + scalar(k_difference,velocity)*c_val/local_v_thermal << '\n';
  
  double new_freq_obs = (x + scalar(k_difference,velocity)*c_val/local_v_thermal);

   
  return x_to_frequency(new_freq_obs,local_v_thermal);
}



void LymanAlphaLine::Lorentz_Transform_in(const double velocity[],const double in[],double out[])
{
  int i;
  double dotproduct;
  double abs = sqrt(velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]);
  double velocity_n[3],vec[3];
  for(i=0;i<3;i++)
    velocity_n[i]=velocity[i]/std::max(abs,1e-30);
  dotproduct=scalar(velocity_n,in);
  double gamma=1.0/sqrt(1.0-dotproduct*dotproduct);
  for(i=0;i<3;i++)
    vec[i] = (in[i]-velocity_n[i]*dotproduct)*gamma;
  double fraction = abs/c_val;
  double dcosp = (dotproduct-fraction)/(1.0-dotproduct*fraction);
  double dcosp_oth = (sqrt(1.0-dcosp*dcosp));
  for(i=0;i<3;i++)
    out[i] = dcosp*velocity_n[i]+dcosp_oth*vec[i];
}

void LymanAlphaLine::Lorentz_Transform_out(const double velocity[],const double in[],double out[])
{
  int i;
  double dotproduct;
  double abs = sqrt(velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]);
  double velocity_n[3],vec[3];
  for(i=0;i<3;i++)
    velocity_n[i]=velocity[i]/std::max(abs,1e-30);
  dotproduct=scalar(velocity_n,in);
  double gamma=1.0/sqrt(1.0-dotproduct*dotproduct);
  for(i=0;i<3;i++)
    vec[i] = (in[i]-velocity_n[i]*dotproduct)*gamma;
  double fraction = abs/c_val;
  double dcosp = (dotproduct+fraction)/(1.0+dotproduct*fraction);
  double dcosp_oth = (sqrt(1.0-dcosp*dcosp));
  for(i=0;i<3;i++)
  {  out[i] = dcosp*velocity_n[i]+dcosp_oth*vec[i];
     if(out[i]!=out[i])
       std::cout << dcosp << " vel_n "<< velocity_n[i] <<" abs "<< abs <<" dcosp_oth "<< dcosp_oth << " vel "<< vec[i] << " dot "<< dotproduct << " gamma " << gamma << " in "<< in[i] <<'\n';
  } 
  
  

}

double LymanAlphaLine::frequency_to_wavelength(double frequency)
{
 
  return (c_val/(frequency+nu_0)-lambda_0/1.0e8)*1.0e8; 
//    return frequency;
}
double LymanAlphaLine::convert_wavelength(double frequency) const
{
    return frequency_to_wavelength(frequency);
    
}

double  LymanAlphaLine::shift_to_frequency(double shift_in_kms) const//converts v in km/s to delta nu in 1/s
{
  return nu_0/(shift_in_kms/(c_val/1.e5)+1)-nu_0;
  
}

