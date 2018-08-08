#include "NeutralFractionModule.H"
#include "math.h"
#include <vector>
#include "LightParmParse/LParmParse.H"
#include "Parallel.H"


NeutralFractionRahmati13::NeutralFractionRahmati13(double z)  : NeutralFractionModule()
{
    redshift = z;
    init();
    init_read();
}
NeutralFractionRahmati13::~NeutralFractionRahmati13()
{
    if(interp_gamma)
        delete interp_sigma;
    if(interp_gamma)
        delete interp_gamma;
}

int NeutralFractionRahmati13::init_read()
{
    
    FILE *fp = fopen( "HM12_sigma_gamma.dat","r");
    if(!fp)
    {
        std::cout << "NeutralFractionRahmati13:init(): Could not open data file!" << std::endl;
        std::abort();
    }
    char line[1000];
    //remove the header line
    fgets(line,1000,fp);
    std::vector<double> redshifts;
    std::vector<double> gamma;
    std::vector<double> sigma;
    while(fgets(line,1000,fp))
    {
        double x,y,z;
        if(sscanf(line,"%le %le %le\n",&x,&y,&z)!=3)
        {
            std::cout << "NeutralFractionRahmati13:init(): Wrong data in line!" << std::endl;
            std::abort();
        }
        redshifts.push_back(x);
        sigma.push_back(y);
        gamma.push_back(z);
    
    }
    fclose(fp);
    interp_gamma = new Interpolate(redshifts.data(),gamma.data(),gamma.size());
    interp_sigma = new Interpolate(redshifts.data(),sigma.data(),sigma.size());
    return 0;
    
}
double NeutralFractionRahmati13::neutralFraction(double density, double temperature)
{
    //eq A8
    //INPUT: density in number density (cgs), temperature in K
    //OUTPUT: neutral fraction 
    double alpha = alpha_A(temperature);
    double nssh = nH_ssh(temperature,redshift);
    double lambda = Lambda_T(temperature);
    double gamma = Gamma_Phot(redshift,density,nssh);
    double A = alpha+lambda;
    double B = 2*alpha+lambda+gamma/density;
    double C = alpha;
    
    double eta = (B - sqrt(B*B - 4*A*C))/(2*A);
    return eta;
}




double NeutralFractionRahmati13::alpha_A(double T)
{
    // eq. A3 in R+13
    // case A recombination rate coefficient
    double lambda = 315614.0/T;
    double c = 1.269e-13;
    double up = pow(lambda,1.503);
    double down = pow(1+pow(lambda/0.522,0.47),1.923);

    return c*up/down;
}
double NeutralFractionRahmati13::Lambda_T(double T)
{
    //eq A6 in R+13
    double c = 1.17e-10;
    double up = sqrt(T)*exp(-157809.1/T);
    double down = 1+sqrt(T/1.e5);
   
    return c*up/down;
}
double NeutralFractionRahmati13::nH_ssh(double T, double redshift)
{
    //eq 13
    double fg = 0.450000017881393E-01/0.307999998331070E+00; //correct for althaea runs.
    double c = 6.73e-3;
    double T4 = T/1.e4;
    double sigma = sigma_vh(redshift);
    double Gamma = Gamma12(redshift);
    double term1 = pow(sigma/2.49e-18,-2./3.);
    double term2 = pow(T4,0.17)*pow(Gamma,2./3.)*pow(fg/0.17,-1./3.);
    return c*term1*term2;
}
double NeutralFractionRahmati13::Gamma_Phot(double redshift,double nh,double nh_ssh)
{
    double f = nh/nh_ssh;
    double term1 = 0.98*pow((1+pow(f,1.64)),-2.28);
    double term2 = 0.02*pow(1+f,-0.84);
    double Gamma = Gamma12(redshift)*1.e-12;
    return (term1+term2)*Gamma;
}
double NeutralFractionRahmati13::Gamma12(double redshift)
{
    return interp_gamma->y(redshift)/1.0e-12;
}
double NeutralFractionRahmati13::sigma_vh(double redshift)
{
    return interp_sigma->y(redshift); 
}


NeutralFractionChardin17::NeutralFractionChardin17(double z) : NeutralFractionRahmati13(z)
{
    init();
    if(z<3.00)
    {
        std::cout << "NeutralFractionChardin17::NeutralFractionChardin17(): supporting only z>3.0" << std::endl;
        std::abort();
    }
    init_read();
}

NeutralFractionChardin17::~NeutralFractionChardin17()
{
    if(interp_n0)
        delete interp_n0;
    if(interp_alpha1)
        delete interp_alpha1;
    if(interp_alpha2)
        delete interp_alpha2;
    if(interp_beta)
        delete interp_beta;
    if(interp_f)
        delete interp_f;
}

int NeutralFractionChardin17::init_read()
{
    //read gamma.
    NeutralFractionRahmati13::init_read();
    //this is the hardcoded data from table A1 in Chardin+17 
    
    double z[] = {3,4,5,6,7,8,9,10};
    double n0[] = {0.0090346,0.009346,0.010379,0.006955,0.002658,0.003974,0.004598,0.004694};
    double alpha1[] = {-1.115653,-0.950010,-1.294665,-0.941372,-0.866887,-0.742237,-0.642645,-0.386421};
    double alpha2[] = {-1.648170,-1.503310,-1.602099,-1.507124,-1.272957,-1.397100,-1.212105,-0.857335};
    double beta[] = {5.324900,5.872852,5.056422,6.109438,7.077863,7.119987,9.988387,12.94204};
    double f[] = {0.018174,0.015308,0.024356,0.028830,0.040894,0.041213,0.028581,0.005710};
    
    interp_n0 = new Interpolate(z,n0,8);
    interp_alpha1 = new Interpolate(z,alpha1,8);
    interp_alpha2 = new Interpolate(z,alpha2,8);
    interp_beta = new Interpolate(z,beta,8);
    interp_f = new Interpolate(z,f,8);
    return 0;    
}


double NeutralFractionChardin17::nH_ssh(double T, double redshift)
{
    return interp_n0->y(redshift);
    
}
double NeutralFractionChardin17::Gamma_Phot(double redshift,double nh,double nh_ssh)
{
    double f = interp_f->y(redshift);
    double alpha1 = interp_alpha1->y(redshift);
    double alpha2 = interp_alpha2->y(redshift);
    double beta = interp_beta->y(redshift);
    double Gamma = Gamma12(redshift)*1.e-12;
    double nr = nh/nh_ssh;
    
    double term1 = (1.-f)*pow(1.+pow(nr,beta),alpha1);
    double term2 = f*pow(1.+nr,alpha2);
    
    return (term1+term2)*Gamma;
}

#ifdef TEST_STANDALONE
int main(int argc,char* args[])
{
    double redshift = 6.01;
    NeutralFractionRahmati13 ra(redshift);
    NeutralFractionChardin17 ch(redshift);

    double nh_ssh = ra.nH_ssh(1e4,redshift);
    std::cout << "#log(nh) f_R+13 f_C+17 " << nh_ssh << " " << ra.Gamma12(redshift)*1.e-12 << "\n";
    for(int i=-60;i<40;i++)
    {
        double n = pow(10.0,(double)i/10.0);
        double T = 1.e4;
        std::cout << n << " " << ra.neutralFraction(n,T) << " " << ra.Gamma_Phot(redshift,n,nh_ssh) <<" " << ch.neutralFraction(n,T)  << std::endl;
    }
    return 0;
}
#endif


