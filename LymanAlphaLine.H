#include "BaseEmissionLine.H"
//defines an emission line object for the lyman alphe line.
class BaseDataset;
//the following defines control the acceleration scheme used.
//define to enable acceleration:
#define USE_ACCELERATION
//define to use the "full" scheme, with a fixed value for the critical frequency. undefine to use the Laursen scheme with an adaptive critical frequency.
#define FULL_ACC
//the fixed value for the critical frequency.
#define ACC_XCRIT 3.0
//define to use dipole scattering above |x|>x_crit_dipole:
//#define DIPOLE_SCATTERING
//set the value of x_crit_dipole
#define X_CRIT_DIPOLE 3.0

//enum type for Neufeld acceleration scheme
enum NeufeldAccStatus{RegularScattering, Success, Absorbed};
class LymanAlphaLine : public BaseEmissionLine{
public:
    
    void setup();
    //get the Voigt profile
    static double get_voigt(double x, double a);
    //get the atom velocity of the scattering atom; uses rejection method internally
    static int get_atom_velocity(double vel[], double b, const double k[],double freq, const double atimestau0);
    //convert frequency to dimensionless frequency
    static double frequency_to_x(double freq, double vthermal);
    //and back
    static double x_to_frequency(double x, double vthermal);
    double shift_to_frequency(double shift_in_kms) const;//converts v in km/s to delta nu in 1/s
    double Nu0() const {return nu_0;};
    double convert_wavelength(double frequency) const;
    

    
    double get_optical_depth(const BaseCell* cell, const BaseParticle* p, double pathlength, double *fraction_in_dust) const;
    double get_optical_depth_line_center(const BaseCell* cell, const BaseParticle* p, double pathlengt, double *fraction_in_dust) const;
    static double frequency_to_wavelength(double frequency);

    //scatter on hydrogen atom. atom velocity, new k vector and new frequency are determined in the function. if k_and_atom_velocity_given is set, then both the new k vector and the atom velocity are read from the arrays provided as arguments
    double scatter(const BaseCell* cell,const BaseParticle* p, double newk[], const double tau0,double atom_velocity[],bool k_and_atom_velocity_given);
    //leaves cell, skipping scattering within, and sets frequency according to analytic (Neufeld) solution for homogeneous sphere.
    double scatter_and_leave(const BaseCell *cell,const BaseParticle* p,const BaseDataset* ds,double newk[],double new_position[3],NeufeldAccStatus *success);
    
    static const double lambda_0;//in angstrom
    static const double nu_0;//1/s
    
    
    static void Lorentz_Transform_in(const double velocity[],const double in[],double out[]);
    static void Lorentz_Transform_out(const double velocity[],const double in[],double out[]);
    //draw a random frequency consistent with the Neufeld solution for the homogeneous sphere.
    double draw_neufeld_frequency(double a, double tau0);
    //draw from the distribution of mu = cos theta for the direction of exit of a cell, see Tasitsiomi 2006
    double draw_exit_direction();
    //draw from dipole distribution 
    double draw_mu_dipole();

private:
    const double eta=0.71;//parameter used in the Neufeld acceleration scheme adopted from Laursen 08/09
   
    
};
