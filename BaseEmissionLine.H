#ifndef __BASEEMISSIONLINE_
#define __BASEEMISSIONLINE_
#include "DustModule.H"

//defines the base class for an emission line. Virtual functions need to be implemented in each derived class. Besides, some helper functions are defined.


class BaseCell;
class BaseParticle;
class BaseEmissionLine{
public:
    BaseEmissionLine();
    //setup the Emission line object. Some functionality might be moved to the constructor.
    virtual void setup();
    //convert frequency to (delta) wavelength. Actually we would not need to make this virtual. TODO implement this for the base class
    virtual double convert_wavelength(double frequency)const = 0;
    //convert a line shift in km/s to a physical frequency. TODO implement this for the base class, no need for it to be virtual
    virtual double shift_to_frequency (double shift_in_kms)const= 0;
    //return the line center frequency of the line
    virtual double Nu0() const = 0;
    //given a hydro cell, a particle, and a path length, return the optical depth seen in this cell. returns the fraction steming from dust as well.
    virtual double get_optical_depth(const BaseCell* cell, const BaseParticle* p, double pathlength, double *fraction_in_dust) const = 0;
    //given a cell and a particle, execute the scattering. the boolean flag is necessary since for the peeling off photons, we do know the outgoing k vector and the atom velocity.
    virtual double scatter(const BaseCell* cell,const BaseParticle* p,double newk[], const double tau0,double atom_velocity[],bool k_and_atom_velocity_given) = 0;
    //returns the local thermal velocity, including any turbulent contribution, in cm/s
    static double get_local_thermal_velocity(const BaseCell* cell);
    //given cell, particle, and pathlength, returns the optical depth in the line center, and the fraction of optical depth stemming from dust. 
    virtual double get_optical_depth_line_center(const BaseCell* cell, const BaseParticle* p, double pathlength, double *fraction_in_dust) const = 0;
    //returns the Hubble flow given a particle by calculating H(z)*r, where r is the distance between the particle and its last point of interaction.
    static void getLocalHubbleFlow(const BaseParticle& p, double v[3]);
    //same, but evaluating at midpoint between last point of interaction and current location of particle
    static void getLocalHubbleFlowAtMidpoint(const BaseParticle& p, double v[3], double pathlength);
    //velocity of light in cm/s TODO: remove and replace with the Constants namespace
    static const double c_val;//cm/s
    //returns the boxsize in cm
    static double get_boxsize();
    //sets the boxsize in cm
    static void set_boxsize(double size);
    //holds the DustModule
    static DustModule *dm;


    
protected:
    //boxsize in cm TODO: probably should be moved to a different place, but dunno where.
    static double boxsize;//cm
    //the hubble flow rate in cm/s/cm
    static double hubble_flow_value;//cm/s/cm
    

    
};
#endif
