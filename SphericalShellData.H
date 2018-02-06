#ifndef __SPHERICAL_SHELL_H
#define __SPHERICAL_SHELL_H
#include "BaseDataset.H"

//implements a spherical shell geometry.

struct ShellParameters{
  double boxsize;
  double inner_radius;
  double outer_radius;
  double column_density;
  double column_density_dust;
  double temperature;
  double outflow_velocity;
  
  double density;
  double density_dust;
    
};

class SphericalShellData : public BaseDataset{
public:
    const BaseCell* data_at(const double x[3], const double k[3],double &pathlength,DatasetStatus &status) const;
    double get_dx(const double x[3],DatasetStatus &status) const {return (params->outer_radius-params->inner_radius)/1e4;}
    int setup();
    int set_parameters(ShellParameters sh);
    int set_domain(BaseParticleVector *vp) const {return 0;}
private:
    ShellParameters* params = NULL;
    BaseCell* cells = NULL;
    //calculates intersections between spherical shell and particles.
    double intersection(const double x[3], const double k[3], double radius) const;
};
#endif
