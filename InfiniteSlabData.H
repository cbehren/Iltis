#ifndef __INFINITESLAB_H
#define __INFINITESLAB_H
#include "BaseDataset.H"

//this defines an "infinite slab" geometry, that is, a slab infinitely extended in y and z, but finite in x.

struct SlabParameters{
  double boxsize;
  double column_density;
  double column_density_dust;
  double temperature;
  
  double density;
  double density_dust;
    
};

class InfiniteSlabData : public BaseDataset{
public:
    const BaseCell* data_at(const double x[3], const double k[3],double &pathlength,DatasetStatus &status) const;
    double get_dx(const double x[3],DatasetStatus &status) const {return 1./1e4;}
    int setup();
    int set_parameters(SlabParameters sh);
    int set_domain(BaseParticleVector *vp) const {return 0;}
private:
    SlabParameters* params = NULL;
    BaseCell* cells = NULL;
    double intersection(const double x[3], const double k[3]) const;
    
};
#endif
