#ifndef __UNIGRID_DATASET_MPI_H
#define __UNIGRID_DATASET_MPI_H
#include "BaseDataset.H"

//implements a unigrid geometry, read from a text file.

class UnigridDataset : public BaseDataset{
public:
    const BaseCell* data_at(const double x[3], const double k[3],double &pathlength,DatasetStatus &status) const;
    int setup();
    double get_dx(const double x[3],DatasetStatus &status) const {return dx;}
    int set_domain(BaseParticleVector *vp) const {return 0;}
    void get_optical_depth_stats();
    const int get_nearest_face(const double x[3],double *orth_distance, double distance[3]) const;
private:
    BaseCell* cells = NULL;
    int ngrid;
    double dx;
    
};
#endif
