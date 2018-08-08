#ifndef __BASEDATASET_
#define __BASEDATASET_

//defines the base class for datasets. Each type of dataset needs to implement the routines specified below.

#include <cmath>
#include <iostream>
#include "BaseCell.H"
#include "BaseParticleVector.H"
#ifdef USE_MPI
//include MPI headers
#endif
//a enum datatype for signaling whether a position was within the data. ALL_GOOD = data has been returned, NOT_YOUR_DOMAIN = position is valid, but data is on a different process, DS_NOT_IN_DOMAIN = position is outside of computational domain.
enum DatasetStatus {DS_ALL_GOOD, DS_NOT_IN_DOMAIN, DS_NOT_YOUR_DOMAIN};

class BaseDataset{
public:
    //given a position x and a k vector, get the hydro data at x. Also returns the pathlength in this cell given x and k, and the status.
    virtual const BaseCell* data_at(const double x[3], const double k[3],double &pathlength,DatasetStatus &status) const = 0;
    virtual double get_dx(const double x[3],DatasetStatus &status) const = 0;
    //setup the Dataset, e.g. loading files, converting units, etc. Some things might be better off in the constructor.
    virtual int setup() = 0 ;
    //relevant only for datasets with distributed data; given a vector with particles, update each particles domain attribute to indicate which process should get this particle.
    virtual int set_domain(BaseParticleVector *vp) const =0;
    //a routine relevant only for the Neufeld acceleration scheme implemented for the Lyman alpha line. returns the distance vector and surface normal of the nearest face of the cell in which x lies.
    virtual const int get_nearest_face(const double x[3],double *orth_distance, double distance[3]) const {std::cout << "Not implemented" << std::endl;std::abort();}
    //the maximum pathlength (in code units) allowed in a single step.
    double max_step=1e-6;
    //the minimum pathlength allowed for a single step (in code units)
    double min_step=1e-10;

};
#endif
