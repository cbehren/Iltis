#ifndef __BASEPARTICLEVECTOR_
#define __BASEPARTICLEVECTOR_
#include <vector>
#include "BaseParticle.H"

//defines a vector of BaseParticles and adds some helper functions.

class BaseParticleVector : public std::vector<BaseParticle>{
public:
    BaseParticleVector();
    //sort the particles by particle status.
    void sort_by_status();
    //sort by order. currently unused.
    void sort_by_order();
    //key  function to exchange particles between processes for distributed data sets
    void exchange();
    //debugging function to check whether exchange went well.
    int check_consistency();
    //prints the load balance, i.e. the number of particles per process.
    int print_balance();
    
protected:
    //true if we have a datset with distributed data
    static bool distributed_data;
    
};
#endif
