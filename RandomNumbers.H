#include <random>
//this is just a thread-safe wrapper around the random templates from C++11. We make sure that each process/thread gets its own, unique RNG with a unique seed.
class RNG{
public:
    //return a uniformly distributed pseudo random number
    static double uniform(int tid=-1);
    //return a gaussian distributed pseudo random number
    static double gaussian(double mean,double sigma,int tid=-1);
    //return a uniformly distributed point on the unit sphere
    static void point_on_sphere(double x[3]);
    //return a exponentially distributed pseudo random number
    static double exponential(int tid=-1);
    //return a random integer between L and R
    static int uniform_integer(int L,int R,int tid=-1);
    //initialize the RNG
    static void initialize();
private:
    static std::default_random_engine* generators;
    static std::uniform_real_distribution<double>* distributions;
    static int nthreads;
    static int nprocs;
};
