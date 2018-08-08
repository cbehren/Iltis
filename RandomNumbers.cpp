#include "RandomNumbers.H"
#include <cmath>
#ifdef _OPENMP
#include "omp.h"
#endif
#include "Parallel.H"
#include "LightParmParse/LParmParse.H"


int RNG::nthreads = 1;
int RNG::nprocs = 1;

std::uniform_real_distribution<double>* RNG::distributions = NULL;
std::default_random_engine* RNG::generators = NULL;

 double RNG::uniform(int tid)
{
    if(tid>-1)
        return distributions[tid](generators[tid]);
    else
    {
        tid = 0;
    #ifdef _OPENMP
        tid = omp_get_thread_num(); 
    #endif
        return distributions[tid](generators[tid]);
    }
}

 void RNG::point_on_sphere(double x[3])
{
 int tid = 0;
#ifdef _OPENMP
 tid = omp_get_thread_num();
#endif
  double a,b,c;
  do{
    a=2.0*RNG::uniform(tid)-1.0;
    b=2.0*RNG::uniform(tid)-1.0;
    c=2.0*RNG::uniform(tid)-1.0;
  } while (a*a+b*b+c*c>1.0 || ((a==0.0)&&(b==0.0)&&(c==0.0)));
  double norm = sqrt(a*a+b*b+c*c);
  x[0]=a/norm;
  x[1]=b/norm;
  x[2]=c/norm;
  return;
}
     

 void RNG::initialize()
{
    LParmParse pp("rng");
    long seed;
    pp.get("seed",seed);
    int threads = 1;
#ifdef _OPENMP
   threads = omp_get_max_threads(); 
#endif
   int myproc= Parallel::MyProc();
   int procs = Parallel::NProcs();
   
   int base_seed = myproc*threads+seed;
   
   nthreads = threads;
   nprocs = procs;
   std::uniform_real_distribution<double> temp(0.0,1.0);
   generators = new std::default_random_engine[nthreads];
   distributions = new std::uniform_real_distribution<double>[nthreads];
   for(int i=0;i<nthreads;i++)
   {
       generators[i].seed(base_seed+i);
       distributions[i].param(temp.param());
   }
}


double RNG::exponential(int tid)
{
    return -1.0*log(uniform(tid));
}

double RNG::gaussian(double mean,double sigma,int tid)
{
 double u1=uniform(tid);
 double u2=uniform(tid);
 double z=sqrt(-2*log(u1))*cos(2*M_PI*u2);
 return z*sigma+mean;
}

int RNG::uniform_integer(int L,int R,int tid)
{
    std::uniform_int_distribution<int> distribution(L,R);
    if(tid>-1)
        return distribution(generators[tid]);
    else
    {
        tid = 0;
    #ifdef _OPENMP
        tid = omp_get_thread_num(); 
    #endif
        return distribution(generators[tid]);
    }
}
