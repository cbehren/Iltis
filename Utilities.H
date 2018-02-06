#ifndef __UTILITIES_
#define __UTILITIES_


namespace Constants{
namespace CGS{
    const double Msun = 1.989e33;//g
    const double cval = 299792458.0*100;//cm/s
    const double mp = 1.6726219e-24;//g
    
    
}
    
    
}



void crossproduct(const double x[3],const double y[3],double result[]);
void normalize(double* vec);
double norm(const double a[3]);
double scalar(const double a[3],const double b[3]);
struct line_of_sight{
 double x[3];   
};

#endif
