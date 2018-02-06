#include <iostream>
enum BOUNDARY_COND {CONSTANT,ABORT,PERIODIC,LINEAR};

class Interpolate{
public:
  Interpolate(double *x,double *y,long n, BOUNDARY_COND leftedge_bc=ABORT, BOUNDARY_COND rightedge_bc=ABORT);
  Interpolate();
  ~Interpolate();
  void Init(double *x,double *y,long n, BOUNDARY_COND leftedge_bc=ABORT, BOUNDARY_COND rightedge_bc=ABORT);
  double y(double x);
private:
  double *X;
  double *Y;
  long N;
  BOUNDARY_COND bc_left;
  BOUNDARY_COND bc_right;
  
  void init();
  double left_boundary(double x);
  double right_boundary(double x);
  double get_boundary(double x);

  
  
  
};