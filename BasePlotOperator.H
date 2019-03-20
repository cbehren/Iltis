#ifndef BASE_PLOT_OPERATOR_H
#define BASE_PLOT_OPERATOR_H

#include "BaseCell.H"
#include "BaseEmissionLine.H"
#include "BaseParticle.H"
#include "Utilities.H"
#include <string>
#include "LightParmParse/LParmParse.H"

class BasePlotOperator{
public:
  BasePlotOperator(const BaseEmissionLine* l);
  virtual ~BasePlotOperator(){}
  virtual void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx) = 0;
  double get_value_normed(int ipix);
  virtual void exchange();
  virtual void write_results(int index);
  virtual bool is_slice (){return false;}
  //this IS bad, but in case we use the point sampling we need to be able to tell the operators what the total number of points we use is.
  static int override_tot_pixels;

protected:
  static const BaseEmissionLine* line;
  static bool is_point_sampling;
  bool has_only_constant_normalization=false;
  std::vector<double> contributions;
  std::vector<double> normalization;  
  double constant_normalization = 1.0;
  std::string name;
  int npixels = 0;
};
#endif
