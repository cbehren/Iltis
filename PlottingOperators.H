#include "BaseCell.H"
#include "BaseEmissionLine.H"
#include "BaseParticle.H"
#include "Utilities.H"
#include <string>
#include "LightParmParse/LParmParse.H"
//defines a number of functions for plotting. 
double columnDensity(const BaseCell* cell,double pathlength, double dx, double& weight);
double MWDensity(const BaseCell* cell,double pathlength, double dx, double& weight);


class BasePlotOperator{
public:
  BasePlotOperator(const BaseEmissionLine* l);
  virtual void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx) = 0;
  double get_value_normed(int ipix);
  virtual void write_results();
protected:
  static const BaseEmissionLine* line;
  std::vector<double> contributions;
  std::vector<double> normalization;  
  double constant_normalization = 1.0;
  std::string name;
  int npixels = 0;
};

class columnDensityOperator : public BasePlotOperator{
public:
    columnDensityOperator(const BaseEmissionLine* l) : BasePlotOperator(l){constant_normalization=BaseEmissionLine::get_boxsize(); name="columnDensity";}
    void increment(int ipix, const BaseParticle* p,const BaseCell* cell, double length, double dx);
};

class MWDensityOperator : public BasePlotOperator{
public:
    MWDensityOperator(const BaseEmissionLine* l) : BasePlotOperator(l){constant_normalization=1.0; name="massWeightedDensity";}
    void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx);
};
class MWDustDensityOperator : public BasePlotOperator{
public:
    MWDustDensityOperator(const BaseEmissionLine* l) : BasePlotOperator(l){constant_normalization=1.0; name="massWeightedDustDensity";}
    void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx);
};
class MWTemperatureOperator : public BasePlotOperator{
public:
    MWTemperatureOperator(const BaseEmissionLine* l) : BasePlotOperator(l){constant_normalization=1.0; name="massWeightedTemperature";}
    void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx);
};
class opticalDepthOperator : public BasePlotOperator{
public:
    opticalDepthOperator(const BaseEmissionLine* l,bool wavelength_set_in_params=true) : BasePlotOperator(l)
    {
        constant_normalization=1.0;
        name="opticalDepth";
        if(wavelength_set_in_params)
        {
            LParmParse pp;
            double wavelength=0.0;
            pp.get("plotter.opticalDepth.wavelength",wavelength);
            //this should be in angstrom relative to the line center.
            this->set_frequency(wavelength);
        }
        
    }
    void set_frequency(double wavelength);//wavelength is offset wrt line center in Angstrom
    void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx);
private:
    double observed_frequency=0.0;
};

class opticalDepthGridOperator : public BasePlotOperator{
public:
    opticalDepthGridOperator(const BaseEmissionLine* l) : BasePlotOperator(l)
    {
        constant_normalization=1.0;
        name="opticalDepthGrid";
        LParmParse pp;
        pp.get("plotter.opticalDepthGrid.nlambda",nlambda);
        //this should be in angstrom relative to the line center.
        pp.get("plotter.opticalDepthGrid.lambda_left",lambda_left);
        pp.get("plotter.opticalDepthGrid.lambda_right",lambda_right);
        dlambda = (lambda_right-lambda_left)/nlambda;
        taus.resize(nlambda);
        for(int i=0;i<nlambda;i++)
            taus[i].resize(npixels*npixels);
    }
    
    void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx);
    void write_results();
private:
    std::vector<std::vector<double> > taus;
    double dlambda=0.0;
    double nlambda=0;
    double lambda_left;
    double lambda_right;
};


class LoSVelocityOperator : public BasePlotOperator{
public:
    LoSVelocityOperator(const BaseEmissionLine* l) : BasePlotOperator(l){
    constant_normalization=1e-5;//to convert cm/s to km/s 
    name="LoSVelocity";
    LParmParse pp;
    pp.get_line_of_sight("plotter.los",los);
    normalize(los.x);
    }
    void increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx);
private:
    line_of_sight los;
    
};
BasePlotOperator* get_plot_operator(std::string name,const BaseEmissionLine *l);
