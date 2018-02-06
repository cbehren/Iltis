#ifndef __DUST_MODULE_H
#define __DUST_MODULE_H
  #ifdef _OPENMP
#include "omp.h"
  #endif
#include <iostream>
#include "BaseParticle.H"
#include "BaseCell.H"
//This defines a number of implementations for the treatment of dust.
//TODO this is all for the Lyman-alpha line, more or less. Make it more general?


enum DUST_MODULE_TYPE{DM_NULL,DM_VERHAMME12,DM_DAHLIA};
class DustModule{
protected:
  virtual void init() = 0;
  static double G_Greenstein;
  static double albedo;
  static DUST_MODULE_TYPE type;
public:
  virtual ~DustModule() {}
  virtual double scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given=false) = 0;
  static double GreensteinPhase(const double R);
  static double GreensteinProbability(const double mu);
  static void AnisotropicReemission(const double oldk[],double newk[]);
  virtual double ConvertDust(double density, double temperature, double metallicity) = 0;
  virtual double DustConversionFactor() = 0;
  virtual double tauConstFactor() = 0;
  virtual double getAlbedo() = 0;

};

class DustModuleNull : public DustModule{
protected: 
  void init();
public:
  DustModuleNull() {init();}
  double scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given=false);
  double DustConversionFactor();
  double ConvertDust(double density, double temperature, double metallicity) {return 0;}
  //int ConvertDust(MultiFab& fab,const double dx);  
  double tauConstFactor();
  double getAlbedo();
};


class DustModuleVerhamme12 : public DustModule{
protected: 
  void init();
public:
  DustModuleVerhamme12() {init();}
  double scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given=false);
  double DustConversionFactor();
  double ConvertDust(double density, double temperature, double metallicity);
  //int ConvertDust(MultiFab& fab,const double dx);
  double tauConstFactor();
  double getAlbedo();
};


class DustModuleDahlia : public DustModule{
protected: 
  void init();
  double tauConst=-1;
public:
  DustModuleDahlia();
  double scatter(const BaseCell* cell,const BaseParticle* p,double newk[],bool k_and_atom_velocity_given=false);
  double DustConversionFactor();
  double ConvertDust(double density, double temperature, double metallicity);
  //int ConvertDust(MultiFab& fab,const double dx);
  double tauConstFactor();
  double getAlbedo();
};




#endif
