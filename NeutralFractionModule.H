#ifndef __NEUTRALFRACTION_MODULE_H
#define __NEUTRALFRACTION_MODULE_H
  #ifdef _OPENMP
#include "omp.h"
  #endif
#include <iostream>
#include "stdlib.h"
#include "stdio.h"
#include "Interpolate.H"
#include <cstdlib>

//This file declares classes for the calculation of the neutral hydrogen fraction using different models. 

enum NEUTRALFRACTION_MODULE_TYPE{NM_NULL,NM_FIELD,NM_CIE,NM_RAHMATI13,NM_CHARDIN17};
class NeutralFractionModule{
protected:
    
    //INPUT: density in number density (cgs), temperature in K
    //OUTPUT: neutral fraction 
  virtual void init() = 0;
  NEUTRALFRACTION_MODULE_TYPE type;
public:
  virtual double neutralFraction(double density, double temperature) = 0;
  virtual ~NeutralFractionModule() {};

};
//null implementation. Does nothing and fails.
class NeutralFractionNull : public NeutralFractionModule{
protected: 
    void init()
    {
        type = NM_NULL;
        std::cout << "NeutralFractionNull::init(): Invalid selection for NeutralFractionModule\n"; 
        std::abort();
    }
public:  
    NeutralFractionNull(){init();}
    ~NeutralFractionNull() {}
    double neutralFraction(double density, double temperature){init(); return 0;}
};
//This is a stub. At some point we want to implement it for a case in which the neutral fraction is explicitly given by one of the field in the multifab.
class NeutralFractionField : public NeutralFractionModule{
protected: 
    void init()
    {
        type = NM_FIELD;
        std::cout << "NeutralFractionField::init(): Not implemented!\n"; 
        std::abort();
    }
public:  
    NeutralFractionField(){init();}
    double neutralFraction(double density, double temperature){init(); return 0;}
};


//calculate the neutral fraction with the model used in Rahmati 2013
// http://adsabs.harvard.edu/abs/2013MNRAS.430.2427R
// includes CIE+UVB from H+M2012
class NeutralFractionRahmati13 : public NeutralFractionModule{
protected: 
    void init()
    {
        type = NM_RAHMATI13;
 
    }
    double redshift;
    int init_read();
    double alpha_A(double T);
    double Lambda_T(double T);
    
    
    double sigma_vh(double redshift);
    
    Interpolate* interp_gamma = NULL;
    Interpolate* interp_sigma = NULL;

public:  
    NeutralFractionRahmati13(double z);
    ~NeutralFractionRahmati13();
    double neutralFraction(double density, double temperature);

    virtual double Gamma_Phot(double redshift,double nh,double nh_ssh);
    virtual double nH_ssh(double T, double redshift);
    double Gamma12(double redshift);

};


//calculate the neutral fraction with the model used in Chardin 2013
// formalism is the same as in Rahmati, parameters are different
// http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1707.06993
// includes CIE+UVB
class NeutralFractionChardin17 : public NeutralFractionRahmati13{
protected: 
    void init()
    {
        type = NM_CHARDIN17;
    }
    double redshift;
    double init_Chardin17();
    Interpolate* interp_n0 = NULL;
    Interpolate* interp_alpha1 = NULL;
    Interpolate* interp_alpha2 = NULL;
    Interpolate* interp_beta = NULL;
    Interpolate* interp_f = NULL;
    

public:  
    NeutralFractionChardin17(double z);
    ~NeutralFractionChardin17();
    int init_read();
    double Gamma_Phot(double redshift,double nh,double nh_ssh);
    double nH_ssh(double T, double redshift);
     
};


#endif
