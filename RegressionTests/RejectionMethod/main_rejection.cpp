#include "../../LymanAlphaLine.H"
#include "../../LightParmParse/LParmParse.H"
#include "../../BaseParticle.H"
#include "../../BaseCell.H"
#include "../../RandomNumbers.H"
#include "../../RejectionMethod.H"
#include <string> 
#include "stdio.h"
#include <iostream>
//tests the rejection method for different values of x (the frequency of the photon)
int main(int argc,char* args[])
{
#ifdef USE_ACCELERATION
    std::cout << "you need to switch off USE_ACCELERATION in LymanAlphaLine.H!";
    std::abort();
#endif
    std::string inputsfile("inputs");
    LParmParse::Initialize(argc-2,args+2,inputsfile.c_str());
    RNG::initialize();

    int ntries = 500000;
    int nins = 6;
    double xins[] = {0.0,1.0,2.0,3.0,4.0,5.0};

    double temperature = 2.74;
    double local_v_thermal = 12.85*sqrt(temperature/1e4); //thermal velocity in km/s
    double a=4.7e-4*(12.85/(local_v_thermal));
    
    char fname[100];
    
    for(int i=0;i<nins;i++)
    {
        double xin = xins[i];
        sprintf(fname,"rejection_method_x=%1.1f.txt",xin);
        FILE *fp = NULL;
        fp = fopen(fname,"w");
        fprintf(fp,"#for xin=%le\n",xin);
        for(int j=0;j<ntries;j++)
        {
            double u = rejection_method(xin,a);
            fprintf(fp,"%le\n",u);
        }
        fclose(fp);
        fp = NULL;
    }

    return 0;

}

