#include "../../LymanAlphaLine.H"
#include "../../LightParmParse/LParmParse.H"
#include "../../BaseParticle.H"
#include "../../BaseCell.H"
#include "../../RandomNumbers.H"
#include <string> 
#include "stdio.h"
#include <iostream>
//tests drawing random frequencies from the Neufeld solution for a static, homogeneous sphere.
//tests drawing random exit directions from the Tasitsiomi 2006 CDF. 
int main(int argc,char* args[])
{
    std::string inputsfile("inputs");
    LParmParse::Initialize(argc-2,args+2,inputsfile.c_str());
    RNG::initialize();

    int ntries = 1000000;
    LymanAlphaLine line;
    line.setup();
    const double T = 2e4;
    const double local_v_thermal = 12.85*sqrt(T/1e4)*1000.0*100.0;//cm/s ;
    const double a_param = 4.7e-4*(12.85/(local_v_thermal/1.e5));
    const double tau0 = 1e7;
    FILE *fp = NULL;
    char fname[100];
    
    sprintf(fname,"drawn_frequencies.txt");
    fp = fopen(fname,"w");
    for(int j=0;j<ntries;j++)
    {
        double newfreq = line.draw_neufeld_frequency(a_param,tau0);
        fprintf(fp,"%le\n",newfreq);
    }
    fclose(fp);
    fp = NULL;
    
    sprintf(fname,"drawn_direction.txt");
    fp = fopen(fname,"w");
    for(int j=0;j<ntries;j++)
    {
        double newfreq = line.draw_exit_direction();
        fprintf(fp,"%le\n",newfreq);
    }
    fclose(fp);
    fp = NULL;

    return 0;

}
