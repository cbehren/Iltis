#include "../../LymanAlphaLine.H"
#include "../../LightParmParse/LParmParse.H"
#include "../../BaseParticle.H"
#include "../../BaseCell.H"
#include "../../RandomNumbers.H"
#include <string> 
#include "stdio.h"
#include <iostream>

int main(int argc,char* args[])
{
#ifdef USE_ACCELERATION
    std::cout << "you need to switch off USE_ACCELERATION in LymanAlphaLine.H!";
    std::abort();
#endif
    std::string inputsfile(args[1]);
    LParmParse::Initialize(argc-2,args+2,inputsfile.c_str());
    RNG::initialize();

    int ntries = 100000;

    int nin = 6;
    double xin[] = {0,1,2,3,4,5};

    double temperature = 2.74;

    BaseCell cell;
    cell.temperature = temperature;
    cell.velocity[0] = 0.0;
    cell.velocity[1] = 0.0;
    cell.velocity[2] = 0.0;

    BaseParticle p;
    p.k[0] = 1.0;
    p.k[1] = 0;
    p.k[2] = 0;
    
    LymanAlphaLine line;
    line.setup();
    FILE *fp = NULL;
    char fname[200];
    for(int i=0;i<nin;i++)
    {
        double freq = LymanAlphaLine::x_to_frequency(xin[i],BaseEmissionLine::get_local_thermal_velocity(&cell));
        p.frequency = freq;
        double bla[3];
        sprintf(fname,"x_%d.txt",i);
        fp = fopen(fname,"w");
        for(int j=0;j<ntries;j++)
        {
            double newfreq = line.scatter(&cell,&p,bla,1e9,bla,false);
            newfreq = LymanAlphaLine::frequency_to_x(newfreq,BaseEmissionLine::get_local_thermal_velocity(&cell));
            fprintf(fp,"%le\n",newfreq);
        }
        fclose(fp);
        fp = NULL;

        

    }

    return 0;

}
