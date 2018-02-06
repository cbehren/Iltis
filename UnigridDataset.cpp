#include "UnigridDataset.H"
#include <cmath>
#include "Utilities.H"
#include "LightParmParse/LParmParse.H"
#include "TraversalLength.H"

#ifdef _OPENMP
#include "omp.h"
#endif
#include "Parallel.H"
int UnigridDataset::setup()
{
    LParmParse pp;
    pp.query("max_step",max_step);
    std::string filename;
    pp.get("unigrid.filename",filename);
    FILE* fp = fopen(filename.c_str(),"r");
    char line[1000];
    fgets(line,1000,fp);
    sscanf(line,"# %d",&ngrid);
    cells = new BaseCell[ngrid*ngrid*ngrid];
    dx = 1./(ngrid);
    for(int ix=0;ix<ngrid;ix++)
        for(int iy=0;iy<ngrid;iy++)
            for(int iz=0;iz<ngrid;iz++)
            {
                fgets(line,1000,fp);
                double density,temperature,velx,vely,velz,dust_density;
                sscanf(line,"%le %le %le %le %le %le\n",&density,&temperature,&velx,&vely,&velz,&dust_density); 
                BaseCell& mycell = cells[ix*ngrid*ngrid+iy*ngrid+iz];
                mycell.density = density;
                mycell.temperature = temperature;
                mycell.velocity[0] = velx;
                mycell.velocity[1] = vely;
                mycell.velocity[2] = velz;
                mycell.dust_density = dust_density;
                
            }
    
    fclose(fp);
    return 0;
}

const BaseCell* UnigridDataset::data_at(const double x[3], const double k[3],double &pathlength,DatasetStatus &status) const 
{
    if(x[0]<0.0 || x[1]<0.0 || x[2]<0.0 || x[0]>=1.0 || x[1]>=1.0 || x[2]>=1.0)
    {
        status = DS_NOT_IN_DOMAIN;
        return NULL;
    }
    int index[3];
    double low[3];
    for(int i=0;i<3;i++)
    {
        index[i] = x[i]/dx;
        low[i] = index[i]*dx;
    }
    BaseCell* cell = &cells[index[0]*ngrid*ngrid+index[1]*ngrid+index[2]];
    
    double subx[3] = {0,0,0};
    for(int i=0;i<3;i++)
    {
        subx[i] = (x[i] - low[i])/dx;  
    }

    int correction=0;

    double len = TraversalLength(subx,k,&correction);
    pathlength = std::max(len,5e-8);
    pathlength *= dx;
    pathlength = std::min(pathlength,max_step);
    
    status = DS_ALL_GOOD;
    return cell;
    
    
    
    
    return NULL;
    
    
}
