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
    bool output_tau_stats = false;
     pp.query("unigrid.output_tau_stats",output_tau_stats);
     if(output_tau_stats)
         get_optical_depth_stats();
    return 0;
}
#include "LymanAlphaLine.H"
void UnigridDataset::get_optical_depth_stats()
{
    if(Parallel::IOProcessor())
    {
        std::cout << "UnigridDataset::get_optical_depth_stats(): Collecting data..." << std::endl;
        std::vector<double> tau0;
        std::vector<double> taud;
        std::vector<double> atau0;
        LymanAlphaLine line;
        line.setup();
        for(int ix=0;ix<ngrid;ix++)
        {
            for(int iy=0;iy<ngrid;iy++)
            {
                for(int iz=0;iz<ngrid;iz++)
                {
                    BaseCell& mycell = cells[ix*ngrid*ngrid+iy*ngrid+iz];
                    double pathlength = dx;
                    double fraction_in_dust;
                    double dtau0 = line.get_optical_depth_line_center(&mycell,NULL,pathlength,&fraction_in_dust);
                    double dtaud = fraction_in_dust*dtau0;
                    dtau0 -= dtaud;
                    const double local_v_thermal = line.get_local_thermal_velocity(&mycell);
                    const double a_param = 4.7e-4*(12.85/(local_v_thermal/1.e5));
                    
                    tau0.push_back(dtau0);
                    taud.push_back(dtaud);
                    atau0.push_back(a_param*dtau0);
                }
            }
        }
        unsigned long n = tau0.size();
        std::ofstream of("optical_depth_cells.dat") ;
        int max_log_tau = 15;
        long hist[max_log_tau];
        for(int i=0;i<max_log_tau;i++)
            hist[i] = 0;
        long high_values=0;
        for(unsigned long i=0;i<n;i++)
        {
            of << tau0[i] << " " << taud[i] << " " << atau0[i] << std::endl;
            int index = log10(tau0[i]);
            if(index>0)
            {
                if(index>max_log_tau-1)
                    high_values++;   
                else
                    hist[index]++;
            }
        }
        std::cout << "UnigridDataset::get_optical_depth_stats(): " << high_values <<  " cells have tau>1e10." << std::endl;
        for(int i=0;i<max_log_tau;i++)
            std::cout << "log tau = " << i << ": " << hist[i] << std::endl;
    } 
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
    pathlength = dx*len;
    pathlength = std::min(pathlength,max_step);
    pathlength = std::max(pathlength,min_step);
    status = DS_ALL_GOOD;
    return cell;
    
    return NULL;   
}

const int UnigridDataset::get_nearest_face(const double x[3],double *orth_distance, double distance[3]) const
{
    double up[3] = {false,false,false};
    int smallest = 0;   
//figure out the subgrid coordinates.
    int index[3];
    double low[3];
    for(int i=0;i<3;i++)
    {
        index[i] = x[i]/dx;
        low[i] = index[i]*dx;
    }
    
    double subx[3] = {0,0,0};
    for(int i=0;i<3;i++)
    {
        subx[i] = (x[i] - low[i])/dx;  
    }
    for(int i=0;i<3;i++)
    {
        //figure out which surface is closest.
        if( subx[i]>0.5)
        {
             subx[i]-=0.5;
             up[i] = true;
        }
        if(subx[i]<subx[smallest])
            smallest = i;
    }
    double surface_normal[3];
    //set the surface normal vector
    for(int i=0;i<3;i++)
        surface_normal[i] = 0.0;
    surface_normal[smallest] = 1.0;
    if(!up[smallest])
        surface_normal[smallest] *= -1;
    
    //get the distance to the _center_ of that cell surface.
    double cellC[3];
    for(int i=0;i<3;i++)
        cellC[i] = dx*(index[i]+0.5);
    double center[3];
    for(int i=0;i<3;i++)
        center[i]=cellC[i]+subx[smallest]*dx*surface_normal[i];
    for(int i=0;i<3;i++)
        distance[i] = center[i] - x[i];
    *orth_distance = subx[smallest]*dx;
    return 0;
}
