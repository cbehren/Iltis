#include "hdf5_generic.H"
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>     /* atof */
#include <cstdlib>
#include "stdio.h"
#include <ctype.h>
#include "../LightParmParse/LParmParse.H"
#include "../Ramses/RAMSES_info.cpp"
#include "utils.H"

void read_meta(metadata& meta)
{
    
    bool has_shifted = false;
    meta.redshift = 0.0;
    meta.redshift_shifted = 0.0;
    meta.lengthOfBox = -1.0;
    meta.number_of_instruments = 0;
    meta.cutOffLength = -1.0;
    
    
    LParmParse pp;
    pp.query("redshift",meta.redshift);
    if(pp.query("redshift_shifted",meta.redshift_shifted))
        has_shifted = true;
    else
        meta.redshift_shifted = meta.redshift;
    std::string type("SphericalShell");
    pp.query("dataset_type",type);
    if(type.compare("Octet")==0 ||type.compare("SphericalShell")==0)
    {
        pp.get("boxsize",meta.lengthOfBox);
    }
    else if(type.compare("Ramses")==0)
    {
        std::string info_file;
        pp.get("ramses.root",info_file);
        RAMSES::snapshot rsnap(info_file);
        std::cout << rsnap.m_header.unit_l << " " << rsnap.m_header.aexp << std::endl;
        meta.lengthOfBox = rsnap.m_header.unit_l/rsnap.m_header.aexp;
        
    }
    
    
    bool use_peeling_off=false;
    pp.query("use_peeling_off",use_peeling_off);
    if(use_peeling_off)
    {
        meta.number_of_instruments = 1;
        pp.query("number_of_instruments",meta.number_of_instruments);
        for(int i=0;i<meta.number_of_instruments;i++)
        {
            std::vector<double> temp;
            if(!pp.queryktharr("line_of_sight",i,temp))
            {
                std::cout << "Failed to read the line of sight for instrument " << i << std::endl;
                std::abort();
            }
            line_of_sight l;
            for(int j=0;j<3;j++)
            {
                std::cout << temp[j]<< " ";
                l.x[j] = temp[j];
            }
            meta.observationDirections.push_back(l);
            
        }
    }
    
}

void write_meta_to_hdf5(hdf5_file& file,metadata& meta,int i)
{
    file.create_attribute("redshift",meta.redshift);
    file.create_attribute("lengthOfBox",meta.lengthOfBox);
    if(meta.redshift_shifted>0.0)
        file.create_attribute("redshift_shifted",meta.redshift_shifted);
    file.create_attribute("number_of_instruments",meta.number_of_instruments);
    file.create_attribute("cutOffLength",meta.cutOffLength);
    std::string name("observationDirection");
    char arr[3] = {'X','Y','Z'};
    if(meta.number_of_instruments!=0)        
        for(int j=0;j<3;j++)
        {
            std::string myname = name+arr[j];
            char bla[100];
            sprintf(bla,"%s",myname.c_str());
            file.create_attribute(bla,meta.observationDirections[i].x[j]);
            
                
        }
}


