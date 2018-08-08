#include "Selectors.H"
#include <string>
#include <iostream>

BaseDataset * get_dataset(std::string name)
{
    //translates string name in parameter file to actual dataset object and returns a new instance of the dataset.
    BaseDataset* ds;
    if(name.compare(std::string("SphericalShell"))==0)
        ds = new SphericalShellData;
    else if(name.compare(std::string("InfiniteSlab"))==0)
        ds = new InfiniteSlabData;
    
    else
    {
        std::cout << "Unknown dataset_type" << std::endl;
        std::abort();
        
    }
    return ds;
}

DustModule* get_dust_module(std::string name)
{
    DustModule* dm = NULL;
    if(name.compare(std::string("Null"))==0)
        {
            dm = new DustModuleNull;
        }
        if(name.compare(std::string("Verhamme12"))==0)
        {
        dm = new DustModuleVerhamme12;
        }
        if(name.compare(std::string("Dahlia"))==0)
        {
        dm = new DustModuleDahlia;
        }
        if(dm == NULL)
        {
        std::cout << "get_dust_module(): invalid choice of dust model." << std::endl;
        std::abort();
        }
        return dm;  
}

NeutralFractionModule* get_neutral_fraction_module(std::string name)
{
    LParmParse pp;
    double redshift = 0.0;
    NeutralFractionModule* dm = NULL;
    if(name.compare(std::string("Null"))==0)
        {
            if(Parallel::IOProcessor())
                std::cout << "Neutral Fraction Model: Null" << std::endl;
            dm = new NeutralFractionNull;
        }
        if(name.compare(std::string("Chardin17"))==0)
        {
            if(Parallel::IOProcessor())
                std::cout << "Neutral Fraction Model: Chardin17" << std::endl;
            pp.get("redshift",redshift);
            dm = new NeutralFractionChardin17(redshift);
        }
        if(name.compare(std::string("Rahmati13"))==0)
        {
            if(Parallel::IOProcessor())
                std::cout << "Neutral Fraction Model: Rahmati13" << std::endl;

            pp.get("redshift",redshift);
            dm = new NeutralFractionRahmati13(redshift);
        }
        if(dm == NULL)
        {
        std::cout << "get_neutral_fraction_module(): invalid choice of neutral fraction model." << std::endl;
        std::abort();
        }
        return dm;  
}
