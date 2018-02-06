#ifndef __SELECTORS_H
#define __SELECTORS_H
#include "NeutralFractionModule.H"
#include "DustModule.H"
#include "SphericalShellData.H"
#include "InfiniteSlabData.H"
#include "Parallel.H"
#include "LightParmParse/LParmParse.H"

NeutralFractionModule* get_neutral_fraction_module(std::string name);
DustModule* get_dust_module(std::string name);
BaseDataset * get_dataset(std::string name);
#endif
