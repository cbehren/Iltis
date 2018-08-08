#include "PlottingOperators.H"
#include "stdio.h"
#include <iostream>
#include <cmath>

const BaseEmissionLine* BasePlotOperator::line=NULL;

BasePlotOperator* get_plot_operator(std::string name,const BaseEmissionLine* l)
{
    if(l==NULL)
    {
        std::cout << "get_plot_operator: Must provide EmissionLine object!" << std::endl;
        std::abort();
    }
    if(name.compare(std::string("columnDensity"))==0)
        return new columnDensityOperator(l);
    if(name.compare(std::string("massWeightedDensity"))==0)
        return new MWDensityOperator(l);
    if(name.compare(std::string("massWeightedDustDensity"))==0)
        return new MWDustDensityOperator(l);
    if(name.compare(std::string("massWeightedTemperature"))==0)
        return new MWTemperatureOperator(l);
    if(name.compare(std::string("LoSVelocity"))==0)
        return new LoSVelocityOperator(l);
    if(name.compare(std::string("opticalDepth"))==0)
        return new opticalDepthOperator(l);
    if(name.compare(std::string("opticalDepthGrid"))==0)
        return new opticalDepthGridOperator(l);
    if(name.compare(std::string("slice"))==0)
        return new sliceOperator(l);
    std::cout << "Unknown plot operator with name " << name << std::endl;
    std::abort();
}

BasePlotOperator::BasePlotOperator(const BaseEmissionLine* l)
{
    LParmParse pp;
    pp.get("plotter.npixels",npixels);
    contributions.resize(npixels*npixels);
    normalization.resize(npixels*npixels);
    line = l;
}

void BasePlotOperator::write_results()
{
    LParmParse pp;
    line_of_sight los,center,up;
    double width,depth;
    std::string output_prefix;
    pp.get_line_of_sight("plotter.los",los);
    pp.get_line_of_sight("plotter.center",center);
    pp.get_line_of_sight("plotter.up",up);
    pp.get("plotter.width",width);
    pp.get("plotter.depth",depth);
    pp.get("plotter.output_prefix",output_prefix);

    std::string fname;
    
    //put data into array.
  
    fname = output_prefix + name;
    FILE* fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
    
    for(unsigned int j=0;j<contributions.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", get_value_normed(j));
    }
    fclose(fp);
}
double BasePlotOperator::get_value_normed(int ipix)
{
    return contributions[ipix]/normalization[ipix]*constant_normalization;
}

void columnDensityOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    contributions[ipix] += length*cell->density;
    //std::cout << p->x[0] << " " << p->x[1] << p->x[2] << " "  << length << " density " << cell->density << std::endl;
    normalization[ipix] = 1.0;
}

void opticalDepthOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    double fraction_in_dust=0.0;
    BaseParticle p2 = *p;
    p2.frequency = observed_frequency;

    double dtau = line->get_optical_depth(cell,&p2,length,&fraction_in_dust);
    contributions[ipix] += dtau;
    normalization[ipix] = 1.0;
}

void opticalDepthGridOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    double lambda0 = Constants::CGS::cval/line->Nu0();//cm
    for(int i=0;i<nlambda;i++)
    {
        double w = (lambda_left+i*dlambda)/1e8 + lambda0;
        //std::cout << i << "wavelength " << w-lambda0 << std::endl;
        double observed_frequency = Constants::CGS::cval/w;
        double fraction_in_dust=0.0;
        BaseParticle p2 = *p;
        p2.frequency = observed_frequency-line->Nu0();
        //std::cout << "lambda " << line->convert_wavelength(p2.frequency) << std::endl; 
        double dtau = line->get_optical_depth(cell,&p2,length,&fraction_in_dust);
        if(divide_by_density)
            dtau/= cell->density;
        taus[i][ipix] += dtau;
        
    }
}

void opticalDepthGridOperator::write_results()
{
    LParmParse pp;
    line_of_sight los,center,up;
    double width,depth;
    std::string output_prefix;
    pp.get_line_of_sight("plotter.los",los);
    pp.get_line_of_sight("plotter.center",center);
    pp.get_line_of_sight("plotter.up",up);
    pp.get("plotter.width",width);
    pp.get("plotter.depth",depth);
    pp.get("plotter.output_prefix",output_prefix);

    std::string fname;
    for(int i=0;i<nlambda;i++)
    {
        char slambda[10];
        sprintf(slambda,"_%d",i);
        //put data into array.
    
        fname = output_prefix + name+ std::string(slambda);
        FILE* fp = fopen(fname.c_str(),"w");
        fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
        fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
        fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
        fprintf(fp,"#width %le\n",width);
        fprintf(fp,"#depth %le\n",depth);
        fprintf(fp,"#lambda %le\n",lambda_left+i*dlambda);
        
        for(unsigned int j=0;j<taus[i].size();j++)
        {
            if(j % npixels == 0)
                fprintf(fp,"\n");
            fprintf(fp,"%le ", taus[i][j]);
        }
        fclose(fp);
    }
}

void opticalDepthOperator::set_frequency(double wavelength)
{
    double lambda0 = Constants::CGS::cval/line->Nu0();//cm
    double w = wavelength/1.0e8 + lambda0;
    observed_frequency = Constants::CGS::cval/w;
}

void MWDensityOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    double mass = pow(dx,3.0)*cell->density;
    contributions[ipix] += cell->density*mass;
    normalization[ipix] += mass;
}
void MWDustDensityOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    double mass = pow(dx,3.0)*cell->dust_density;
    contributions[ipix] += cell->dust_density*mass;
    normalization[ipix] += mass;
}
void MWTemperatureOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    double mass = pow(dx,3.0)*cell->density;
    contributions[ipix] += cell->temperature*mass;
    normalization[ipix] += mass;
}
void LoSVelocityOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    double losvel = scalar(los.x,cell->velocity)*length*cell->density;
    contributions[ipix] += losvel;
    normalization[ipix] += length*cell->density;
}

void sliceOperator::increment(int ipix, const BaseParticle* p, const BaseCell* cell, double length, double dx)
{
    if(cells[ipix]==NULL)
    {
        cells[ipix]=cell;
        dxs[ipix]=dx;
    }
}


void sliceOperator::write_results()
{
    LParmParse pp;
    line_of_sight los,center,up;
    double width,depth;
    std::string output_prefix;
    pp.get_line_of_sight("plotter.los",los);
    pp.get_line_of_sight("plotter.center",center);
    pp.get_line_of_sight("plotter.up",up);
    pp.get("plotter.width",width);
    pp.get("plotter.depth",depth);
    pp.get("plotter.output_prefix",output_prefix);
    
    std::string var("_density");
    std::string fname = output_prefix + name + var;
    FILE* fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", cells[j]->density);
    }
    fclose(fp);
    
    var ="_temperature";
    fname = output_prefix + name + var;
    fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", cells[j]->temperature);
    }
    fclose(fp);
    
    var ="_velocity_x";
    fname = output_prefix + name + var;
    fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", cells[j]->velocity[0]);
    }
    fclose(fp);
    
    var ="_velocity_y";
    fname = output_prefix + name + var;
    fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", cells[j]->velocity[1]);
    }
    fclose(fp);
    
    var ="_velocity_z";
    fname = output_prefix + name + var;
    fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", cells[j]->velocity[2]);
    }
    fclose(fp);
    
    
    var ="_dx";
    fname = output_prefix + name + var;
    fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", dxs[j]);
    }
    fclose(fp);
    
    
    var ="_dust_density";
    fname = output_prefix + name + var;
    fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
    fprintf(fp,"#width %le\n",width);
    fprintf(fp,"#depth %le\n",depth);
            
    for(unsigned int j=0;j<cells.size();j++)
    {
        if(j % npixels == 0)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", cells[j]->dust_density);
    }
    fclose(fp);

    
}
