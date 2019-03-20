#include "BasePlotOperator.H"
#include "Parallel.H"

const BaseEmissionLine* BasePlotOperator::line=NULL;
bool BasePlotOperator::is_point_sampling = false;

int BasePlotOperator::override_tot_pixels = 0;

BasePlotOperator::BasePlotOperator(const BaseEmissionLine* l)
{
    LParmParse pp;
    pp.get("plotter.npixels",npixels);
    pp.query("plotter.isPointSampling",is_point_sampling);
    int tot_pixels = npixels*npixels;
    if(override_tot_pixels>0)
        tot_pixels = override_tot_pixels;       
    contributions.resize(tot_pixels);
    normalization.resize(tot_pixels);
    line = l;
}

void BasePlotOperator::exchange()
{
    if(has_only_constant_normalization)
        for(unsigned long i=0;i<normalization.size();i++)
            normalization[i]=1.0;
    else
        Parallel::ReduceDoubleSum(normalization.data(),normalization.size());  
    Parallel::ReduceDoubleSum(contributions.data(),contributions.size());
}

void BasePlotOperator::write_results(int index)
{
    LParmParse pp;
    line_of_sight los,center,up;
    double width,depth;
    std::string output_prefix;
    pp.get_line_of_sight("plotter.los",los,index);
    pp.get_line_of_sight("plotter.up",up,index);
    pp.get_line_of_sight("plotter.center",center);
    pp.get("plotter.width",width);
    pp.get("plotter.depth",depth);
    pp.get("plotter.output_prefix",output_prefix);

    std::string fname;
    
    //put data into array.
  
    fname = output_prefix + name + std::string("_") + std::to_string(index);
    FILE* fp = fopen(fname.c_str(),"w");
    fprintf(fp,"#los %le %le %le\n",los.x[0],los.x[1],los.x[2]);
    fprintf(fp,"#up %le %le %le\n",up.x[0],up.x[1],up.x[2]);
    if(!is_point_sampling)
    {
        fprintf(fp,"#center %le %le %le\n",center.x[0],center.x[1],center.x[2]);
        fprintf(fp,"#width %le\n",width);
    }
    fprintf(fp,"#depth %le\n",depth);
    
    for(unsigned int j=0;j<contributions.size();j++)
    {
        if(j % npixels == 0 && !is_point_sampling)
            fprintf(fp,"\n");
        fprintf(fp,"%le ", get_value_normed(j));
    }
    fclose(fp);
}
double BasePlotOperator::get_value_normed(int ipix)
{
    if(contributions[ipix] != contributions[ipix])
        std::cout << "contrib nan" << std::endl;
    if(normalization[ipix] != normalization[ipix])
        std::cout << "norm nan" << std::endl;
    double v = contributions[ipix]/normalization[ipix]*constant_normalization;
//     if(v!=v)
//     {
//         std::cout << "contrib/norm*const nan" << std::endl;
//         std::cout << contributions[ipix] << " " << normalization[ipix] << " " << constant_normalization << std::endl;
//     }
    return v;
}
