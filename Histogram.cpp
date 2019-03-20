#ifndef __HISTOGRAM_H_
#define __HISTOGRAM_H_
#include "Histogram.H"
#include <iostream>
#include <fstream>
#include <cmath>
#include "Parallel.H"

Histogram::Histogram(double left, double right, double nbin)
{
    left_edge = left;
    right_edge = right;
    nbins = nbin;
    dx = (right_edge-left_edge)/nbins;
    
    values = new double[nbins];
    for(int i=0;i<nbins;i++)
    {
        values[i]=0.0;
    }
    ignored = 0.0;
    
}
Histogram::~Histogram()
{
    if(values!=NULL)
    {
        delete[] values;
        values=NULL;
    }
   
}


void Histogram::add_histogram(const Histogram& hist)
{
    if(fabs(hist.left_edge-left_edge)>1e-9)
    {
        std::cout <<  "add_histogram(): left edge does not match!" << std::endl;
        std::abort();
    }
    if(fabs(hist.right_edge-right_edge)>1e-9)
    {
        std::cout <<  "add_histogram(): right edge does not match!" << std::endl;
        std::abort();
    }
    if(hist.nbins!=nbins)
    {
        std::cout <<  "add_histogram(): nbins does not match!" << std::endl;
        std::abort();
    }    
    for(int i=0;i<nbins;i++)
    {
        values[i]+=hist.values[i];
        ignored+=hist.ignored;
    }
 
}
    
void Histogram::reduce()
{
    Parallel::ReduceDoubleSum(values,nbins);

}
int Histogram::write(std::string directory,std::string filename)
{
    
    std::ofstream of(directory+"/"+filename);
    //write header
    of << "#nbins " << nbins << std::endl;
    of << "#left_edge " << left_edge << std::endl;
    of << "#right_edge " << right_edge << std::endl;
    of << "#dx " << dx << std::endl;
    for(int i=0;i<nbins;i++)
    {
        of << left_edge+(i+0.5)*dx << " " << values[i] << std::endl;
    }
    
    return 0;
}
void Histogram::get_data(std::vector<double>& bin_centers,std::vector<double>& data)
{
    for(int i=0;i<nbins;i++)
    {
        bin_centers[i]=left_edge+(i+0.5)*dx;
        data[i]=values[i];
    }

}




#endif
