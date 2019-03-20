#include "SlimOutput.H"
#include "LightParmParse/LParmParse.H"
#include <iostream>
#include "LymanAlphaLine.H"
#include "Parallel.H"
#include <cmath>


SlimOutput::~SlimOutput()
{
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    if(escaped_spectrum != NULL)
        delete escaped_spectrum;
    if(destroyed_spectrum != NULL)
        delete destroyed_spectrum;
    if(peeling_spectra != NULL)
    {
        for(int i=0;i<nfiles*nthreads;i++)
            delete peeling_spectra[i];
        delete peeling_spectra;
    }
        
}

void SlimOutput::writePhoton(BaseParticle &p)
{
    double l = line->convert_wavelength(p.frequency);
    escaped_spectrum->add(l,p.weight);
    totalDone++;
    totalScatterings += p.number_of_scatterings;
}
void SlimOutput::writeDeadPhoton(BaseParticle &p)
{
    double l = line->convert_wavelength(p.frequency);
    destroyed_spectrum->add(l,p.weight);
    totalDeadDone++;
    totalScatterings += p.number_of_scatterings; 
}
void SlimOutput::writePeelingOff(BaseParticle &p,int index)
{
    double l = line->convert_wavelength(p.frequency);
    int ithread = 0;
#ifdef _OPENMP
    ithread = omp_get_thread_num();
#endif
    peeling_spectra[ithread*nfiles+index]->add(l,p.weight);
#pragma omp atomic
    contributions=contributions+1;
}
void SlimOutput::open_files()
{
    LParmParse pp;
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    int nbins;
    double ledge,redge;
    pp.get("output.slim.left_edge",ledge);
    pp.get("output.slim.right_edge",redge);
    pp.get("output.slim.nbins",nbins);
    
    destroyed_spectrum = new Histogram(ledge,redge,nbins);
    escaped_spectrum = new Histogram(ledge,redge,nbins);
    if(use_peeling_off)
        peeling_spectra = new Histogram*[nthreads*nfiles];
    if(use_peeling_off)
        {
        
            for(int i=0;i<nthreads;i++)
            {
                    for(int j=0;j<nfiles;j++)
                    {
                        peeling_spectra[i*nfiles+j] = new Histogram(ledge,redge,nbins);
                    }
            }
        }
}

void SlimOutput::open_peeling_files()
{
    //does nothing
}
void SlimOutput::close_files()
{
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif
    //concat the thread-wise spectra; after this operation, the full spectra are in the array entries 0 to nfiles.
    if(use_peeling_off)
    {
        for(int i=1;i<nthreads;i++)
            {
                for(int j=0;j<nfiles;j++)
                    {
                        peeling_spectra[0*nfiles+j]->add_histogram(*peeling_spectra[i*nfiles+j]);
                    }
            }
    }
}
void SlimOutput::cleanup()
{
    //close and combine files.
    //TODO write directly to HDF5.
    int answer=42;
    long scatterings = totalScatterings;
    long nphotons  = totalDone+totalDeadDone;
    long nphotons_alive = totalDone;
    long nphotons_dead = totalDeadDone;
    long totalcontribs = contributions;
    long total_number_of_events = number_of_events;
    
    Parallel::ReduceLongSum(scatterings);
    Parallel::ReduceLongSum(nphotons);
    Parallel::ReduceLongSum(nphotons_alive);
    Parallel::ReduceLongSum(nphotons_dead);
    Parallel::ReduceLongSum(totalcontribs);
    Parallel::ReduceLongSum(total_number_of_events);
    
    escaped_spectrum->reduce();
    destroyed_spectrum->reduce();


    if(Parallel::IOProcessor())
        std::cout << "Total number of escaped photons " << nphotons << " with " << (double)scatterings/nphotons << " scatterings on average. " << nphotons_dead << " photons were absorbed, escape fraction is " << ((double)nphotons_alive)/(nphotons)*100.0 << " %."  << std::endl;
    std::string path(".");
    Parallel::Barrier();
    //concatenate the photon output
    int procs = Parallel::NProcs();
    if(Parallel::IOProcessor())
    {
        if(write_photons)
        {
            
            std::string fnout("Output.txt");
            LParmParse pp;
            pp.query("output.prefix",fnout);
            escaped_spectrum->write(path,fnout);
            destroyed_spectrum->write(path,fnout+"_destroyed");
        }
    }
    //concatenate the peeling off output
    if(Parallel::IOProcessor() && use_peeling_off)
    {
        std::string fnout("peeling.txt");
        LParmParse pp;
        pp.query("output.peeling_off_prefix",fnout);
        
        for(int j=0;j<nfiles;j++)
        {
            std::string fname = fnout+std::string("_")+std::to_string(j);
            peeling_spectra[0*nfiles+j]->write(path,fname);
        }
        
        
    }
    //concatenate the retrace output

    if(Parallel::IOProcessor() && use_retrace)
    {
        LParmParse pp;
        std::string fnout("scatterings");
        pp.query("output.prefix_retrace",fnout);
        std::ofstream out(fnout, std::ios_base::binary);
        
        std::cout << "Wrote " << total_number_of_events << " events for retracing." << std::endl;
        out.write((char*)&total_number_of_events,sizeof(long));
        for(int i=0;i<procs;i++)
        {
            std::string fn = fnout + std::string("_") + std::string(std::to_string(i));
            std::ifstream in(fn, std::ios_base::binary);
            //dumping of an empty file will fail, so we need to check if the file is empty.
            if(in.peek() == std::ifstream::traits_type::eof())
                continue;
            out << in.rdbuf();
            in.close();
        }
        out.write((char*)&answer,sizeof(int));
        
    }
}

void SlimOutput::writeInputPhotons(BaseParticleVector &vp)
{
    
}
