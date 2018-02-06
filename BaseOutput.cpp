#include "BaseOutput.H"
#include "LightParmParse/LParmParse.H"
#include <iostream>
#include "LymanAlphaLine.H"
#include "Parallel.H"
#include <cmath>

//TODO write the output in a smarter way; right now we have the same format string around three times (for the escaped, absorbed, and input photons), which sooner or later will result in errors.

void BaseOutput::setup(BaseEmissionLine* l)
{
    int myproc= Parallel::MyProc();
    //figure out how many files we need to open, and what their names are. Each process will have its own output file.
    LParmParse pp;
    bool use_peeling_off = false;
    pp.query("use_peeling_off",use_peeling_off);
    if(use_peeling_off)
    {
        nfiles = 1;
        pp.query("number_of_instruments",nfiles);
        fname_peeling = "peeling.txt";
        pp.query("output.peeling_off_prefix",fname_peeling);
        fname_peeling += std::string("_")+std::to_string(myproc)+std::string("_");
    }
    fname = std::string("Output.txt");
    pp.query("output.prefix",fname);
    pp.query("output.binary",binary);
    if(binary && Parallel::IOProcessor())
        std::cout << "BaseOutput::setup(): Generate binary output." << std::endl;
    fname_dead = fname+"_destroyed";

    fname = fname+std::string("_")+ std::to_string(myproc);
    fname_dead = fname_dead+std::string("_")+ std::to_string(myproc);
    line = l;
    open_files();
    
    
}
void BaseOutput::writePhoton(BaseParticle &p)
{
    int type = 0;
    //file << p;
    if(binary)
    {
#pragma omp critical
        {
            double l = line->convert_wavelength(p.frequency);
            file.write((char*)&l,sizeof(double));
            file.write((char*)&p.number_of_scatterings,sizeof(int));
            file.write((char*)&type,sizeof(int));
            file.write((char*)&p.x[0],sizeof(double));
            file.write((char*)&p.x[1],sizeof(double));
            file.write((char*)&p.x[2],sizeof(double));
            file.write((char*)&p.k[0],sizeof(double));
            file.write((char*)&p.k[1],sizeof(double));
            file.write((char*)&p.k[2],sizeof(double));
            file.write((char*)&p.id,sizeof(long));
            file.write((char*)&p.emitter,sizeof(int));
            file.write((char*)&p.path_length,sizeof(double));
            file.write((char*)&p.weight,sizeof(double));
            file.write((char*)&p.lsp[0],sizeof(double));
            file.write((char*)&p.lsp[1],sizeof(double));
            file.write((char*)&p.lsp[2],sizeof(double));
              
            totalDone++;
            totalScatterings += p.number_of_scatterings;}
        
    }
    else
    {
#pragma omp critical
        {
        file << line->convert_wavelength(p.frequency) << " " << p.number_of_scatterings << " "<< "0" << " "  << p.x[0] <<  " " << p.x[1] << " " <<  p.x[2] << " " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " " << p.id << " " << p.emitter << " " << p.path_length << " " << p.weight << " " << p.lsp[0] << " " << p.lsp[1] << " "  << p.lsp[2] << " " << p.optical_depth_seen << std::endl;  
        totalDone++;
        totalScatterings += p.number_of_scatterings;}
    }
}
void BaseOutput::writeDeadPhoton(BaseParticle &p)
{
    
        int type = 0;
    //file << p;
    if(binary)
    {
#pragma omp critical
        {
          double l = line->convert_wavelength(p.frequency);
            file_dead.write((char*)&l,sizeof(double));
            file_dead.write((char*)&p.number_of_scatterings,sizeof(int));
            file_dead.write((char*)&type,sizeof(int));
            file_dead.write((char*)&p.x[0],sizeof(double));
            file_dead.write((char*)&p.x[1],sizeof(double));
            file_dead.write((char*)&p.x[2],sizeof(double));
            file_dead.write((char*)&p.k[0],sizeof(double));
            file_dead.write((char*)&p.k[1],sizeof(double));
            file_dead.write((char*)&p.k[2],sizeof(double));
            file_dead.write((char*)&p.id,sizeof(long));
            file_dead.write((char*)&p.emitter,sizeof(int));
            file_dead.write((char*)&p.path_length,sizeof(double));
            file_dead.write((char*)&p.weight,sizeof(double));
            file_dead.write((char*)&p.lsp[0],sizeof(double));
            file_dead.write((char*)&p.lsp[1],sizeof(double));
            file_dead.write((char*)&p.lsp[2],sizeof(double));
              
            
            totalDeadDone++;
            totalScatterings += p.number_of_scatterings;}
        
    }
    else
    {
#pragma omp critical
        {
        file_dead << line->convert_wavelength(p.frequency) << " " << p.number_of_scatterings << " "<< "0" << " "  << p.x[0] <<  " " << p.x[1] << " " <<  p.x[2] << " " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " " << p.id << " " << p.emitter << " " << p.path_length << " " << p.weight << " " << p.lsp[0] << " " << p.lsp[1] << " "  << p.lsp[2] << " " << p.optical_depth_seen << std::endl;  
        totalDeadDone++;
        totalScatterings += p.number_of_scatterings;}
    }

}
void BaseOutput::writePeelingOff(BaseParticle &p,int index)
{
    if(binary)
    {
        pfiles[index].write((char*)p.x,3*sizeof(double));
        double l = line->convert_wavelength(p.frequency);
        double w = p.weight*exp(-p.optical_depth)*p.bias;
        pfiles[index].write((char*)&l,sizeof(double));
        pfiles[index].write((char*)&w,sizeof(double));
        pfiles[index].write((char*)&p.id,sizeof(long));
        pfiles[index].write((char*)&p.emitter,sizeof(int));
        pfiles[index].write((char*)p.lsp,3*sizeof(double));
        pfiles[index].write((char*)&p.optical_depth,sizeof(double));
        pfiles[index].write((char*)&p.number_of_scatterings,sizeof(double));
        
    }
    else
    {
        pfiles[index] << std::scientific << p.x[0] << " " << p.x[1] << " " << p.x[2] << " " << line->convert_wavelength(p.frequency) << " " << p.weight*exp(-p.optical_depth)*p.bias << " " << " "  << p.id << " " << p.emitter << " " << p.lsp[0] << " " << p.lsp[1] << " "  << p.lsp[2] << " " << p.optical_depth << " " << p.number_of_scatterings << "\n";
    }
    contributions++;
    
}
void BaseOutput::open_files()
{
    auto flag = std::ofstream::out;
    if(binary)
        flag = std::ofstream::out & std::ios_base::binary;
    file.open(fname.c_str(),flag);
    file_dead.open(fname_dead.c_str(),flag);
    if(Parallel::IOProcessor() && !binary)
    {
        file << "#photonfile v0.4 \"delta_lambda scatteringcount type posx posy posz kx ky kz id cpu pathlength weight lspx lspy lspz\"\n";
        file_dead << "#photonfile v0.4 \"delta_lambda scatteringcount type posx posy posz kx ky kz id cpu pathlength weight lspx lspy lspz\"\n";
        
    }
    if(nfiles>0)
    {
        pfiles = new std::ofstream[nfiles];
        for(int i=0;i<nfiles;i++)
        {
            std::string myfn(fname_peeling+std::to_string(i));
            pfiles[i].open(myfn.c_str(),flag);
        }
    }
}
void BaseOutput::close_files()
{
    file.close();
    file_dead.close();
    if(nfiles>0)
    {
        for(int i=0;i<nfiles;i++)
        {
            pfiles[i].close();
        }
    }
    
}
void BaseOutput::cleanup()
{
    //close and combine files.
    //TODO write directly to HDF5.
    int answer=42;
    long scatterings = totalScatterings;
    long nphotons  = totalDone+totalDeadDone;
    long nphotons_alive = totalDone;
    long nphotons_dead = totalDeadDone;
    long totalcontribs = contributions;
    Parallel::ReduceLongSum(scatterings);
    Parallel::ReduceLongSum(nphotons);
    Parallel::ReduceLongSum(nphotons_alive);
    Parallel::ReduceLongSum(nphotons_dead);
    Parallel::ReduceLongSum(totalcontribs);

    if(Parallel::IOProcessor())
        std::cout << "Total number of escaped photons " << nphotons << " with " << (double)scatterings/nphotons << " scatterings on average. " << nphotons_dead << " photons were absorbed, escape fraction is " << ((double)nphotons_alive)/(nphotons)*100.0 << " %."  << std::endl;
    close_files();
    Parallel::Barrier();
    //concatenate the photon output
    int procs = Parallel::NProcs();
    if(Parallel::IOProcessor())
    {
        std::string fnout("Output.txt");
        LParmParse pp;
        pp.query("output.prefix",fnout);
        std::ofstream out(fnout, std::ios_base::binary);
        if(binary)
            out.write((char*)&nphotons_alive,sizeof(long));
        for(int i=0;i<procs;i++)
        {
            std::string fn = fnout + std::string("_") + std::string(std::to_string(i));
            std::ifstream in(fn, std::ios_base::binary);
            out << in.rdbuf();
            in.close();
        }
        if(binary)
        {
            out.write((char*)&answer,sizeof(int));
        }
        std::ofstream out2(fnout+"_destroyed", std::ios_base::binary);
        if(binary)
            out2.write((char*)&nphotons_dead,sizeof(long));
        for(int i=0;i<procs;i++)
        {
            std::string fn = fnout + std::string("_destroyed_") + std::string(std::to_string(i));
            std::ifstream in(fn, std::ios_base::binary);
            out2 << in.rdbuf();
            in.close();
        }
        if(binary)
        {
            out2.write((char*)&answer,sizeof(int));
        }
        
    }
    //concatenate the peeling off output
    if(Parallel::IOProcessor() && nfiles>0)
    {
        std::string fnout("peeling.txt");
        LParmParse pp;
        pp.query("output.peeling_off_prefix",fnout);
        
        for(int j=0;j<nfiles;j++)
        {
            
            std::ofstream out(fnout+std::string("_")+std::to_string(j), std::ios_base::binary);
            if(binary)
                out.write((char*)&totalcontribs,sizeof(long));
            for(int i=0;i<procs;i++)
            {
                std::string fn = fnout + std::string("_") + std::to_string(i) + std::string("_") + std::to_string(j);
                std::ifstream in(fn, std::ios_base::binary);
                out << in.rdbuf();
                in.close();
                if(binary && i==procs-1)
                {
                    out.write((char*)&answer,sizeof(int));
                }
            }
        }
        
        
    }
}
void BaseOutput::writeSlicesParallel(const BaseDataset& ds)
{   
    LParmParse pp;
    bool do_slices = false;
    pp.query("output.do_slices",do_slices);
    if(!do_slices)
        return;
    int resolution = 256;
    pp.query("output.slices.resolution",resolution);
    std::vector<double> densities(pow(resolution,2));
    std::vector<double> temperatures(pow(resolution,2));
    std::vector<double> dust_densities(pow(resolution,2));
    double x[3];
    double dummy;
    DatasetStatus status;

    x[2] = 0.51;
        double k[3] = {1,0,0};
        for(int ix=0;ix<resolution;ix++)
        {
            x[0] = (ix+0.5)/((double)resolution);
            for(int iy=0;iy<resolution;iy++)
            {
                x[1] = (iy+0.5)/((double)resolution);
                const BaseCell* cell = ds.data_at(x,k,dummy,status);
                if(status==DS_ALL_GOOD)
                {
                    
                    temperatures[ix*resolution+iy] = cell->temperature;
                    densities[ix*resolution+iy] =  cell->density;
                    dust_densities[ix*resolution+iy] = cell->dust_density;
                }

                
            }
        }
        Parallel::ReduceDoubleSum(temperatures.data(),pow(resolution,2));
        Parallel::ReduceDoubleSum(densities.data(),pow(resolution,2));
        Parallel::ReduceDoubleSum(dust_densities.data(),pow(resolution,2));
        if(Parallel::IOProcessor())//write a basic slice in the middle of the box.
        {
            std::ofstream ftemperatures("temperature_slice.dat");
            std::ofstream fdensities("density_slice.dat");
            std::ofstream fdust_densities("dust_density_slice.dat");
            for(int ix=0;ix<resolution;ix++)
            {
                for(int iy=0;iy<resolution;iy++)
                {
                    ftemperatures <<  temperatures[ix*resolution+iy] << " ";
                    fdensities << densities[ix*resolution+iy] << " ";
                    fdust_densities << dust_densities[ix*resolution+iy] << " ";
                }
                ftemperatures << std::endl;
                fdensities << std::endl;
                fdust_densities << std::endl;
            }
        }
    
    
    
}



void BaseOutput::writeSlices(const BaseDataset& ds)
{
    LParmParse pp;
    bool do_slices = false;
    pp.query("output.do_slices",do_slices);
    std::string dataset_type("SphericalShell");
    pp.query("dataset_type",dataset_type);
    //if the data is distributed, we need to call the parallel version.
    if(dataset_type.compare("RamsesMPI")==0 ||dataset_type.compare("OctetMPI")==0)
    {
        writeSlicesParallel(ds);
        return;
    }
    
    if(do_slices && Parallel::IOProcessor())//write a basic slice in the middle of the box.
    {
        std::cout << "Writing slices..." << std::endl;
        std::ofstream temperatures("temperature_slice.dat");
        std::ofstream densities("density_slice.dat");
        std::ofstream dust_densities("dust_density_slice.dat");

        //TODO this will not work with a distributed domain.
        
        int resolution = 256;
        pp.query("output.slices.resolution",resolution);
        double x[3];
        double dummy;
        DatasetStatus status;
        x[2] = 0.51;
        double k[3] = {1,0,0};
        for(int ix=0;ix<resolution;ix++)
        {
            x[0] = (ix+0.5)/((double)resolution);
            for(int iy=0;iy<resolution;iy++)
            {
                x[1] = (iy+0.5)/((double)resolution);
                const BaseCell* cell = ds.data_at(x,k,dummy,status);
                if(status!=DS_ALL_GOOD)
                {
                    std::cout << "BaseOutput::writeSlices(): bad status!" << std::endl;
                    
                }
                temperatures << cell->temperature << " ";
                densities << cell->density << " ";
                dust_densities << cell->dust_density << " ";

                
            }
            temperatures << std::endl;
            densities << std::endl;
            dust_densities << std::endl;

        }
    }
    
    
}


void BaseOutput::writeInputPhotons(BaseParticleVector &vp)
{
    int answer = 42;
    int type = 0;
    auto flag = std::ofstream::out;
    if(binary)
        flag = std::ofstream::out & std::ios_base::binary;
    int myproc= Parallel::MyProc();
    long totalinput = (signed long)vp.size();
    Parallel::ReduceLongSum(totalinput);
    std::string fname("input_photons.txt_"+std::to_string(myproc));
    std::ofstream f(fname,flag);
    if(Parallel::IOProcessor() && !binary)
    {
        f << "#photonfile v0.4 \"delta_lambda scatteringcount type posx posy posz kx ky kz id cpu pathlength weight lspx lspy lspz\"\n";
    }
    for(unsigned long i=0;i<vp.size();i++)
    {
        BaseParticle &p = vp[i];
        if(binary)
        {
            double l = line->convert_wavelength(p.frequency);
            f.write((char*)&l,sizeof(double));
            f.write((char*)&p.number_of_scatterings,sizeof(int));
            f.write((char*)&type,sizeof(int));
            f.write((char*)&p.x[0],sizeof(double));
            f.write((char*)&p.x[1],sizeof(double));
            f.write((char*)&p.x[2],sizeof(double));
            f.write((char*)&p.k[0],sizeof(double));
            f.write((char*)&p.k[1],sizeof(double));
            f.write((char*)&p.k[2],sizeof(double));
            f.write((char*)&p.id,sizeof(long));
            f.write((char*)&p.emitter,sizeof(int));
            f.write((char*)&p.path_length,sizeof(double));
            f.write((char*)&p.weight,sizeof(double));
            f.write((char*)&p.lsp[0],sizeof(double));
            f.write((char*)&p.lsp[1],sizeof(double));
            f.write((char*)&p.lsp[2],sizeof(double));
            
        }
        else
        {
            f << line->convert_wavelength(p.frequency) << " " << p.number_of_scatterings << " "<< "0" << " "  << p.x[0] <<  " " << p.x[1] << " " <<  p.x[2] << " " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " " << p.id << " " << p.emitter << " " << p.path_length << " " << p.weight << " " << p.lsp[0] << " " << p.lsp[1] << " "  << p.lsp[2] << " " << p.optical_depth_seen << std::endl;  
        }
    }
    
    f.close();
    Parallel::Barrier();
    if(Parallel::IOProcessor())
    {
        std::string fnout("input_photons.txt");
        int procs = Parallel::NProcs();
            
        std::ofstream out(fnout, std::ios_base::binary);
        if(binary)
        {
            out.write((char*)&totalinput,sizeof(long));
        }
        for(int i=0;i<procs;i++)
        {
            std::string fn = "input_photons.txt_" + std::to_string(i);
            std::ifstream in(fn, std::ios_base::binary);
            out << in.rdbuf();
            in.close();
        }
        if(binary)
            out.write((char*)&answer,sizeof(int));
        out.close();
    }
    
}
