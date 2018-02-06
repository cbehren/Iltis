#include "hdf5_generic.H"
#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>     /* atof */
#include <cstdlib>
#include "stdio.h"
#include <ctype.h>
#include "../LightParmParse/LParmParse.H"
#include "utils.H"
//this tool converts the text output of Sedona into a hdf5 file for the NON PEELING OFF output, i.e. the files named "Output.txt"
//TODO: Add metadata for the simulation run.


void read_meta(metadata& meta);
void write_meta_to_hdf5(hdf5_file& file,metadata& meta,int i);





int main(int argc,char* args[])
{
    
    if(argc!=2)
    {
        std::cout << "Usage: " << args[0] << " <name of inputs file>\n";
        return 0;
        
    }
    std::string inrootname("Output.txt");
    std::string inputsname(args[1]);
    LParmParse::Initialize(argc-2,args+2,inputsname.c_str());
    
    LParmParse pp;
    pp.query("output.prefix",inrootname);
    std::string outname(inrootname+std::string(".hdf5"));
    std::string outname_input(std::string("input_photons.hdf5"));
    std::string outname_destroyed(inrootname+std::string("_destroyed.hdf5"));
    bool binary = false;
    pp.query("output.binary",binary);
    char flag[5];
    sprintf(flag,"r");
    if(binary)
        sprintf(flag,"rb");
        


    
    metadata meta;
    
    //read the meta data from the inputs file if desired.
    read_meta(meta);
    //we should scan every file to find the total length, rendering allocation less awkward.
    

    char fname[1000];
    char fname_input[1000];
    char fname_destroyed[1000];
    sprintf(fname,"%s",inrootname.c_str());
    sprintf(fname_input,"%s","input_photons.txt");
    sprintf(fname_destroyed,"%s_destroyed",inrootname.c_str());
    char* names[3] = {fname,fname_input,fname_destroyed};
    std::string outnames[3] = {outname,outname_input,outname_destroyed};
    for(int ifile=0;ifile<3;ifile++)
    {
        long count = 0;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        std::vector<double> lspx;
        std::vector<double> lspy;
        std::vector<double> lspz;
        std::vector<double> lambda;
        std::vector<double> weight;
        std::vector<long> pindex;
        std::vector<int> eindex;
        std::vector<double> plength;
        
        std::vector<double> kx;
        std::vector<double> ky;
        std::vector<double> kz;
        std::vector<int> scount;
        std::vector<int> type;
        char line[2000];
        FILE *fp = fopen(names[ifile],flag);
        if(!fp)
        {
                std::cout << "Could not open file " << names[ifile] << std::endl;
                std::abort();

        }
        if(binary)
        {
            fread((char*)&count,1,sizeof(long),fp);
//             std::cout << "read " << count << std::endl;
            
            
        }
        else
        {
            
            double fdummy;
            fgets(line,2000,fp);
            while(fgets(line,2000,fp))
            {
        //	  #delta_lambda scatteringcount type posx posy posz kx ky kz id cpu pathlength weight lspx lspy lspz"
            if(line[0]=='#')
                break;

            if(sscanf(line,"%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le",&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&fdummy)!=16)
                {
                    std::cout << "wrong format in line: " << line << std::endl;
                    std::abort();
                }
                count++;
                
            }
        }
        fclose(fp);
    
        std::cout << "Found " << count << " photons.\n";

        x.resize(count);
        y.resize(count);
        z.resize(count);
        lspx.resize(count);
        lspy.resize(count);
        lspz.resize(count);
        lambda.resize(count);
        weight.resize(count);
        pindex.resize(count);
        eindex.resize(count);
        plength.resize(count);
        kx.resize(count);
        ky.resize(count);
        kz.resize(count);
        scount.resize(count);
        type.resize(count);

        
        long c=0;

       
        fp = fopen(names[ifile],flag);
        if(binary)
        {
            long dummy=0;
            fread((char*)&dummy,1,sizeof(long),fp);
                
            for(long j=0;j<count;j++)
            {
                fread((char*)&lambda[j],1,sizeof(double),fp);
                fread((char*)&scount[j],1,sizeof(int),fp);
                fread((char*)&type[j],1,sizeof(int),fp);
                fread((char*)&x[j],1,sizeof(double),fp);
                fread((char*)&y[j],1,sizeof(double),fp);
                fread((char*)&z[j],1,sizeof(double),fp);
                fread((char*)&kx[j],1,sizeof(double),fp);
                fread((char*)&ky[j],1,sizeof(double),fp);
                fread((char*)&kz[j],1,sizeof(double),fp);
                fread((char*)&pindex[j],1,sizeof(long),fp);
                fread((char*)&eindex[j],1,sizeof(int),fp);
                fread((char*)&plength[j],1,sizeof(double),fp);
                fread((char*)&weight[j],1,sizeof(double),fp);
                fread((char*)&lspx[j],1,sizeof(double),fp);
                fread((char*)&lspy[j],1,sizeof(double),fp);
                fread((char*)&lspz[j],1,sizeof(double),fp);
                pindex[j]*=-1;
                c++;
            }
            
        }
        else
        {
            fgets(line,2000,fp);
            while(fgets(line,2000,fp))
            {
                if (line[0]=='#')
                break;
                if(sscanf(line,"%le %d %d %le %le %le %le %le %le %ld %d %le %le %le %le %le",&lambda[c],&scount[c],&type[c],&x[c],&y[c],&z[c],&kx[c],&ky[c],&kz[c],&pindex[c],&eindex[c],&plength[c],&weight[c],&lspx[c],&lspy[c],&lspz[c])!=16)
                {
                    std::cout << "wrong format in line: " << line << std::endl;
                    std::abort();
                }
                //fix the pindex - it is negative in the output file.
                pindex[c]*=-1;
                c++;
            }
        }
        fclose(fp);

        char temp[100];
        sprintf(temp,"%s",outnames[ifile].c_str());
        hdf5_file hf(temp,'o');
        hf.create_dataset("x",count,x.data());
        hf.create_dataset("y",count,y.data());
        hf.create_dataset("z",count,z.data());
        hf.create_dataset("kx",count,kx.data());
        hf.create_dataset("ky",count,ky.data());
        hf.create_dataset("kz",count,kz.data());
        hf.create_dataset("lambda",count,lambda.data());
        hf.create_dataset("weight",count,weight.data());
        hf.create_dataset("lspx",count,lspx.data());
        hf.create_dataset("lspy",count,lspy.data());
        hf.create_dataset("lspz",count,lspz.data());
        hf.create_dataset("plength",count,plength.data());
        
        std::vector<int> typeint(type.begin(),type.end());
        hf.create_dataset("type",count,typeint.data());
        std::vector<int> scountint(scount.begin(),scount.end());
        hf.create_dataset("scount",count,scountint.data());
        std::vector<int> pindexint(pindex.begin(),pindex.end());
        hf.create_dataset("pindex",count,pindexint.data());
        std::vector<int> eindexint(eindex.begin(),eindex.end());
        hf.create_dataset("eindex",count,eindexint.data());
        
        
        write_meta_to_hdf5(hf,meta,0);
    }
    std::cout << "Done.\n";
    return 0;
}


