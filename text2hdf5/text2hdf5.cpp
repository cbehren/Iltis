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
//this tool converts the text output of Sedona into a hdf5 file
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
    std::string inputsname(args[1]);
    LParmParse::Initialize(argc-2,args+2,inputsname.c_str());
    metadata meta;
    //read the meta data from the inputs file if desired.
    read_meta(meta);
    std::string inrootname("peeling.txt");
    LParmParse pp;
    pp.query("output.peeling_off_prefix",inrootname);
    bool binary = false;
    pp.query("output.binary",binary);
    char flag[5];
    sprintf(flag,"r");
    if(binary)
        sprintf(flag,"rb");
    
    std::string outnameprefix(inrootname);
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> lspx;
    std::vector<double> lspy;
    std::vector<double> lspz;
    std::vector<double> lambda;
    std::vector<double> weight;
    std::vector<int> pindex;
    std::vector<int> eindex;
    std::vector<double> tau;
    std::vector<int> number_of_scatterings;
    
 
    //we should scan every file to find the total length, rendering allocation less awkward.
    std::vector<long> counts;
    counts.resize(meta.number_of_instruments);
    for(int i=0;i<meta.number_of_instruments;i++)
    {
        char fname[1000];
        sprintf(fname,"%s_%d",inrootname.c_str(),i);
        if(binary)
        {
            FILE *fp = fopen(fname,flag);
            fread((char*)&counts[i],1,sizeof(long),fp);
            fclose(fp);
            
        }
        else
        {
            FILE *fp = fopen(fname,flag);
                    if(!fp)
                    {
                        std::cout << "Could not open file " << fname << std::endl;
                    std::abort();

                    }
            char line[2000];
            int idummy;
            long ldummy;
            double fdummy;
            while(fgets(line,2000,fp))
            {
    //              peelingOffFile << std::scientific << p.m_pos[0] << " " << p.m_pos[1] << " " << p.m_pos[2] << " " << freqdata << " " << weight << " " << " "  << p.m_id << " " << p.m_cpu << " " << p.m_data[INDEX_LSP] << " " << p.m_data[INDEX_LSP+1] << " "  << p.m_data[INDEX_LSP+2] << " " << p.m_data[INDEX_PATHLENGTH] << "\n";
                if(sscanf(line,"%le %le %le %le %le %d %d %le %le %le %le",&fdummy,&fdummy,&fdummy,&fdummy,&fdummy,&idummy,&idummy,&fdummy,&fdummy,&fdummy,&fdummy)!=11)
                {
                    std::cout << "wrong format in line: " << line << std::endl;
                    std::abort();
                }
                counts[i]++;
                
            }
            fclose(fp);
        };
        std::cout << "Found " << counts[i] << " contributions in instrument " <<  i << ".\n";
    }
    

    

    
   
    for(int i=0;i<meta.number_of_instruments;i++)
    {
        long c=0;
        long count = counts[i];
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
        tau.resize(count);
        number_of_scatterings.resize(count);
        char fname[1000];
        sprintf(fname,"%s_%d",inrootname.c_str(),i);
        FILE *fp = fopen(fname,flag);
        char line[2000];
        double fdummy;
        long ldummy;
        if(binary)
        {
            fread((char*)&ldummy,1,sizeof(long),fp);
            for(long j=0;j<count;j++)
            {
                fread((char*)&x[c],1,sizeof(double),fp);
                fread((char*)&y[c],1,sizeof(double),fp);
                fread((char*)&z[c],1,sizeof(double),fp);
                
                fread((char*)&lambda[c],1,sizeof(double),fp);
                fread((char*)&weight[c],1,sizeof(double),fp);
                fread((char*)&pindex[c],1,sizeof(long),fp);
                fread((char*)&eindex[c],1,sizeof(int),fp);
                fread((char*)&lspx[c],1,sizeof(double),fp);
                fread((char*)&lspy[c],1,sizeof(double),fp);
                fread((char*)&lspz[c],1,sizeof(double),fp);
                fread((char*)&tau[c],1,sizeof(double),fp);
                fread((char*)&number_of_scatterings[c],1,sizeof(int),fp);
//                 std::cout << x[c] << " " << y[c] << " " << z[c] << " " << lambda[c] << " " << weight[c] << " " << pindex[c] << " " << eindex[c] << " " << lspx[c] << " " << lspy[c] << " " << lspz[c] << " " << tau[c] << std::endl;
                c++;
            }
            
            
        }
        else
        {
            while(fgets(line,2000,fp))
            {
                if(sscanf(line,"%le %le %le %le %le %d %d %le %le %le %le %d" ,&(x[c]),&(y[c]),&(z[c]),&(lambda[c]),&(weight[c]),&(pindex[c]),&(eindex[c]),&(lspx[c]),&(lspy[c]),&(lspz[c]),&(tau[c]),&(number_of_scatterings[c]))!=12)
                {
                    std::cout << "wrong format in line: " << line << std::endl;
                    std::abort();
                }
                c++;
            }
        }
        int check=0;
        if(binary)
        {
            fread((char*)&check,1,sizeof(int),fp);
            if(check!=42)
            {
                std::cout << "check byte is wrong!" << std::endl;
                std::abort();
            }
            
        }
        fclose(fp);
        
        char temp[100];
        sprintf(temp,"%s_%d.hdf5",inrootname.c_str(),i);
        hdf5_file hf(temp,'o');
        hf.create_dataset("x",count,x.data());
        hf.create_dataset("y",count,y.data());
        hf.create_dataset("z",count,z.data());
        hf.create_dataset("lambda",count,lambda.data());
        hf.create_dataset("weight",count,weight.data());
        hf.create_dataset("pindex",count,pindex.data());
        hf.create_dataset("eindex",count,eindex.data());
        hf.create_dataset("lspx",count,lspx.data());
        hf.create_dataset("lspy",count,lspy.data());
        hf.create_dataset("lspz",count,lspz.data());
        hf.create_dataset("tau",count,tau.data());
        hf.create_dataset("scount",count,number_of_scatterings.data());
        
        write_meta_to_hdf5(hf,meta,i);
        c=0;
            
        
    }
    std::cout << "Done.\n";
    return 0;
}


