#include "BaseParticleVector.H"
#include <algorithm>
#include "LightParmParse/LParmParse.H"
#include "Parallel.H"

#define PARTICLE_CHUNK_SIZE 1000000UL
bool BaseParticleVector::distributed_data = false;

bool by_status(const  BaseParticle& p1, const BaseParticle& p2)
{return p1.status > p2.status;}
bool by_order(const  BaseParticle& p1, const BaseParticle& p2)
{return p1.order < p2.order;}
void BaseParticleVector::sort_by_status()
{
    std::sort(begin(),end(),by_status);
}
void BaseParticleVector::sort_by_order()
{
    std::sort(begin(),end(),by_order);
}
void BaseParticleVector::exchange()
{
    Parallel::Barrier();
    if(!distributed_data)
      return;
    //figure out how many particles we have to transfer from which proc to which proc.
    sort_by_order();
    if(Parallel::IOProcessor())
        std::cout << "Entering exchange()" << std::endl;
    int nprocs = Parallel::NProcs();
    int myproc = Parallel::MyProc();
    std::vector<std::vector<BaseParticle> > vps;
    vps.resize(nprocs);
    for(auto it = rbegin(); it != rend();it++)
    {
        vps[it->order].push_back(*it);
        pop_back();        
    }
    std::vector<int> counts_for_procs(nprocs);
    for(int i=0;i<nprocs;i++)
    {
        counts_for_procs[i] = vps[i].size();
//         std::cout << myproc << " has " <<counts_for_procs[i] << "belonging to proc " << i <<  std::endl;
    }
    for(int i=0;i<nprocs;i++)
    {
        bool has_particles=true;
        long offset = 0;
        long size_iproc = vps[i].size();
        Parallel::ReduceLongSum(size_iproc);
        if(i==myproc)
            resize(size_iproc);
        while(has_particles)
        {
            std::vector<int> counts_for_proc_i(nprocs);
//             if(Parallel::IOProcessor())
//                 std::cout << "trying to exchange particle for proc " << i << std::endl;
            counts_for_proc_i[myproc]=std::min(vps[i].size(),PARTICLE_CHUNK_SIZE);
            Parallel::ReduceIntSum(counts_for_proc_i.data(),nprocs);
            int total_for_iproc_this_chunk=0;
            total_for_iproc_this_chunk=counts_for_proc_i[myproc];
            Parallel::ReduceIntSum(total_for_iproc_this_chunk);
//             if(Parallel::IOProcessor())
//             {
//                 std::cout << "total: " << total_for_iproc_this_chunk << " in detail:" << std::endl;
//                 for(int j=0;j<nprocs;j++)
//                     std::cout << counts_for_proc_i[j] << std::endl;
//             }
            
            //get the floating point data.
            std::vector<double> sndfloats(counts_for_proc_i[myproc]*BaseParticle::nfloats);
            std::vector<long> sndlongs(counts_for_proc_i[myproc]*BaseParticle::nlongs);
            long findex=0;
            long lindex=0;
            int iparticle=0;
            for(auto it = vps[i].rbegin(); it != vps[i].rend();it++)
            {
                if(iparticle==PARTICLE_CHUNK_SIZE)
                    break;
                it->aggregate_floats(sndfloats,findex);
                it->aggregate_longs(sndlongs,lindex);
                vps[i].pop_back();
                iparticle++;
                
            }
            double *dbuf=NULL;
            long *lbuf=NULL;
            std::vector<double> rcvfloat;
            std::vector<long> rcvlong;
            if(myproc==i)
            {
                rcvfloat.resize(total_for_iproc_this_chunk*BaseParticle::nfloats);
                rcvlong.resize(total_for_iproc_this_chunk*BaseParticle::nlongs);
                dbuf = rcvfloat.data();
                lbuf = rcvlong.data();
            }
            std::vector<int> fcounts = counts_for_proc_i;
            std::vector<int> lcounts = counts_for_proc_i;
            for(unsigned int j=0;j<fcounts.size();j++)
            {
                fcounts[j]*=BaseParticle::nfloats;
                lcounts[j]*=BaseParticle::nlongs;
            }
            
            
//             std::cout << myproc << ": trying gatherv: " << sndfloats.size() << " " << counts_for_proc_i[myproc] << std::endl;
            
            Parallel::Gatherv(sndfloats.data(),sndfloats.size(),fcounts.data(),dbuf,i);
            Parallel::Gatherv(sndlongs.data(),sndlongs.size(),lcounts.data(),lbuf,i);
            if(myproc==i)
            {
                //unpack
                findex=0;
                lindex=0;
                for(long j=offset;j<offset+total_for_iproc_this_chunk;j++)
                {
                    at(j).read_floats(rcvfloat,findex);
                    at(j).read_longs(rcvlong,lindex);
                }
                
            }
            offset+=total_for_iproc_this_chunk;
            has_particles = !vps[i].empty();
            Parallel::ReduceBoolOr(has_particles);
        }
    }   
    Parallel::Barrier();
    
    
    
}
BaseParticleVector::BaseParticleVector()  : std::vector<BaseParticle>()
{
    LParmParse pp;
    distributed_data = false;
    std::string name;
    pp.query("dataset_type",name);
    if(name.compare("RamsesMPI")==0 || name.compare("OctetMPI")==0)
        distributed_data=true;
    
    
    
}

int BaseParticleVector::check_consistency()
{
    int myproc = Parallel::MyProc();
    long nparticles = (signed long)size();
    
    Parallel::ReduceLongSum(nparticles);
    
    std::cout << "We have " << nparticles << " in total" << std::endl;
    for(auto it = rbegin(); it != rend();it++)
        if(myproc!=it->order)
        {
            std::cout << "something went wrong " << std::endl;
            std::abort();
        }
    return 0;   
}

int BaseParticleVector::print_balance()
{
    int nprocs = Parallel::NProcs();
    int myproc = Parallel::MyProc();
    std::vector<long> nparticles(nprocs);
    nparticles[myproc] = size();
    Parallel::ReduceLongSum(nparticles.data(),nprocs);
    if(Parallel::IOProcessor())
    {
        std::cout << "process id\t# of particles" << std::endl;
        for(int i=0;i<nprocs;i++)
            std::cout << i << "\t" << nparticles[i] << std::endl;
      std::cout << std::endl;   
    }
    return 0;
    
}
