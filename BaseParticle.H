#ifndef __BASEPARTICLE_
#define __BASEPARTICLE_

//this defines the default class of particles.

#include <iostream>
#include <fstream>
#include <vector>
//enum variable for the status of a particle. DONE = the particle has escaped, DEAD = the particle has been absorbed, OKAY = particle is in transfer, BP_UNINITIALIZED = particle has not be initialized correctly yet.
enum ParticleStatus {BP_DONE,BP_DEAD,BP_OKAY,BP_UNINITIALIZED};
class BaseParticle{
public:
    BaseParticle();
    //position in code units
    double x[3];
    //last point of scattering (LSP) in code units
    double lsp[3];
    //statistical weight of this particle
    double weight;
    //unit vector in direction of traveling
    double k[3];
    //not used yet (supposed to speed up code by grouping particles spatially)
    int order;
    //index of the emitter from which this particle originates
    int emitter;
    //index of the photon - unique for emitter and process TODO: Probably smarter to make them unique per emitter.
    long id;
    //the optical depth this photon will penetrate before the next interaction occurs
    double optical_depth;
    //the total optical depth it has seen
    double optical_depth_seen;
    //the total path length this photon has travelled in code units
    double path_length;
    //its physical frequency in the observers frame
    double frequency;
    //bias - experimental right now.
    double bias;
    //number of scatterings this particle has seen
    int number_of_scatterings;
    //particle status, see above
    ParticleStatus status;
    //helper for pretty printing
    friend std::ostream& operator<<(std::ostream& stream,const BaseParticle& p);
    
    //inline function for quickly turning particles into array of floats/long and back
    inline void aggregate_floats(std::vector<double>& v,long &index){
        for(int i=0;i<3;i++)
        {
            v[index+i] = x[i];
            v[index+3+i] = k[i];
            v[index+6+i] = lsp[i];
        }
        v[index+9]=weight;
        v[index+10]=optical_depth;
        v[index+11]=optical_depth_seen;
        v[index+12]=path_length;
        v[index+13]=frequency;
        index+=nfloats;
    }
    inline void read_floats(std::vector<double>& v,long &index){
        for(int i=0;i<3;i++)
        {
            x[i] = v[index+i];
            k[i]= v[index+3+i];
            lsp[i] = v[index+6+i];
        }
        weight=v[index+9];
        optical_depth=v[index+10];
        optical_depth_seen=v[index+11];
        path_length=v[index+12];
        frequency=v[index+13];
        index+=nfloats;
    }
    
    inline void aggregate_longs(std::vector<long>& v,long &index){
        v[index]=order;
        v[index+1]=emitter;
        v[index+2]=id;
        v[index+3]=number_of_scatterings;
        index+=nlongs;
    }
    inline void read_longs(std::vector<long>& v,long &index){
        order=v[index];
        emitter=v[index+1];
        id=v[index+2];
        number_of_scatterings=v[index+3];
        index+=nlongs;
    }
    //move a particle by pathlengh (in code units) in direction of k
    void move(double pathlength);
    //move a particle to position and update the path accordingly
    void move_to(const double pos[3]);
    //stores the number of floats/ints that a particle holds.
    static const int nfloats;
    static const int nlongs;
protected:
    //stuff to generate unique IDs
    static long n;
    static long nextID();
    
    
    
    
    
};
#endif
