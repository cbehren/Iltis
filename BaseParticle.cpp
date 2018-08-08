#include "BaseParticle.H"
#include <cmath>
long BaseParticle::n=0;

const int BaseParticle::nfloats=15;
const int BaseParticle::nlongs=4;

BaseParticle::BaseParticle()
{
    for(int i=0;i<3;i++)
    {
     x[i] = 0.0;
     k[i] = 0.0;
     lsp[i] = 0.0;
    }
    order = 0.0;
    number_of_scatterings = 0;
    optical_depth = 0.0;
    optical_depth_seen = 0.0;
    weight = 1.0;
    
    path_length = 0.0;
    status = BP_UNINITIALIZED;
    frequency = 0;
    id = nextID();
    
    emitter = 0;
    bias = 1.0;

    
}
// #define DEBUG_MOVE
void BaseParticle::move(double pathlength)
{
#ifdef DEBUG_MOVE
    double xtemp[3];
    for(int i=0;i<3;i++)
        xtemp[i]=x[i];
    for(int i=0;i<3;i++)
        xtemp[i] += k[i]*pathlength;
    if(xtemp[0]==x[0] && xtemp[1]==x[1] && xtemp[2]==x[2])
        std::cout << "WARNING: step is empty" << std::endl;
#endif
    for(int i=0;i<3;i++)
        x[i] += k[i]*pathlength;
    this->path_length += pathlength;   
    
}
void BaseParticle::move_to(const double pos[3])
{
    double d=0;
    for(int i=0;i<3;i++)
        d += pow(x[i]-pos[i],2);
    d = sqrt(d);
    for(int i=0;i<3;i++)
        x[i] = pos[i];
    path_length += d;
    
}


std::ostream& operator<<(std::ostream& stream,const BaseParticle& p)
{
    stream << p.id << " : " << p.frequency << " " << p.number_of_scatterings << " "  << p.x[0] << " " << p.x[1] << " " << p.x[2] << " " << p.k[0] << " " << p.k[1] << " " << p.k[2] << " " << p.status << " " << p.order <<  std::endl;
    return stream;
    
}
long BaseParticle::nextID()
{
    return n++;
}

