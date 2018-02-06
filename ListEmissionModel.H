#include "BaseEmissionModel.H"
//defines a emission model for emitting from a list of emitters provided as a text file.

class ListEmissionModel : public BaseEmissionModel{
public:
    void setup(BaseParticleVector& particles,const BaseEmissionLine& line, const BaseDataset& ds);
    
private:

       
};
