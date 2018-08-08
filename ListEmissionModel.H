#include "BaseEmissionModel.H"
#include "EmissionList.H"
#include <vector>
//defines a emission model for emitting from a list of emitters provided as a text file.

class ListEmissionModel : public BaseEmissionModel{
public:
    void setup();
    void launch_bunch(BaseParticleVector& particles, const BaseEmissionLine& line, const BaseDataset& ds, double dt);
protected:
    EmissionList emitters;
    std::vector<double> relative_emissivity;
    std::vector<double> absolute_emissivity;
       
};
