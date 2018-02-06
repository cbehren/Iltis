#include "BaseCell.H"

//implements BaseCell class. Nothing spectactular here.
BaseCell::BaseCell()
{
    density = 0.0;
    for(int i=0;i<3;i++)
        velocity[i] = 0.0;
    temperature = 0.0;
    dust_density = 0.0;
    turbulent_velocity = 0.0;
    
    
}

std::ostream& operator<<(std::ostream& stream,const BaseCell& cell)
{
    stream << cell.density << " " << cell.velocity[0] << " " << cell.velocity[1] << " " << cell.velocity[2] << " " << cell.temperature << " " << cell.dust_density << std::endl;
    return stream;
    
}
