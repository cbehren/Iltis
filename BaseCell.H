#ifndef __BASE_CELL_
#define __BASE_CELL_
#include <iostream>

//defines base class for hydro cells. 

class BaseCell{
public:
  double density;
  double dust_density;
  double velocity[3];
  double temperature;
  double turbulent_velocity;
  friend std::ostream& operator<<(std::ostream& stream,const BaseCell& cell);//for pretty printing cells
  BaseCell();
    
};
#endif
