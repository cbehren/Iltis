#ifndef __ANISO_SHELL_H
#define __ANISOL_SHELL_H
#include "SphericalShellData.H"

//implements a spherical shell geometry with two holes in direction -z, +z in the shell, parametrized by an opening angle.


class AnisoShellData : public SphericalShellData{
public:
    virtual const BaseCell* data_at(const double x[3], const double k[3],double &pathlength,DatasetStatus &status) const;
    virtual int setup();
private:
    double opening_angle;//opening angle of the hole in the shell, in deg.
};
    










#endif
