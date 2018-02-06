#ifndef __EMISSION_LIST_H
#define __EMISSION_LIST_H

#include <vector>
#include <iostream>

enum EmissionModelEnum{EM_DEFAULT};
enum SpectrumModelEnum{SP_ZZ10,SP_DELTA,SP_GAUSSIAN_FIXED};


class Emitter{
	public:
	Emitter(){for(int i=0;i<3;i++) position[i]=0.0;emissivity=0.0;index=-1;width=0.0;mass=0.0;hindex=0;findex=0;sfr=0.0;}
	double position[3];
	double mass;
	double emissivity;
        double width;
        int index;
	int findex;
	int hindex;
	double sfr;

	private:


};

class EmissionList{

	public:
	EmissionList(std::string fname,EmissionModelEnum m=EM_DEFAULT,SpectrumModelEnum s=SP_DELTA,double minLuminosity=0.0);
	EmissionList();
	void readFromFile(std::string fname,EmissionModelEnum m=EM_DEFAULT,SpectrumModelEnum s=SP_DELTA,double minLuminosity=0.0,bool applyBoundingBox=false);
	void writeToFile(std::string fname);
	Emitter& operator[](int index){return e[index];}
	std::vector<double> relativeEmissivity();
        std::vector<double> absoluteEmissivity();
	bool isInside(double x,double y,double z);
	void setUpBoundingBox(double lox,double loy,double loz,double hix,double hiy,double hiz);
	void setFixedWidth(double w);

        unsigned int size();


	private:
	void init();
	double boundingBoxLo[3];
	double boundingBoxHi[3];
	double fixedWidth;
	std::vector<Emitter> e;




};
#endif
