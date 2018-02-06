#include "EmissionList.H"
#include "stdio.h"
#include "math.h"

EmissionList::EmissionList(std::string fname,EmissionModelEnum m,SpectrumModelEnum s,double minLuminosity)
{
	init();
	readFromFile(fname,m,s,minLuminosity);

}
EmissionList::EmissionList()
{
	init();
}
void EmissionList::readFromFile(std::string fname,EmissionModelEnum m,SpectrumModelEnum s,double minLuminosity,bool applyBoundingBox)
{
	FILE *fp = fopen(fname.c_str(),"r");
        if(!fp)
            std::cout << "Could not open EmissionList file!" << fname.c_str() << "\n";
	char line[1000];
	double x,y,z,mass,sfr;
	int id,fid;

	fgets(line,1000,fp);
	int count=-1;
	sscanf(line,"#%d",&count);
	e.clear();
	e.reserve(count);
	int n=0;
	while(fgets(line,1000,fp))
	{
		if(sscanf(line,"%d %d %le %le %le %le %le\n",
		       &id,&fid,&x,&y,&z,&mass,&sfr)!=7)
                {
                    if(n==count)
                        break;
                    else
                    {
                        std::cout << "Error: Found more emitters in file than expected!\n" << std::endl;
                        std::abort();
                    }
                    
                }
		if(sfr<minLuminosity || (applyBoundingBox && !isInside(x,y,z)))
		{
		  count--;
		  n++;
		  continue;		  
		}
		Emitter entry;
		entry.position[0]  = x;
		entry.position[1]  = y;
		entry.position[2]  = z;
                entry.index = n;
		entry.findex=fid;
		entry.hindex=id;
		entry.mass=mass;
		entry.sfr=sfr;
		switch(m)
		{
			case EM_DEFAULT:
				entry.emissivity = sfr;//emissivitiy is in units of 1e42
		}
		switch(s)
		{
                    case SP_DELTA:
                        entry.width=0.0;                        
                        break;
                    case SP_ZZ10:
			entry.width = 31.9*pow((mass/1.e10),0.333333); //km/s, following ZZ2010
                        break;
                    case SP_GAUSSIAN_FIXED:
                        if(fixedWidth<0.0)
                        {
                            printf("Fixed Width is negative.\n");
                            std::abort();
                            
                        }
			entry.width = fixedWidth;
                        
                        break;
                    default:
                        printf("Unknown model for spectrum!\n");
                        std::abort();
		}
		e.push_back(entry);  
		n++;
	}
	fclose(fp);
	//write the list to a new file
	std::string outname(fname+"_used");
	writeToFile(outname);


}
void EmissionList::writeToFile(std::string fname)
{
  FILE *fp = fopen(fname.c_str(),"w");
  if(!fp)
    std::cout << "Could not open EmissionList file for writing!" << fname.c_str() << "\n";
  fprintf(fp,"# %lu\n",e.size());
  for(unsigned int i=0;i<e.size();i++)
  {
    Emitter& mye = e[i];
    fprintf(fp,"%d %d %le %le %le %le %le\n",mye.hindex,mye.findex,mye.position[0],mye.position[1],mye.position[2],mye.mass,mye.sfr);
    
  }
	fclose(fp);
     
  
}


std::vector<double> EmissionList::relativeEmissivity()
{
    double sum = 0.0;
    std::vector<double> rel;
    rel.resize(e.size());
    for(unsigned int i=0;i<e.size();i++)
    {
        sum+=e[i].emissivity;
        
    }
    for(unsigned int i=0;i<e.size();i++)
    {
        rel[i]=e[i].emissivity/sum;
        
    }
    return rel;
    
}

std::vector<double> EmissionList::absoluteEmissivity()
{
    std::vector<double> rel;
    rel.resize(e.size());
    for(unsigned int i=0;i<e.size();i++)
    {
        rel[i]=e[i].emissivity;
        
    }
    return rel;
    
}


unsigned int EmissionList::size()
{
    return e.size();
}
void EmissionList::init()
{
  for(int i=0;i<3;i++)
    boundingBoxLo[i]=0.0;
  for(int i=0;i<3;i++)
    boundingBoxHi[i]=1.0;
  fixedWidth = -1.0;
  
}
void EmissionList::setUpBoundingBox(double lox,double loy,double loz,double hix,double hiy,double hiz)
{
  boundingBoxLo[0]=lox;
  boundingBoxLo[1]=loy;
  boundingBoxLo[2]=loz;
  
  boundingBoxHi[0]=hix;
  boundingBoxHi[1]=hiy;
  boundingBoxHi[2]=hiz;
  
}
void EmissionList::setFixedWidth(double w)
{
  fixedWidth = w;
  
}

bool EmissionList::isInside(double x,double y,double z)
{
  if(x<boundingBoxLo[0] || x>=boundingBoxHi[0])
    return false;
  if(y<boundingBoxLo[1] || y>=boundingBoxHi[1])
    return false;
  if(z<boundingBoxLo[2] || z>=boundingBoxHi[2])
    return false;

  return true;
  
}


