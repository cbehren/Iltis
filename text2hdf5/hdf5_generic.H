#ifdef OLD_HEADER_FILENAME
#include <iostream.h>
#else
#include <iostream>
#endif
#include <string>
#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
#endif // H5_NO_STD
#endif
#include "H5Cpp.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

//Implements a very simple interface to read hdf5 files. Limited to one dimensional tables.
//You need to know the type of the table you want to open in advance.
//
//
//

class hdf5_file{
private: 

  void init();
  H5File *fp;
public:
  ~hdf5_file();
    hdf5_file();
  hdf5_file(char* fn,char mode);//open a hdf5 file. modes 'r','w' and 'o' (create new file and overwrite)
  void open_dataset(char* name,long *n,double **data);//read dataset <name> into <data>. n holds the number of datapoints.
  void open_dataset(char* name,long *n,float **data);//same but for float tables.
  void open_dataset(char* name,long *n,long **data);//same but for long tables.
  void open_dataset(char* name,long *n,unsigned long **data);//"" for unsigned long

  void open_attribute(char* name, long* data);//open attribute <name> and return its value in <data>
  void open_attribute(char* name, double* data);//same but for long
  void open_attribute(char* name, char* data, int* n);//same but for strings.

  void create_attribute(char *name, char* data);
  void create_attribute(char *name, double data);
  
  void create_dataset(char* name,long n,double *d);//write dataset <name> from <data>. n holds the number of datapoints.
  void create_dataset(char* name,long n,long *d);//same but for long tables.
  void create_dataset(char* name,long n,int *d);//same but for long tables.
  
  H5File* getFile();
  
  
  
};
