#include "hdf5_generic.H"
#include "string.h"
#include <cstdlib>
hdf5_file::hdf5_file()
{
 init(); 
}
hdf5_file::hdf5_file(char* fn,char mode)
{
  //TODO implement mode
  init(); 
 switch(mode)
 {
   case 'r':
    fp = new H5File(fn,H5F_ACC_RDONLY);
    break;
   case 'w':
    fp = new H5File(fn,H5F_ACC_RDWR);
    break;
   case 'o':
    fp = new H5File(fn,H5F_ACC_TRUNC);
    break; 
   default:
     std::cout << "Unknown file mode." << std::endl;
     std::abort();
 }
}
void hdf5_file::init()
{
  fp=NULL;
}

hdf5_file::~hdf5_file()
{
  if(fp)
 delete fp; 
}

void hdf5_file::open_dataset(char* name,long *n,double **data)
{
   std::string sname(name);
 DataSet *d = new DataSet( fp->openDataSet(sname));
 H5T_class_t type_class = d->getTypeClass();
 if( type_class != H5T_FLOAT )
  {
    std::cout << "float dataset has wrong type." << std::endl;
    std::abort();
  }
  FloatType intype = d->getFloatType();
  if(sizeof(double)!=intype.getSize())
  {
    std::cout << "Mismatch of float data type." << std::endl;
    std::abort();
  }
  size_t s = d->getInMemDataSize();
  double *buffer = new double[s/sizeof(double)];
  d->read(buffer,PredType::NATIVE_DOUBLE);//ensure endianness is set to the right value.
 delete d;
 *n = s/sizeof(double);
 *data = buffer;
  
}

void hdf5_file::open_dataset(char* name,long *n,float **data)//same but for float tables.
{
 std::string sname(name);
 DataSet *d = new DataSet( fp->openDataSet(sname));
 H5T_class_t type_class = d->getTypeClass();
 if( type_class != H5T_FLOAT )
  {
    std::cout << "float dataset has wrong type." << std::endl;
    std::abort();
  }
  FloatType intype = d->getFloatType();
  if(sizeof(float)!=intype.getSize())
  {
    std::cout << "Mismatch of float data type." << std::endl;
    std::abort();
  }
  size_t s = d->getInMemDataSize();
  float *buffer = new float[s/sizeof(float)];
  d->read(buffer,PredType::NATIVE_FLOAT);//ensure endianness is set to the right value.
 delete d;
 *n = s/sizeof(float);
 *data = buffer;
  
}

void hdf5_file::open_dataset(char* name,long *n,long **data)
{
 std::string sname(name);
 DataSet *d = new DataSet( fp->openDataSet(sname));
 H5T_class_t type_class = d->getTypeClass();
 if( type_class != H5T_INTEGER )
  {
    std::cout << "long dataset has wrong type." << std::endl;
    std::abort();
  }
  IntType intype = d->getIntType();
  if(sizeof(long)!=intype.getSize())
  {
    std::cout << "Mismatch of long data type." << std::endl;
    std::abort();
  }
  size_t s = d->getInMemDataSize();
  long *buffer = new long[s/sizeof(long)];
  d->read(buffer,PredType::NATIVE_LONG);
  delete d;
  *n = s/sizeof(long);
  *data = buffer;
 
  
}

void hdf5_file::open_dataset(char* name,long *n,unsigned long **data)
{
 std::string sname(name);
 DataSet *d = new DataSet( fp->openDataSet(sname));
 H5T_class_t type_class = d->getTypeClass();
 if( type_class != H5T_INTEGER )
  {
    std::cout << "long dataset has wrong type." << std::endl;
    std::abort();
  }
  IntType intype = d->getIntType();
  if(sizeof(unsigned long)!=intype.getSize())
  {
    std::cout << "Mismatch of long data type." << std::endl;
    std::abort();
  }
  size_t s = d->getInMemDataSize();
  unsigned long *buffer = new unsigned long[s/sizeof(unsigned long)];
  d->read(buffer,PredType::NATIVE_ULONG);
  delete d;
  *n = s/sizeof(unsigned long);
  *data = buffer;
 
  
}


//TODO somethings wrong with this part. I will fix it later.
// void hdf5_file::open_attribute(char* name, long* data)
// {
//    char *d = strrchr(name,'/');
//  int offset = d-name;
//  char dname[1000]; 
//  char *aname;
//  strncpy(dname,name,offset);
//  aname = d+1;
//  std::cout << name << " " << dname << " " << aname << std::endl;
//  std::string dname2(dname);
//  std::string aname2(aname);
//  
//  DataSet *dataset = new DataSet( fp->openDataSet(dname2));
//  Attribute myatt_out = dataset->openAttribute(aname2);
//  H5T_class_t type_class = myatt_out.getTypeClass();
//  if( type_class != H5T_INTEGER)
//   {
//     std::cout << "int attribute has wrong type." << std::endl;
//     std::abort();
//   }
//   IntType intype = myatt_out.getIntType();
//   if(sizeof(long)!=intype.getSize())
//   {
//     std::cout << "Mismatch of int data type." << std::endl;
//     std::abort();
//   }
//  long buffer=0;
//  myatt_out.read(intype,&buffer);
//  delete dataset;
//  *data =  buffer;
//   
// }
// void hdf5_file::open_attribute(char* name, double* data)
// {
//  char *d = strrchr(name,'/');
//  int offset = d-name;
//  char dname[1000]; 
//  char *aname;
//  strncpy(dname,name,offset);
//  aname = d+1;
//  std::cout << name << " " << dname << " " << aname << std::endl;
//  std::string dname2(dname);
//  std::string aname2(aname);
//  
//  DataSet *dataset = new DataSet( fp->openDataSet(dname2));
//  Attribute myatt_out = dataset->openAttribute(aname2);
//  H5T_class_t type_class = myatt_out.getTypeClass();
//  if( type_class != H5T_FLOAT )
//   {
//     std::cout << "float attribute has wrong type." << std::endl;
//     std::abort();
//   }
//   FloatType intype = myatt_out.getFloatType();
//   if(sizeof(double)!=intype.getSize())
//   {
//     std::cout << "Mismatch of float data type." << std::endl;
//     std::abort();
//   }
//  double buffer=0;
//  myatt_out.read(intype,&buffer);
//  delete dataset;
//  *data =  buffer;
//   
// }
void hdf5_file::open_attribute(char* name, char* data, int* n)
{
  std::cout << "Not yet implemented.";
  std::abort();
  
}


void hdf5_file::create_dataset(char* name,long n,double *d)//write dataset <name> from <data>. n holds the number of datapoints.
{
  std::string dname(name);
  hsize_t size = n;
  DataSpace space(1,&size,&size);
  DataSet data = fp->createDataSet(dname.c_str(), PredType::NATIVE_DOUBLE,space);
  data.write((void*)d, PredType::NATIVE_DOUBLE);
  data.close();
  
}
void hdf5_file::create_dataset(char* name,long n,long *d)//same but for long tables.
{
  std::string dname(name);
  hsize_t size = n;
  DataSpace space(1,&size,&size);
  DataSet data = fp->createDataSet(dname.c_str(), PredType::NATIVE_LONG,space);
  data.write((void*)d, PredType::NATIVE_LONG);
  data.close();
 
}

void hdf5_file::create_dataset(char* name,long n,int *d)//same but for int tables.
{
  std::string dname(name);
  hsize_t size = n;
  DataSpace space(1,&size,&size);
  DataSet data = fp->createDataSet(dname.c_str(), PredType::NATIVE_INT,space);
  data.write((void*)d, PredType::NATIVE_INT);
  data.close();
 
}

void hdf5_file::create_attribute(char *name, char* data)
{
    const H5std_string attr_name(name);
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    StrType strdatatype(PredType::C_S1, 256); // of length 256 characters
    const H5std_string strwritebuf (data);
    Attribute myatt_in = fp->createAttribute(attr_name, strdatatype, attr_dataspace);
    myatt_in.write(strdatatype, strwritebuf);
}

void hdf5_file::create_attribute(char *name, double data)
{
    const H5std_string attr_name(name);
    DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    Attribute myatt_in = fp->createAttribute(attr_name, PredType::NATIVE_DOUBLE, attr_dataspace);
    myatt_in.write(PredType::NATIVE_DOUBLE,&data);
}

H5File* hdf5_file::getFile()
{
  return fp;
  
  
}
