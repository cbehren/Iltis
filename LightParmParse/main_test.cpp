#include "LParmParse.H"
#include <string>
#include <iostream>
#include <vector>

int main(int argc,char *args[])
{
LParmParse::Initialize(argc-1,args+2,"inputs");
std::string mystring;
int myint;
double mydouble;
std::vector<int> myint_array;
std::vector<double> mydouble_array;
LParmParse pp;
pp.query("test.string",mystring);
pp.query("test.int",myint);
pp.query("test.double",mydouble);
pp.queryarr("test.int_array",myint_array);
pp.queryarr("test.double_array",mydouble_array);

std::cout << "test.string "<< mystring << std::endl;
std::cout << "test.int "<< myint << std::endl;
std::cout << "test.double "<< mydouble << std::endl;
std::cout << "test.double_array: " << std::endl;
for(unsigned int i=0;i<mydouble_array.size();i++)
    std::cout << "\t"<<mydouble_array[i] << std::endl;
std::cout << "test.int_array: " << std::endl;
for(unsigned int i=0;i<myint_array.size();i++)
    std::cout << "\t"<< myint_array[i] << std::endl;


return 0;
}
