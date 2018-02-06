#include <iostream>
#include <cmath>

double TraversalLength(double SubCell[], const double LoS[],int* mustbecorrected)
{
  double ret = 0.0;
  double l[3] = {0,0,0};
 
  for(int i = 0; i < 3; i++)
  {
    if(LoS[i] == 0.0) l[i] = 100000.0;
    else{
      if(LoS[i] > 0.0){
        l[i] = (1.0 - SubCell[i]) / LoS[i];
      }
      else{
        l[i] = (0.0 - SubCell[i]) / LoS[i];
      }
    }
  }
  int minl = 0;
  if((l[0] >= l[1]) && (l[2] >= l[1])) minl = 1;
  if((l[0] >= l[2]) && (l[1] >= l[2])) minl = 2;
 
  for(int i = 0; i < 3; i++){
    ret += (l[minl] * LoS[i]) * (l[minl] * LoS[i]);
if(ret!=ret)
{
  std::cout << "TraversalLength(): "<< LoS[0] << " " << LoS[1] << " " << LoS[2] << " " << l[minl] << " "; 
  return -1;
}
  
  SubCell[i] = SubCell[i] + l[minl] * LoS[i];
    if(i == minl){
      if(LoS[i] < 0.0) SubCell[i] = 1.0;
      else SubCell[i] = 0.0;
    }
  }
 
    if(LoS[minl] < 0.0)
      *mustbecorrected = minl+1;
    else *mustbecorrected = 0;
  return sqrt(ret);
}