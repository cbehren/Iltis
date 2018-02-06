struct line_of_sight{
 double x[3];   
};
struct metadata{
  double lengthOfBox;//physical boxlength in cm
  double redshift;  
  double redshiftShifted;  
  double cutOffLength;
  int number_of_instruments;
  std::vector<line_of_sight> observationDirections;
};
