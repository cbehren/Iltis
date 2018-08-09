
struct metadata{
  double lengthOfBox;//physical boxlength in cm
  double redshift;  
  double redshift_shifted;  
  double cutOffLength;
  int number_of_instruments;
  std::vector<line_of_sight> observationDirections;
};
