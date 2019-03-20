#include <vector>
#include <string>
class Histogram
{
public:
    Histogram(double left, double right, double nbin);
    ~Histogram();
    inline void add(double value, double weight)
    {
        int index = (value-left_edge)/dx;
        if(index<0 || index>=nbins)
        {
            ignored+=value*weight;
        }
        else
        {
            values[index]+=weight;
        }
        
    };
    void add_histogram(const Histogram& hist);
    int write(std::string directory,std::string filename);
    void get_data(std::vector<double>& bin_centers,std::vector<double>& data);
    inline int get_nbins(){return nbins;};
    inline int get_left(){return left_edge;};
    inline int get_right(){return right_edge;};
    void reduce();
private:
    double *values = NULL;
    double left_edge=0.;
    double right_edge=0.;
    double ignored=0.;
    int nbins=0;
    double dx=0.;
    
};
