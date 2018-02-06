#include "main.H"
#include <iostream> 
#include <ctime>
int main(int argc,char* args[])
{
    Parallel::StartParallel();
    if(Parallel::IOProcessor())
        std::cout << "Light-weight Line Transfer Code - LLTC\n";
    //Initialize the parameter file reader
    if(argc<2)
    {
        std::cout << "You need to provide an inputs file."<< std::endl;
        std::abort();

    }
    std::string inputsfile(args[1]);

    LParmParse::Initialize(argc-2,args+2,inputsfile.c_str());
    //initialize the random number generator
    RNG::initialize();
    //create the simulation object
    BaseSimulation bs;
    //read the parameters
    bs.read_params();
    //setup the ingredients of the simulation
    bs.setup();
    //do steps until we are done
    double start = Parallel::second();
    while(bs.is_done()!=true)
        bs.do_step();
    //clean up and write output
    bs.cleanup();
    double end = Parallel::second();
    
    int number_of_photons;
    LParmParse pp;
    pp.get("number_of_photons",number_of_photons);
    if(Parallel::IOProcessor())
        std::cout << "Time needed for simulation " << end-start << " seconds, " << number_of_photons/(double)(end-start) << " photons per second.\n"; 
    
    Parallel::EndParallel();
    return 0;
}
