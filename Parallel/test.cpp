#include "Parallel.H"
#include <iostream>
int main(int argc, char* args[])
{
    Parallel::StartParallel();
    std::vector<double> vec(Parallel::NProcs());
    vec[Parallel::MyProc()] = Parallel::MyProc();
    Parallel::ReduceDoubleSum(vec.data(),Parallel::NProcs());
    std::cout << "This is process " << Parallel::MyProc() << std::endl;
    if(Parallel::IOProcessor())
    {
        for(int i=0;i<Parallel::NProcs();i++)
            std::cout << vec[i] << " ";
        std::cout << std::endl;

    }
    Parallel::EndParallel();

}
