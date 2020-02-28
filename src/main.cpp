#include <iostream>
#include "observables.hpp"

int main(int argc, char* argv[]){
    parm p_(argv);
    create_output_directory(p_);
    Ising2dMC MC_(p_);
    std::cout << "Starting Monte Carlo using Metropolis algorithm.." << std::endl;
    std::cout << "Monte Carlo log to be found in -> " << p_.log << std::endl;
    MC_.InitMC(p_);
    std::cout << "Monte Carlo runs over. Calculating observables.." << std::endl;
    observables ob_(p_);
    ob_.CalculateConfigAvgs(p_);
    std::cout << "Code executed successfully." << std::endl;
    return 0;
}
