#pragma once

#include <armadillo>
#include <omp.h>
#include "RandNum.hpp"
#include "parm.hpp"

class Ising2dMC {
    protected:
        arma::mat S;
        arma::mat H;
    public:
        Ising2dMC(parm&);
        void GenerateHamiltonian(parm&);
        void InitMC(parm&);
        void MonteCarlo(parm&);
        double CalculateLocalEnergyDiff(parm&, double&, int&);
        ~Ising2dMC();
};
