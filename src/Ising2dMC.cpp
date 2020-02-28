#include "Ising2dMC.hpp"

Ising2dMC::Ising2dMC(parm& p_) : S(p_.L, p_.L, arma::fill::zeros), H(p_.L, p_.L, arma::fill::zeros) {
}

void Ising2dMC::GenerateHamiltonian(parm& p_) {
    /** Bulk Hamiltonian **/
    for (int i = 1; i < p_.L-1; i++) {
        for (int j = 1; j < p_.L-1; j++) {
            H(i, j) = (-p_.J * S(i, j) * (S(i, j-1) + S(i, j+1) + S(i-1, j) + S(i+1, j))) + (-p_.h * S(i, j));
        }
    }

    /** Periodic Boundary Conditions **/
    for (int j = 1; j < p_.L-1; j++) {
        H(0, j) = (-p_.J * S(0, j) * (S(0, j-1) + S(0, j+1) + S(p_.L-1, j) + S(1, j))) + (-p_.h * S(0, j));
        H(j, 0) = (-p_.J * S(j, 0) * (S(j-1, 0) + S(j+1, 0) + S(j, p_.L-1) + S(j, 1))) + (-p_.h * S(j, 0));
        H(j, p_.L-1) = (-p_.J * S(j, p_.L-1) * (S(j-1, p_.L-1) + S(j+1, p_.L-1) + S(j, p_.L-2) + S(j, 0))) + (-p_.h * S(j, p_.L-1));
        H(p_.L-1, j) = (-p_.J * S(p_.L-1, j) * (S(p_.L-1, j-1) + S(p_.L-1, j+1) + S(p_.L-2, j) + S(0, j))) + (-p_.h * S(p_.L-1, j));
    }
    H(0, 0) = (-p_.J * S(0, 0) * (S(0, p_.L-1) + S(0, 1) + S(p_.L-1, 0) + S(1, 0))) + (-p_.h * S(0, 0));
    H(p_.L-1, 0) = (-p_.J * S(p_.L-1, 0) * (S(p_.L-1, p_.L-1) + S(p_.L-1, 1) + S(p_.L-2, 0) + S(0, 0))) + (-p_.h * S(p_.L-1, 0));
    H(p_.L-1, p_.L-1) = (-p_.J * S(p_.L-1, p_.L-1) * (S(p_.L-1, p_.L-2) + S(p_.L-1, 0) + S(p_.L-2, p_.L-1) + S(0, p_.L-1))) + (-p_.h * S(p_.L-1, p_.L-1));
    H(0, p_.L-1) = (-p_.J * S(0, p_.L-1) * (S(0, p_.L-2) + S(0, 0) + S(p_.L-1, p_.L-1) + S(1, p_.L-1))) + (-p_.h * S(0, p_.L-1));
    /*********************************/
}

void Ising2dMC::InitMC(parm& p_){
    // Random initial state
    for (int i = 0; i < p_.L*p_.L; i++) {
        int i_x = i / p_.L;
        int i_y = i % p_.L;
        S(i_x, i_y) = std::pow(-1, rand()%2);
    }

    //Start Monte-Carlo
    MonteCarlo(p_);
}

double Ising2dMC::CalculateLocalEnergyDiff(parm& p_, double& S_local_prev, int& i){
    double res(0.0);
    arma::Col<int> xcords(2);
    arma::Col<int> ycords(2);

    int x = i / p_.L;
    int y = i % p_.L;
    /** Periodic Boundary **/
    if ((x > 0) && (x < p_.L-1))
        xcords << x-1 << x+1;
    else {
        if (x == 0)
            xcords << p_.L-1 << 1;
        if (x == p_.L-1)
            xcords << p_.L-2 << 0;
    }
    if ((y > 0) && (y < p_.L-1))
        ycords << y-1 << y+1;
    else{
        if (y == 0)
            ycords << p_.L-1 << 1;
        if (y == p_.L-1)
            ycords << p_.L-2 << 0;
    }
    res = ((-p_.J * (S(x, ycords(0)) + S(x, ycords(1)) + S(xcords(0), y) + S(xcords(1), y))) - p_.h) * (S(x, y) - S_local_prev);
    return res;
}

void Ising2dMC::MonteCarlo(parm& p_){

    std::ofstream logfile(p_.log.c_str(), std::ios::out | std::ios::app);

    logfile << "Monte Carlo annealing.." << std::endl;
    srand(time(NULL));
    long init = -(long)time(NULL);
    ran2(&init); // Initialise the random number generator.
    GenerateHamiltonian(p_);

    for (unsigned int b = 0; b < p_.Temperature.size(); b++) {
        p_.T = p_.Temperature[b];
        logfile << p_.T << " :\t";
        std::ostringstream fn;
        fn << p_.out_dir << "config_T" << p_.T << ".txt";
        std::ofstream config(fn.str().c_str());
        long seed = (long)rand()/(RAND_MAX);
        //int cnt = 0;
        /** Equilibriate **/
        for (int s = 0; s < p_.sweep; s++) { // MC sweep
            for (int i = 0; i < p_.L*p_.L; i++) {// System sweep
                // Select a random site
                int j = floor(ran2(&seed) * (p_.L*p_.L));
                int x = j / p_.L;
                int y = j % p_.L;
                // Remember the current local configuration and free energy
                double S_local_prev = S(x, y);
                // Update site (local update)
                S(x, y) = std::pow(-1, floor(ran2(&seed)*2));
                // Update the local Hamiltonian
                double delta_F = CalculateLocalEnergyDiff(p_, S_local_prev, j);
                //Metropolis decision tree
                if(delta_F > 0) {
                    double toss = ran2(&seed);
                    double beta = 1.0 / p_.T;
                    if(exp(-beta * delta_F) < toss) {
                        S(x, y) = S_local_prev;
                    }
                }
            }
            /* Save configuration */
            //int save = (p_.sweep/2) + (cnt * (p_.sweep/(2*p_.Nconf)));
            //if((s >= save) && (cnt < p_.Nconf)) {
                //cnt ++;
                for (int i = 0; i < p_.L*p_.L; i++) {
                    int x = i / p_.L;
                    int y = i % p_.L;
                    //config << x << '\t' << y << '\t' << S(x, y) << '\t' << F << std::endl;
                    //config << x << '\t' << y << '\t' << S(x, y) << std::endl;
                    config << S(x, y) << '\t';
                    //if(y == (p_.L-1))
                        //config << std::endl;
                }
                config << std::endl;
            //}
        }
        config.close();
        logfile << "Sweeps over.. Configs saved." << std::endl;
    }
    logfile.close();
}

Ising2dMC::~Ising2dMC(){};
