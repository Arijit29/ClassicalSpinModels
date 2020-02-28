#include "observables.hpp"

observables::observables(parm& p_) : Ising2dMC(p_), S_conf(p_.L, p_.L, p_.Nconf, arma::fill::zeros), E_conf(p_.Nconf, arma::fill::zeros), M_conf(p_.Nconf, arma::fill::zeros), E_av(0.0), M_av(0.0), Cv(0.0), MagSusp(0.0) {}

void observables::CalculateConfigAvgs(parm& p_) {

    std::ofstream logfile(p_.log.c_str(), std::ios::out | std::ios::app);
    std::ofstream obs (p_.out_dir+"observables_h0.dat");

    //double junk;
    for (unsigned int b = 0; b < p_.Temperature.size(); b++) {
        p_.T = p_.Temperature[b];
        logfile << "Temperature:\t" << p_.T << std::endl;
        std::ostringstream fn;
        fn << p_.out_dir << "config_T" << p_.T << ".txt";
        std::ifstream file (fn.str().c_str());
        Cv = 0.0;
        MagSusp = 0.0;
        //for (int n = 0; n < p_.Nconf; n++) {
        int s(0), n(0), save(p_.sweep/(2*p_.Nconf));
        std::string line;
        logfile << "Reading config file.." << std::endl;
        while (std::getline(file, line)) {
            if ((s >= (p_.sweep/2)) && ((s % save) == 0)) {
                std::istringstream iss(line);
                std::vector<double> temp = std::vector<double>{std::istream_iterator<double>(iss), std::istream_iterator<double>()};
                for (int i = 0; i < p_.L; i++) {
                    S_conf(i/p_.L, i%p_.L, n) = temp[i];
                }

                S = S_conf.slice(n);
                std::cout << S.t() << std::endl;
                GenerateHamiltonian(p_);
                E_conf(n) = arma::accu(H);
                M_conf(n) = arma::accu(S_conf.slice(n));
                n++;
            }
            s++;
        }
        file.close();
        E_av = arma::sum(E_conf) / ((double)p_.Nconf);
        M_av = arma::sum(M_conf) / ((double)p_.Nconf);
        for (int n = 0; n < p_.Nconf; n++) {
            Cv += (pow(E_conf(n), 2) / ((double)p_.Nconf));
            MagSusp += (pow(M_conf(n), 2) / ((double)p_.Nconf));
        }
        Cv = (Cv - (E_av*E_av)) / ((double)p_.T*p_.T);
        MagSusp = (MagSusp - (M_av*M_av)) / ((double)p_.T);

        E_av /= ((double)p_.L*p_.L);
        M_av /= ((double)p_.L*p_.L);
        Cv /= ((double)p_.L*p_.L);
        MagSusp /= ((double)p_.L*p_.L);
        obs << p_.T << '\t' << std::abs(M_av) << '\t' << MagSusp << '\t' << E_av << '\t' << Cv << std::endl;
    }
    obs.close();
    logfile.close();
    std::cout << std::endl;
}

observables::~observables() {}
