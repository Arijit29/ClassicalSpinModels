#include <iterator>
#include "Ising2dMC.hpp"

class observables : protected Ising2dMC {
    protected:
        arma::cube S_conf;
        arma::vec E_conf;
        arma::vec M_conf;
        double E_av;
        double M_av;
        double Cv;
        double MagSusp;
    public:
        observables(parm&);
        void CalculateConfigAvgs(parm&);
        ~observables();
};
