#include "../src/parm.hpp"
#include "StrucFac.hpp"
#include <armadillo>
#include <vector>
#include <iterator>
#include <fstream>
#include <sstream>

int main(int argc, char* argv[]) {
    parm p_(argv);
    create_output_directory(p_);

    for (unsigned int b = 0; b < p_.Temperature.size(); b++) {
        p_.T = p_.Temperature[b];
        std::cout << "Temperature:\t" << p_.T << std::endl;
    	std::ostringstream in;
    	in << p_.out_dir << "/config_T" << p_.T << ".txt";
    	/*  open file  */
    	std::fstream ifs;
    	ifs.open(in.str().c_str());
    	if(ifs.fail()) {
    	    std::cout << "File not found at\t" << in.str().c_str() << std::endl << "Aborting!" << std::endl;
    	    exit(EXIT_FAILURE);
    	}
    	int s(0), n(0), save(p_.sweep/(2*p_.Nconf));
	std::string line;
    	arma::mat Sq(p_.Nconf,p_.L*p_.L);
    	//int flag = 0;
    	while (std::getline(ifs, line)) {
    	    if ((s >= (p_.sweep/2)) && ((s % save) == 0)) {
    	    	std::istringstream iss(line);
    	        std::vector<double> temp = std::vector<double>{std::istream_iterator<double>(iss), std::istream_iterator<double>()};
    	    	arma::mat S_conf(p_.L*p_.L,1);
    	        for (int i = 0; i < p_.L*p_.L; i++) {
    	            S_conf(i,0) = temp[i];
    	        }
    	    	Sq.row(n) = StrucFac(S_conf,p_.L,p_.L,1).t();
    	    	n++;
    	    }
    	    s++;
    	}
    	ifs.close();
    	std::cout << "Configs loaded.. No. of configs = " << s << std::endl;
    	std::cout << "No. of configs used for averaging = " << n << std::endl;
    	std::cout << "Saving data.." << std::endl;

    	arma::rowvec sfac = arma::sum(Sq,0) / n;
    	std::ostringstream sf;//, pd;
    	sf << p_.out_dir << "sfac_T" << p_.T << ".txt";
    	std::ofstream sff(sf.str().c_str());
    	sff << sfac.t();
    	sff.close();
    }

    return 0;
}

