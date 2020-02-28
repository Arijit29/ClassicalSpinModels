#include "parm.hpp"

parm::parm (char *argv[]) : L(atoi(argv[1])), J(1.0), h(0.0), sweep(atoi(argv[2])), Nconf(atoi(argv[3])) {
    Temperature = {5.0, 4.5, 4.0, 3.75, 3.5, 3.25, 3.0, 2.8, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.8, 1.6, 1.4, 1.0 , 0.01};
}

void create_output_directory(parm& p_){
    /** Set paths for output directory and logfile **/
    std::ostringstream fn, fn1;
    fn << "./L" << p_.L << "/";

    fn1 << fn.str() << "log_h" << p_.h << ".dat";
    p_.log = fn1.str();

    /** Create logfile **/
    std::ofstream logfile(p_.log.c_str(), std::ios::out | std::ios::app);

    /** Create output directory **/
    int err = mkdir(fn.str().c_str(),0777);
    if (err == 0)
	    logfile << "Output directory -> " << fn.str() << " created" << std::endl;
    else {
        if (errno == EEXIST) {
            std::cout << "Output directory ->" << fn.str() << " already exists." << std::endl;
        } else {
            std::cout << "Cannot create output directory " << fn.str() << std::endl;
            exit(EXIT_FAILURE);
        }
    }

	logfile.close();
    p_.out_dir = fn.str();
}

parm::~parm() {}
