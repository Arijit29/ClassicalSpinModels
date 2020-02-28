#pragma once

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cerrno>       // To handle error numbers
#include <iomanip>
#include <sstream>	    // To dynamically name files and folders
#include <sys/stat.h>	// To create folders from inside the program
#include <unistd.h>	    // To access directories in the filesystem

struct parm {
    public:
    const int L;
    const double J;
    const double h;
    double T;
    std::vector<double> Temperature;
    const int sweep;
    const int Nconf;
    std::string out_dir;
    std::string log;
    parm(char *argv[]);
    ~parm();
};

void create_output_directory(parm&);

