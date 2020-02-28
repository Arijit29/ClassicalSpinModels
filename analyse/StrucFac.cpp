#include "StrucFac.hpp"

arma::vec StrucFac(arma::mat& in, const unsigned int& Lx, const unsigned int& Ly, const unsigned int& Lz) {
    const unsigned int L =  Lx * Ly * Lz;
    int comp = in.n_cols;
    std::cout << comp << std::endl;
    arma::cx_mat out(L,comp,arma::fill::zeros);
    arma::vec sfac(L);
    for(unsigned int n = 0; n < L; n++) {
	int z = n / (Lx*Ly);
	int perp = n % (Lx*Ly);
	int x = perp % Lx;
	int y = perp / Lx;
	double kx = ((2.0 * x / Lx)) * M_PI;
	double ky = ((2.0 * y / Ly)) * M_PI;
	double kz = ((2.0 * z / Lz)) * M_PI;
	for(unsigned int s = 0; s < comp; s++) {
	    for(unsigned int i = 0; i < L; i++) {
	    	int iz = i / (Lx*Ly);
	    	int iperp = i % (Lx*Ly);
	    	int ix = iperp % Lx;
	    	int iy = iperp / Lx;
		double phase = (kx * ix) + (ky * iy) + (kz * iz);
	    	out(n,s) += std::polar(in(i,s), phase);
	    }
	    out (n,s) /= ((double)L);
	}
	sfac(n) = std::real(arma::cdot(out.row(n),out.row(n)));
    }
    return sfac;	
}

