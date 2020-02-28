#include <iostream>
#include <sstream>
#include <fstream>
#include <armadillo>
#include "../src/RandNum.hpp"

const int N = 100000;
const int Nbin = 10;
const int L = 2;

void FreqTest() { /** Check the distribution of the random numbers **/

    long seed = (long)time(NULL);

    arma::vec sample(N, arma::fill::zeros);
    for (int n = 0; n < N; n++) {
        sample(n) = ran2(&seed);
    }

    arma::vec hist(Nbin, arma::fill::zeros);
    for (int n = 1; n <= Nbin; n++) {
        double x1 = (n - 1) / ((double)Nbin);
        double x2 = n / ((double)Nbin);
        for (int i = 0; i < N; i++) {
            double x = sample(i);
            if((x >= x1) && (x <= x2))
                hist(n-1)++;
        }
    }

    double norm = sum(hist); // Normalise
    hist /= norm;

    std::ostringstream fn;
    fn << "histogram_N" << N << "_Nbin" << Nbin << ".txt";
    std::ofstream file (fn.str().c_str());
    for (int n = 0; n < Nbin; n++) {
        double x = (n + 1) / ((double)Nbin);
        file << x << '\t' << hist(n) << std::endl;
    }
    file.close();
}

void intSampling() { /* Test for extraction of integers */

    long seed = (long)time(NULL);
    arma::Col<int> int_sample(N, arma::fill::zeros);
    for (int n = 0; n < N; n++) {
        int_sample(n) = floor(ran2(&seed)*L);
    }

    arma::Col<double> int_hist(L, arma::fill::zeros);
    for (int l = 0; l < L; l++) {
        for (int i = 0; i < N; i++) {
            if(int_sample(i) == l)
                int_hist(l)++;
        }
    }

    int int_norm = sum(int_hist); // Normalise
    int_hist /= (double)int_norm;

    std::ostringstream fn1, fn2;
    fn1 << "int_sample_N" << N << "_L" << L << ".txt";
    std::ofstream file1 (fn1.str().c_str());
    for (int n = 0; n < N; n++) {
        file1 << n << '\t' << int_sample(n) << std::endl;
    }
    file1.close();
    fn2 << "intHist_N" << N << "_L" << L << ".txt";
    std::ofstream file2 (fn2.str().c_str());
    for (int n = 0; n < L; n++) {
        file2 << n << '\t' << int_hist(n) << std::endl;
    }
    file2.close();
}

int main() {
    /* Test ran2 */
    long init = -(long)time(NULL);
    ran2(&init);

    FreqTest();
    intSampling();
    return 0;

}
