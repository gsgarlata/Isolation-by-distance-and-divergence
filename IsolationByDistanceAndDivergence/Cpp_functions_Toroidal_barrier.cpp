#include <Rcpp.h>
using namespace Rcpp;
//using namespace std;

inline double deltaFunc(int x){
    const double deltaD = x==0 ? 1 : 0.5;
    return deltaD;
}


// [[Rcpp::export]]
double sumSij(int i, int j, int d1, int d2, double migr) {

        const int x = d1 / 2;
        const int y = d2 / 2;
        const double pi = 2 * acos(0.0);

            double total = 0;
        for(int k = 0; k <= x; ++k) {

            const double pi_k = pi * k;
            double deltaD1 = deltaFunc(k);

            for(int l = 0; l <= y; ++l) {

                const double pi_l = pi * l;

                double fkl;
                if (k == 0 && l == 0) {
                      fkl = 0;
                    }else{
                      fkl = pow((1 - (migr * (1 - cos(2.0 * pi_k / d1 )))),2.0) * pow((1 - (migr * (1 - cos(2.0 * pi_l / d2)))),2.0);
                    }

                double a1 = 1 -  ( cos(2.0 * pi_k * i / d1 ) * cos(2.0 * pi_l * j / d2) );

                double deltaD2 = deltaFunc(l);

                double res = (fkl * a1) / (deltaD1 * deltaD2 * (1 - fkl));

                total += res;
          }
          }
                return total;
                      }


// [[Rcpp::export]]
NumericVector series_Sij(NumericVector xvec, NumericVector yvec, int d1, int d2, double migr){

    int size = xvec.length();
    NumericVector series_res(size);

for(double i=0;i<size;i++){
    series_res[i] = sumSij(xvec[i], yvec[i], d1, d2, migr);
}

return series_res;
}

// [[Rcpp::export]]
NumericVector series_walk(NumericVector anc_pos, int t0_pos, NumericMatrix time_matrix){

int size = anc_pos.length();
NumericVector pAt(size);

for(double i=0;i<size;i++){
  pAt[i]  = time_matrix(anc_pos[i], t0_pos);
}

return pAt;
}
