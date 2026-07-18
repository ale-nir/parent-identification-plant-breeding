#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_lod_matrix_cpp(IntegerMatrix bin_matrix) {

    int n = bin_matrix.nrow();
    int m = bin_matrix.ncol();

    NumericVector freqs(m);

    // fréquences alléliques
    for(int j=0;j<m;++j){
        double s = 0.0;
        for(int i=0;i<n;++i)
            s += bin_matrix(i,j);
        freqs[j] = s / n;
    }

    NumericMatrix L(n,n);

    for(int i=0;i<n;++i){

        L(i,i) = NA_REAL;

        for(int j=i+1;j<n;++j){

            double lod = 0.0;

            for(int k=0;k<m;++k){

                double f = freqs[k];

                int r1 = bin_matrix(i,k);
                int r2 = bin_matrix(j,k);

                if(r1==1 && r2==1)
                    lod += log((f+(1+f)*(1-f))/(f*pow(2-f,2)));

                else if(r1!=r2)
                    lod += log(1/(2-f));

                else
                    lod += log(1/(1-f));
            }

            L(i,j)=lod;
            L(j,i)=lod;
        }
    }

    return L;
}
