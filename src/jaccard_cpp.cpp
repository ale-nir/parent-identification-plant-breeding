#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compute_jaccard_matrix_cpp(const NumericMatrix& bin) {

    int n = bin.nrow();
    int p = bin.ncol();

    NumericMatrix J(n, n);

    for (int i = 0; i < n; ++i) {

        J(i, i) = 1.0;

        for (int j = i + 1; j < n; ++j) {

            int m11 = 0;
            int uni = 0;

            for (int k = 0; k < p; ++k) {

                int a = bin(i, k);
                int b = bin(j, k);

                if (a | b)
                    ++uni;

                if (a & b)
                    ++m11;
            }

            double score =
                (uni == 0) ? 0.0 :
                static_cast<double>(m11) / uni;

            J(i, j) = score;
            J(j, i) = score;
        }
    }

    return J;
}