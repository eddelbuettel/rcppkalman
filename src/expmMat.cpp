// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppArmadillo.h>

/* Interface to expm package. */
typedef enum { Ward_2, Ward_1, Ward_buggy_octave } precond_type;
/* Matrix exponential exp(x), where x is an (n x n) matrix. Result z
 * is an (n x n) matrix. Mostly lifted from the core of fonction
 * expm() of package Matrix, which is itself based on the function of
 * the same name in Octave. */
void (*expmat)(double *x, int n, double *z, precond_type precond_kind);

extern "C" void R_init_RcppKalman(DllInfo *dll) { 
    expmat = (void (*) (double*, int, double*, precond_type)) R_GetCCallable("expm", "expm"); 
} 

//' This function computes the exponential of a matrix.
//'
//' This functions calls the \code{expm} function from the eponymous package 
//' \pkg{expm}. This is implemented via a registered function call, and does
//' not required explicit linking at the C level. However, the \pkg{expm} package
//' is imported in order to access its registered function at the C level.
//' 
//' As the documentation of package \pkg{expm} states, the underlying implementation
//' borrows from the \pkg{Matrix} package which itself takes it from GNU Octave.
//' @title Compute the exponential of a matrix
//' @param x An numeric matrix
//' @return A numeric matrix
//' @seealso The \pkg{expm} package and its documentation.
//' @author Dirk Eddelbuettel
//' @examples 
//' 
//' ## example is from the vignette in package expm
//' M <- matrix(c(4, 1, 1, 2, 4, 1, 0, 1, 4), 3, 3)
//' 
//' ## expected output
//' expM <- matrix(c(147.8666, 127.7811, 127.7811, 183.7651, 183.7651, 163.6796, 71.79703, 91.88257, 111.96811), 3, 3)
//' 
//' ## we only have the expected result to about six digits
//' all.equal(expm(M), expM, tolerance=1.0e-6)
// [[Rcpp::export]]
arma::mat expm(arma::mat x) {
    arma::mat z(x.n_rows, x.n_cols);
    (*expmat)(x.begin(), x.n_rows, z.begin(), Ward_2);
    return z;
}
