// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// expm
arma::mat expm(arma::mat x);
RcppExport SEXP RcppKalman_expm(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    __result = Rcpp::wrap(expm(x));
    return __result;
END_RCPP
}
// kfPredict
Rcpp::List kfPredict(const arma::vec& x, const arma::mat& P, const arma::mat& A, const arma::mat& Q, const arma::mat& B, const arma::vec& u);
RcppExport SEXP RcppKalman_kfPredict(SEXP xSEXP, SEXP PSEXP, SEXP ASEXP, SEXP QSEXP, SEXP BSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u(uSEXP);
    __result = Rcpp::wrap(kfPredict(x, P, A, Q, B, u));
    return __result;
END_RCPP
}
// kfUpdate
Rcpp::List kfUpdate(const arma::vec& x, const arma::mat& P, const arma::vec& y, const arma::mat& H, const arma::mat& R);
RcppExport SEXP RcppKalman_kfUpdate(SEXP xSEXP, SEXP PSEXP, SEXP ySEXP, SEXP HSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    __result = Rcpp::wrap(kfUpdate(x, P, y, H, R));
    return __result;
END_RCPP
}
// ltiDisc
Rcpp::List ltiDisc(const arma::mat& F, const arma::mat& L, const arma::mat& Q, const double dt);
RcppExport SEXP RcppKalman_ltiDisc(SEXP FSEXP, SEXP LSEXP, SEXP QSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const double >::type dt(dtSEXP);
    __result = Rcpp::wrap(ltiDisc(F, L, Q, dt));
    return __result;
END_RCPP
}
// rtsSmoother
Rcpp::List rtsSmoother(const arma::mat& M, const arma::cube& P, const arma::mat& A, const arma::mat& Q);
RcppExport SEXP RcppKalman_rtsSmoother(SEXP MSEXP, SEXP PSEXP, SEXP ASEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    __result = Rcpp::wrap(rtsSmoother(M, P, A, Q));
    return __result;
END_RCPP
}
// tfSmoother
Rcpp::List tfSmoother(const arma::mat& M, const arma::cube& P, const arma::mat& Y, const arma::mat& A, const arma::mat& Q, const arma::mat& H, const arma::mat& R, const bool useinf);
RcppExport SEXP RcppKalman_tfSmoother(SEXP MSEXP, SEXP PSEXP, SEXP YSEXP, SEXP ASEXP, SEXP QSEXP, SEXP HSEXP, SEXP RSEXP, SEXP useinfSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type H(HSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const bool >::type useinf(useinfSEXP);
    __result = Rcpp::wrap(tfSmoother(M, P, Y, A, Q, H, R, useinf));
    return __result;
END_RCPP
}
