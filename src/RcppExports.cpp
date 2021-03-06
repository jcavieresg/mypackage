// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// tps
arma::mat tps(int ep, arma::mat r);
RcppExport SEXP _mypackage_tps(SEXP epSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ep(epSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(tps(ep, r));
    return rcpp_result_gen;
END_RCPP
}
// baseSVD
arma::vec baseSVD(const arma::mat& X);
RcppExport SEXP _mypackage_baseSVD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(baseSVD(X));
    return rcpp_result_gen;
END_RCPP
}
// dcSVD
arma::vec dcSVD(const arma::mat& X);
RcppExport SEXP _mypackage_dcSVD(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(dcSVD(X));
    return rcpp_result_gen;
END_RCPP
}
// f1
arma::mat f1(arma::mat x, arma::mat y);
RcppExport SEXP _mypackage_f1(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f1(x, y));
    return rcpp_result_gen;
END_RCPP
}
// f2
arma::mat f2(arma::mat x, arma::mat y);
RcppExport SEXP _mypackage_f2(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f2(x, y));
    return rcpp_result_gen;
END_RCPP
}
// f3
arma::mat f3(arma::mat x, arma::mat y);
RcppExport SEXP _mypackage_f3(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f3(x, y));
    return rcpp_result_gen;
END_RCPP
}
// f4
arma::mat f4(arma::mat x, arma::mat y);
RcppExport SEXP _mypackage_f4(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(f4(x, y));
    return rcpp_result_gen;
END_RCPP
}
// testfunction
arma::mat testfunction(arma::mat x, arma::mat y);
RcppExport SEXP _mypackage_testfunction(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(testfunction(x, y));
    return rcpp_result_gen;
END_RCPP
}
// meshgrid
void meshgrid(arma::mat& x, arma::mat& y, arma::vec& xv, arma::vec& yv);
RcppExport SEXP _mypackage_meshgrid(SEXP xSEXP, SEXP ySEXP, SEXP xvSEXP, SEXP yvSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type xv(xvSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type yv(yvSEXP);
    meshgrid(x, y, xv, yv);
    return R_NilValue;
END_RCPP
}
// radialFunction
arma::mat radialFunction(arma::mat& r, const int RBFtype, const double R);
RcppExport SEXP _mypackage_radialFunction(SEXP rSEXP, SEXP RBFtypeSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type r(rSEXP);
    Rcpp::traits::input_parameter< const int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(radialFunction(r, RBFtype, R));
    return rcpp_result_gen;
END_RCPP
}
// DistanceMatrix
arma::mat DistanceMatrix(const arma::mat dsites, const arma::mat ctrs);
RcppExport SEXP _mypackage_DistanceMatrix(SEXP dsitesSEXP, SEXP ctrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ctrs(ctrsSEXP);
    rcpp_result_gen = Rcpp::wrap(DistanceMatrix(dsites, ctrs));
    return rcpp_result_gen;
END_RCPP
}
// DistanceMatrix2
arma::mat DistanceMatrix2(const arma::mat dsites, const arma::mat ctrs);
RcppExport SEXP _mypackage_DistanceMatrix2(SEXP dsitesSEXP, SEXP ctrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type ctrs(ctrsSEXP);
    rcpp_result_gen = Rcpp::wrap(DistanceMatrix2(dsites, ctrs));
    return rcpp_result_gen;
END_RCPP
}
// RBF_LS
Rcpp::List RBF_LS(arma::mat dsites, arma::mat ctrs, const int RBFtype, const double R, const int neval);
RcppExport SEXP _mypackage_RBF_LS(SEXP dsitesSEXP, SEXP ctrsSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP nevalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ctrs(ctrsSEXP);
    Rcpp::traits::input_parameter< const int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type neval(nevalSEXP);
    rcpp_result_gen = Rcpp::wrap(RBF_LS(dsites, ctrs, RBFtype, R, neval));
    return rcpp_result_gen;
END_RCPP
}
// RBF_LS2
Rcpp::List RBF_LS2(arma::vec y, arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval);
RcppExport SEXP _mypackage_RBF_LS2(SEXP ySEXP, SEXP dsitesSEXP, SEXP ctrsSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP nevalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ctrs(ctrsSEXP);
    Rcpp::traits::input_parameter< int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type neval(nevalSEXP);
    rcpp_result_gen = Rcpp::wrap(RBF_LS2(y, dsites, ctrs, RBFtype, R, neval));
    return rcpp_result_gen;
END_RCPP
}
// RBF_LSP
Rcpp::List RBF_LSP(arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval, const int omega);
RcppExport SEXP _mypackage_RBF_LSP(SEXP dsitesSEXP, SEXP ctrsSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP nevalSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ctrs(ctrsSEXP);
    Rcpp::traits::input_parameter< int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type neval(nevalSEXP);
    Rcpp::traits::input_parameter< const int >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(RBF_LSP(dsites, ctrs, RBFtype, R, neval, omega));
    return rcpp_result_gen;
END_RCPP
}
// RBF_LSP2
Rcpp::List RBF_LSP2(arma::vec y, arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval, const int omega);
RcppExport SEXP _mypackage_RBF_LSP2(SEXP ySEXP, SEXP dsitesSEXP, SEXP ctrsSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP nevalSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ctrs(ctrsSEXP);
    Rcpp::traits::input_parameter< int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type neval(nevalSEXP);
    Rcpp::traits::input_parameter< const int >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(RBF_LSP2(y, dsites, ctrs, RBFtype, R, neval, omega));
    return rcpp_result_gen;
END_RCPP
}
// RBF_LSPSVD
Rcpp::List RBF_LSPSVD(arma::vec y, arma::mat dsites, arma::mat ctrs, int RBFtype, const double R, const int neval, const int omega);
RcppExport SEXP _mypackage_RBF_LSPSVD(SEXP ySEXP, SEXP dsitesSEXP, SEXP ctrsSEXP, SEXP RBFtypeSEXP, SEXP RSEXP, SEXP nevalSEXP, SEXP omegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dsites(dsitesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ctrs(ctrsSEXP);
    Rcpp::traits::input_parameter< int >::type RBFtype(RBFtypeSEXP);
    Rcpp::traits::input_parameter< const double >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type neval(nevalSEXP);
    Rcpp::traits::input_parameter< const int >::type omega(omegaSEXP);
    rcpp_result_gen = Rcpp::wrap(RBF_LSPSVD(y, dsites, ctrs, RBFtype, R, neval, omega));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mypackage_tps", (DL_FUNC) &_mypackage_tps, 2},
    {"_mypackage_baseSVD", (DL_FUNC) &_mypackage_baseSVD, 1},
    {"_mypackage_dcSVD", (DL_FUNC) &_mypackage_dcSVD, 1},
    {"_mypackage_f1", (DL_FUNC) &_mypackage_f1, 2},
    {"_mypackage_f2", (DL_FUNC) &_mypackage_f2, 2},
    {"_mypackage_f3", (DL_FUNC) &_mypackage_f3, 2},
    {"_mypackage_f4", (DL_FUNC) &_mypackage_f4, 2},
    {"_mypackage_testfunction", (DL_FUNC) &_mypackage_testfunction, 2},
    {"_mypackage_meshgrid", (DL_FUNC) &_mypackage_meshgrid, 4},
    {"_mypackage_radialFunction", (DL_FUNC) &_mypackage_radialFunction, 3},
    {"_mypackage_DistanceMatrix", (DL_FUNC) &_mypackage_DistanceMatrix, 2},
    {"_mypackage_DistanceMatrix2", (DL_FUNC) &_mypackage_DistanceMatrix2, 2},
    {"_mypackage_RBF_LS", (DL_FUNC) &_mypackage_RBF_LS, 5},
    {"_mypackage_RBF_LS2", (DL_FUNC) &_mypackage_RBF_LS2, 6},
    {"_mypackage_RBF_LSP", (DL_FUNC) &_mypackage_RBF_LSP, 6},
    {"_mypackage_RBF_LSP2", (DL_FUNC) &_mypackage_RBF_LSP2, 7},
    {"_mypackage_RBF_LSPSVD", (DL_FUNC) &_mypackage_RBF_LSPSVD, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_mypackage(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
