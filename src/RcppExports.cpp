// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_kmeans
List rcpp_kmeans(NumericMatrix x, int nstart);
RcppExport SEXP CASC_rcpp_kmeans(SEXP xSEXP, SEXP nstartSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nstart(nstartSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_kmeans(x, nstart));
    return rcpp_result_gen;
END_RCPP
}