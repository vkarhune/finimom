// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// arma_setdiff
inline arma::vec arma_setdiff(arma::vec x, arma::vec y);
RcppExport SEXP _finimom_arma_setdiff(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(arma_setdiff(x, y));
    return rcpp_result_gen;
END_RCPP
}
// subset_vector
inline arma::vec subset_vector(arma::vec x, arma::uvec pos);
RcppExport SEXP _finimom_subset_vector(SEXP xSEXP, SEXP posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    rcpp_result_gen = Rcpp::wrap(subset_vector(x, pos));
    return rcpp_result_gen;
END_RCPP
}
// set_vector_vals
inline arma::vec set_vector_vals(arma::vec x, arma::uvec pos, arma::vec vals);
RcppExport SEXP _finimom_set_vector_vals(SEXP xSEXP, SEXP posSEXP, SEXP valsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pos(posSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vals(valsSEXP);
    rcpp_result_gen = Rcpp::wrap(set_vector_vals(x, pos, vals));
    return rcpp_result_gen;
END_RCPP
}
// LMarlik
double LMarlik(arma::vec beta, arma::mat sematinv, arma::vec tau, double psi, double r, int d, arma::mat LDmat, double gval);
RcppExport SEXP _finimom_LMarlik(SEXP betaSEXP, SEXP sematinvSEXP, SEXP tauSEXP, SEXP psiSEXP, SEXP rSEXP, SEXP dSEXP, SEXP LDmatSEXP, SEXP gvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sematinv(sematinvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LDmat(LDmatSEXP);
    Rcpp::traits::input_parameter< double >::type gval(gvalSEXP);
    rcpp_result_gen = Rcpp::wrap(LMarlik(beta, sematinv, tau, psi, r, d, LDmat, gval));
    return rcpp_result_gen;
END_RCPP
}
// gf
double gf(arma::vec x, arma::vec z, arma::mat sematinv, arma::mat LDmat, arma::vec tau, int k, int r);
RcppExport SEXP _finimom_gf(SEXP xSEXP, SEXP zSEXP, SEXP sematinvSEXP, SEXP LDmatSEXP, SEXP tauSEXP, SEXP kSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sematinv(sematinvSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LDmat(LDmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(gf(x, z, sematinv, LDmat, tau, k, r));
    return rcpp_result_gen;
END_RCPP
}
// opt_nm
Rcpp::List opt_nm(arma::vec& initval, Rcpp::List& pars);
RcppExport SEXP _finimom_opt_nm(SEXP initvalSEXP, SEXP parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type initval(initvalSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type pars(parsSEXP);
    rcpp_result_gen = Rcpp::wrap(opt_nm(initval, pars));
    return rcpp_result_gen;
END_RCPP
}
// LMarlikApprox
double LMarlikApprox(arma::vec beta, arma::mat sematinv, arma::vec z, arma::vec tau, double psi, double r, int d, arma::mat LDmat, double gval);
RcppExport SEXP _finimom_LMarlikApprox(SEXP betaSEXP, SEXP sematinvSEXP, SEXP zSEXP, SEXP tauSEXP, SEXP psiSEXP, SEXP rSEXP, SEXP dSEXP, SEXP LDmatSEXP, SEXP gvalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sematinv(sematinvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type LDmat(LDmatSEXP);
    Rcpp::traits::input_parameter< double >::type gval(gvalSEXP);
    rcpp_result_gen = Rcpp::wrap(LMarlikApprox(beta, sematinv, z, tau, psi, r, d, LDmat, gval));
    return rcpp_result_gen;
END_RCPP
}
// posterior
Rcpp::List posterior(Rcpp::List dat, arma::vec tau, int maxsize, double r, int p, int niter, arma::vec lpriorval, int approx);
RcppExport SEXP _finimom_posterior(SEXP datSEXP, SEXP tauSEXP, SEXP maxsizeSEXP, SEXP rSEXP, SEXP pSEXP, SEXP niterSEXP, SEXP lpriorvalSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type dat(datSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type maxsize(maxsizeSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lpriorval(lpriorvalSEXP);
    Rcpp::traits::input_parameter< int >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior(dat, tau, maxsize, r, p, niter, lpriorval, approx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_finimom_arma_setdiff", (DL_FUNC) &_finimom_arma_setdiff, 2},
    {"_finimom_subset_vector", (DL_FUNC) &_finimom_subset_vector, 2},
    {"_finimom_set_vector_vals", (DL_FUNC) &_finimom_set_vector_vals, 3},
    {"_finimom_LMarlik", (DL_FUNC) &_finimom_LMarlik, 8},
    {"_finimom_gf", (DL_FUNC) &_finimom_gf, 7},
    {"_finimom_opt_nm", (DL_FUNC) &_finimom_opt_nm, 2},
    {"_finimom_LMarlikApprox", (DL_FUNC) &_finimom_LMarlikApprox, 9},
    {"_finimom_posterior", (DL_FUNC) &_finimom_posterior, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_finimom(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
