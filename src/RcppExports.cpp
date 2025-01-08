// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// c_assert_surv
bool c_assert_surv(const NumericMatrix& mat);
RcppExport SEXP _mlr3proba_c_assert_surv(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(c_assert_surv(mat));
    return rcpp_result_gen;
END_RCPP
}
// c_score_intslogloss
NumericMatrix c_score_intslogloss(const NumericVector& truth, const NumericVector& unique_times, const NumericMatrix& cdf, double eps);
RcppExport SEXP _mlr3proba_c_score_intslogloss(SEXP truthSEXP, SEXP unique_timesSEXP, SEXP cdfSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type unique_times(unique_timesSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_score_intslogloss(truth, unique_times, cdf, eps));
    return rcpp_result_gen;
END_RCPP
}
// c_score_graf_schmid
NumericMatrix c_score_graf_schmid(const NumericVector& truth, const NumericVector& unique_times, const NumericMatrix& cdf, int power);
RcppExport SEXP _mlr3proba_c_score_graf_schmid(SEXP truthSEXP, SEXP unique_timesSEXP, SEXP cdfSEXP, SEXP powerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type unique_times(unique_timesSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< int >::type power(powerSEXP);
    rcpp_result_gen = Rcpp::wrap(c_score_graf_schmid(truth, unique_times, cdf, power));
    return rcpp_result_gen;
END_RCPP
}
// c_weight_survival_score
NumericMatrix c_weight_survival_score(const NumericMatrix& score, const NumericMatrix& truth, const NumericVector& unique_times, const NumericMatrix& cens, bool proper, double eps);
RcppExport SEXP _mlr3proba_c_weight_survival_score(SEXP scoreSEXP, SEXP truthSEXP, SEXP unique_timesSEXP, SEXP censSEXP, SEXP properSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type truth(truthSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type unique_times(unique_timesSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cens(censSEXP);
    Rcpp::traits::input_parameter< bool >::type proper(properSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(c_weight_survival_score(score, truth, unique_times, cens, proper, eps));
    return rcpp_result_gen;
END_RCPP
}
// c_concordance
float c_concordance(const NumericVector& time, const NumericVector& status, const NumericVector& crank, double t_max, const std::string& weight_meth, const NumericMatrix& cens, const NumericMatrix& surv, float tiex);
RcppExport SEXP _mlr3proba_c_concordance(SEXP timeSEXP, SEXP statusSEXP, SEXP crankSEXP, SEXP t_maxSEXP, SEXP weight_methSEXP, SEXP censSEXP, SEXP survSEXP, SEXP tiexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type crank(crankSEXP);
    Rcpp::traits::input_parameter< double >::type t_max(t_maxSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type weight_meth(weight_methSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type cens(censSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type surv(survSEXP);
    Rcpp::traits::input_parameter< float >::type tiex(tiexSEXP);
    rcpp_result_gen = Rcpp::wrap(c_concordance(time, status, crank, t_max, weight_meth, cens, surv, tiex));
    return rcpp_result_gen;
END_RCPP
}
// c_gonen
double c_gonen(const NumericVector& crank, float tiex);
RcppExport SEXP _mlr3proba_c_gonen(SEXP crankSEXP, SEXP tiexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type crank(crankSEXP);
    Rcpp::traits::input_parameter< float >::type tiex(tiexSEXP);
    rcpp_result_gen = Rcpp::wrap(c_gonen(crank, tiex));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mlr3proba_c_assert_surv", (DL_FUNC) &_mlr3proba_c_assert_surv, 1},
    {"_mlr3proba_c_score_intslogloss", (DL_FUNC) &_mlr3proba_c_score_intslogloss, 4},
    {"_mlr3proba_c_score_graf_schmid", (DL_FUNC) &_mlr3proba_c_score_graf_schmid, 4},
    {"_mlr3proba_c_weight_survival_score", (DL_FUNC) &_mlr3proba_c_weight_survival_score, 6},
    {"_mlr3proba_c_concordance", (DL_FUNC) &_mlr3proba_c_concordance, 8},
    {"_mlr3proba_c_gonen", (DL_FUNC) &_mlr3proba_c_gonen, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_mlr3proba(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
