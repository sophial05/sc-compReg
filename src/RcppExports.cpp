// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// compReg
Rcpp::List compReg(const arma::sp_mat& O1, const arma::sp_mat& E1, arma::uvec O1Idx, arma::uvec E1Idx, std::vector<std::string> symbol1, std::vector<std::string> peakName1, const arma::sp_mat& O2, const arma::sp_mat& E2, arma::uvec O2Idx, arma::uvec E2Idx, std::vector<std::string> symbol2, std::vector<std::string> peakName2, std::vector<std::string> peakNameIntersect);
RcppExport SEXP _scCompReg_compReg(SEXP O1SEXP, SEXP E1SEXP, SEXP O1IdxSEXP, SEXP E1IdxSEXP, SEXP symbol1SEXP, SEXP peakName1SEXP, SEXP O2SEXP, SEXP E2SEXP, SEXP O2IdxSEXP, SEXP E2IdxSEXP, SEXP symbol2SEXP, SEXP peakName2SEXP, SEXP peakNameIntersectSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type O1(O1SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type E1(E1SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type O1Idx(O1IdxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type E1Idx(E1IdxSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type symbol1(symbol1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type peakName1(peakName1SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type O2(O2SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type E2(E2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type O2Idx(O2IdxSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type E2Idx(E2IdxSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type symbol2(symbol2SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type peakName2(peakName2SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type peakNameIntersect(peakNameIntersectSEXP);
    rcpp_result_gen = Rcpp::wrap(compReg(O1, E1, O1Idx, E1Idx, symbol1, peakName1, O2, E2, O2Idx, E2Idx, symbol2, peakName2, peakNameIntersect));
    return rcpp_result_gen;
END_RCPP
}
// initializeMatrix
Rcpp::List initializeMatrix(const unsigned int POnRow, const unsigned int POnCol, const unsigned int XnCol, const unsigned int k, const arma::sp_mat& D);
RcppExport SEXP _scCompReg_initializeMatrix(SEXP POnRowSEXP, SEXP POnColSEXP, SEXP XnColSEXP, SEXP kSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type POnRow(POnRowSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type POnCol(POnColSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type XnCol(XnColSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(initializeMatrix(POnRow, POnCol, XnCol, k, D));
    return rcpp_result_gen;
END_RCPP
}
// computeLambda
Rcpp::List computeLambda(const arma::sp_mat& PeakO, const arma::mat& w1, const arma::mat& h1, const arma::sp_mat& X, const arma::mat& w2, const arma::mat& h2, const arma::sp_mat& D, double beta, double alpha, double eps);
RcppExport SEXP _scCompReg_computeLambda(SEXP PeakOSEXP, SEXP w1SEXP, SEXP h1SEXP, SEXP XSEXP, SEXP w2SEXP, SEXP h2SEXP, SEXP DSEXP, SEXP betaSEXP, SEXP alphaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type PeakO(PeakOSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type h1(h1SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w2(w2SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type h2(h2SEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeLambda(PeakO, w1, h1, X, w2, h2, D, beta, alpha, eps));
    return rcpp_result_gen;
END_RCPP
}
// iterateCluster
Rcpp::List iterateCluster(const arma::sp_mat& PeakO, const arma::sp_mat& X, const arma::sp_mat& D, const unsigned int k, const unsigned int maxIter, double lambda1, double lambda2, arma::mat W10, arma::mat H10, arma::mat W20, arma::mat H20, double epsD, double tolX, double tolFun, bool verbose, int loopUpdate);
RcppExport SEXP _scCompReg_iterateCluster(SEXP PeakOSEXP, SEXP XSEXP, SEXP DSEXP, SEXP kSEXP, SEXP maxIterSEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP W10SEXP, SEXP H10SEXP, SEXP W20SEXP, SEXP H20SEXP, SEXP epsDSEXP, SEXP tolXSEXP, SEXP tolFunSEXP, SEXP verboseSEXP, SEXP loopUpdateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type PeakO(PeakOSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type maxIter(maxIterSEXP);
    Rcpp::traits::input_parameter< double >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< double >::type lambda2(lambda2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W10(W10SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H10(H10SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W20(W20SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H20(H20SEXP);
    Rcpp::traits::input_parameter< double >::type epsD(epsDSEXP);
    Rcpp::traits::input_parameter< double >::type tolX(tolXSEXP);
    Rcpp::traits::input_parameter< double >::type tolFun(tolFunSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type loopUpdate(loopUpdateSEXP);
    rcpp_result_gen = Rcpp::wrap(iterateCluster(PeakO, X, D, k, maxIter, lambda1, lambda2, W10, H10, W20, H20, epsD, tolX, tolFun, verbose, loopUpdate));
    return rcpp_result_gen;
END_RCPP
}
// postLapMatMult
Rcpp::List postLapMatMult(arma::mat W1, arma::mat W2, arma::mat H1, arma::mat H2);
RcppExport SEXP _scCompReg_postLapMatMult(SEXP W1SEXP, SEXP W2SEXP, SEXP H1SEXP, SEXP H2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type W1(W1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W2(W2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H1(H1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H2(H2SEXP);
    rcpp_result_gen = Rcpp::wrap(postLapMatMult(W1, W2, H1, H2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scCompReg_compReg", (DL_FUNC) &_scCompReg_compReg, 13},
    {"_scCompReg_initializeMatrix", (DL_FUNC) &_scCompReg_initializeMatrix, 5},
    {"_scCompReg_computeLambda", (DL_FUNC) &_scCompReg_computeLambda, 10},
    {"_scCompReg_iterateCluster", (DL_FUNC) &_scCompReg_iterateCluster, 16},
    {"_scCompReg_postLapMatMult", (DL_FUNC) &_scCompReg_postLapMatMult, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_scCompReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}