// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;


RcppExport SEXP _rcpp_module_boot_stan_fit4homo_hist_Rt_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_Rt_fixed_k_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_Rt_fixed_k_ahead_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_Rt_homo_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4lgp_Rt_homo_ahead_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4homo_hist_Rt_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4homo_hist_Rt_mod, 0},
    {"_rcpp_module_boot_stan_fit4lgp_Rt_fixed_k_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_Rt_fixed_k_mod, 0},
    {"_rcpp_module_boot_stan_fit4lgp_Rt_fixed_k_ahead_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_Rt_fixed_k_ahead_mod, 0},
    {"_rcpp_module_boot_stan_fit4lgp_Rt_homo_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_Rt_homo_mod, 0},
    {"_rcpp_module_boot_stan_fit4lgp_Rt_homo_ahead_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4lgp_Rt_homo_ahead_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_assessEpidemicCurves(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
