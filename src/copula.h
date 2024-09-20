#ifndef COPULA_H
#define COPULA_H

#include "corfun.h"
#include "helpers.h"
#include "distributions.h"
// RcppArmadillo already declared in helpers.h
// Function declarations
inline arma::vec vdstd(arma::vec x, double shape, bool give_log);
inline arma::mat mdstd(arma::mat x, double shape);
inline arma::vec vqstd(arma::vec p, double shape);
inline arma::vec vqnorm(arma::vec p);

// Exported main functions
Rcpp::List copula_constant_normal(const arma::mat u, Rcpp::String method);
Rcpp::List copula_constant_student(double shape, const arma::mat u);
Rcpp::List copula_dynamic_normal(const arma::vec alpha, const arma::vec gamma, const arma::vec beta,
                                 const arma::mat u, Rcpp::IntegerVector dccorder);
Rcpp::List copula_dynamic_student(const arma::vec alpha, const arma::vec gamma, const arma::vec beta,
                                  double shape, const arma::mat u, Rcpp::IntegerVector dccorder);
Rcpp::List copula_constant_normal_filter(const arma::mat u, Rcpp::String method, const int n_update);
Rcpp::List copula_constant_student_filter(double shape, const arma::mat u, const int n_update);
Rcpp::List copula_dynamic_normal_filter(const arma::vec alpha, const arma::vec gamma, const arma::vec beta,
                                        const arma::mat u, Rcpp::IntegerVector dccorder, const int n_update);
Rcpp::List copula_dynamic_student_filter(const arma::vec alpha, const arma::vec gamma,
                                         const arma::vec beta, double shape,
                                         const arma::mat u, Rcpp::IntegerVector dccorder,
                                         const int n_update);
Rcpp::List copula_dynamic_simulate(const arma::vec alpha, const arma::vec gamma,
                                   const arma::vec beta, double shape,
                                   const arma::mat Qbar, arma::mat Nbar,
                                   const arma::cube Qinit,
                                   const arma::mat Zinit,
                                   const arma::mat std_noise,
                                   const int timesteps,
                                   const int burn,
                                   Rcpp::IntegerVector dccorder,
                                   Rcpp::String distribution);
Rcpp::List copula_constant_simulate(const double shape, const arma::mat R,
                                    const arma::mat std_noise, const int timesteps,
                                    Rcpp::String distribution);

// Exported helper functions
arma::mat pit_transform(arma::mat u, double shape, Rcpp::String distribution);
double copula_adcc_constraint(const arma::vec alpha, const arma::vec gamma, const arma::vec beta,
                       double shape, const arma::mat u, Rcpp::IntegerVector dccorder, Rcpp::String distribution);
#endif // COPULA_H
