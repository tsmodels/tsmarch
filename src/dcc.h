#ifndef DCC_H
#define DCC_H

#include "helpers.h"
#include "corfun.h"
#include "distributions.h"
// RcppArmadillo already declared in helpers.h

Rcpp::List dcc_constant_normal(const arma::mat Z, const arma::mat S);
Rcpp::List dcc_constant_student(const arma::mat Z, const arma::mat S, const double shape);
Rcpp::List dcc_dynamic_normal(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder);
Rcpp::List dcc_dynamic_student(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, double shape, const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder);
Rcpp::List dcc_constant_normal_filter(const arma::mat Z, const arma::mat S, const int n_update);
Rcpp::List dcc_constant_student_filter(double shape, const arma::mat Z, const arma::mat S, const int n_update);
Rcpp::List dcc_dynamic_normal_filter(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder, const int n_update);
Rcpp::List dcc_dynamic_student_filter(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, double shape, const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder, const int n_update);
Rcpp::List dcc_dynamic_simulate(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, double shape, const arma::mat Qbar, arma::mat Nbar,
                                const arma::cube Qinit, const arma::mat Zinit, const arma::mat std_noise, const int timesteps, const int burn, Rcpp::IntegerVector dccorder, Rcpp::String distribution);
Rcpp::List dcc_constant_simulate(const double shape, const arma::mat R, const arma::mat std_noise, const int timesteps, Rcpp::String distribution);
double adcc_constraint(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, double shape, const arma::mat z, Rcpp::IntegerVector dccorder);







#endif // DCC_H
