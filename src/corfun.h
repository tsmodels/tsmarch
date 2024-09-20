#ifndef CORFUN_H
#define CORFUN_H

#include <RcppArmadillo.h>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <limits>

arma::mat matrix_sign(arma::mat x);
arma::vec rank_vector(const arma::vec& x);
arma::mat spearman_correlation(const arma::mat& data);
uint64_t insertion_sort(double* arr, size_t len);
uint64_t merge(double* from, double* to, size_t middle, size_t len);
uint64_t merge_sort(double* x, double* buf, size_t len);
uint64_t tied_pairs(double* data, size_t len);
double kendall_tau(arma::vec x, arma::vec y);
arma::uvec sort_index(const arma::vec& v);
arma::mat kendall_correlation(const arma::mat data);
arma::mat transform_correlation(const arma::mat correlation_matrix, Rcpp::String method);
arma::mat make_correlation(const arma::mat data, Rcpp::String method);
arma::mat make_psd(arma::mat x, const double eig_tol = 1e-6, const double conv_tol = 1e-7, const double posd_tol = 1e-8, const int maxit = 100);
arma::mat p2P(arma::vec values, const int m);
bool is_psd(arma::mat x);
arma::mat array_mean(arma::cube x);
#endif // CORFUN_H
