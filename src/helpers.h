#ifndef HELPERS_H
#define HELPERS_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

double quadratic_form(const arma::vec& values, const arma::mat& w, int m, arma::uvec& lower_indices, arma::mat& V);
arma::mat aggregate_mu(arma::cube mu, arma::mat w);
arma::mat aggregate_sigma(arma::cube sigma, arma::mat w);
arma::cube tril2sym(arma::mat values, const int m, const bool diag);
arma::mat sym2tril(arma::cube S, const bool diag);
Rcpp::List generate_constant_covariance(const arma::mat& correlation,
                                        const arma::mat& sigmas,
                                        const arma::mat& residuals);
Rcpp::List generate_dynamic_covariance(arma::cube correlation,
                                       const arma::mat& sigmas,
                                       const arma::mat& residuals);

arma::rowvec cor2covin(const arma::rowvec r, const arma::rowvec sigma, const int m, arma::uvec lower_indices);

struct cor2cov_worker : public RcppParallel::Worker {
    const arma::cube& r;
    const arma::cube& sigma;
    arma::cube& v;
    int m;
    arma::uvec lower_indices;

    cor2cov_worker(const arma::cube& r, const arma::cube& sigma, arma::cube& v, int m, arma::uvec lower_indices);

    void operator()(std::size_t begin, std::size_t end);
};

arma::cube cor2cov(const arma::cube r, const arma::cube sigma, const int m);

struct cor2cov2_worker : public RcppParallel::Worker {
    const arma::rowvec& r;
    const arma::cube& sigma;
    arma::cube& v;
    const int m;
    arma::uvec lower_indices;

    cor2cov2_worker(const arma::rowvec& r, const arma::cube& sigma, arma::cube& v, const int m, arma::uvec lower_indices);

    void operator()(std::size_t begin, std::size_t end);
};

arma::cube cor2cov2(const arma::rowvec r, const arma::cube sigma, const int m);


#endif // HELPERS_H
