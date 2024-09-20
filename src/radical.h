#ifndef RADICAL_H
#define RADICAL_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

struct radicalrot;

static inline void calculate_rotation(arma::mat &rotate, double theta);
inline arma::mat repmat(const arma::mat &a, const int n, const int m);
struct radical_rotate_worker : public RcppParallel::Worker {
    const RcppParallel::RMatrix<double> X;
    const int m;
    const int k;
    const int nr;
    const int nc;
    RcppParallel::RVector<double> entropy;
    radical_rotate_worker(const Rcpp::NumericMatrix &X, int m, int k, int nr, int nc, Rcpp::NumericVector entropy);
    void operator()(std::size_t begin, std::size_t end);
};

arma::vec radical_rotate(const arma::mat &X, const int m, const int k, const int nr, const int nc);
double rotation_fun(const arma::mat &X, const double sigma, const int m,
                    const int replications, const int k,
                    const double range, const int d, const int N, bool trace);
Rcpp::List radical_recursion(const int k, const double sigma, const double samples,
                             const int replications, const arma::mat whitening_signal,
                             const arma::mat whitening_matrix, const arma::mat dewhitening_matrix,
                             const arma::mat mixed_signal, const arma::vec mixed_mean, bool trace);
#endif // RADICAL_H
