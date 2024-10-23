#ifndef GOGARCH_H
#define GOGARCH_H
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>

static inline arma::vec lower_triangular(const arma::mat M, const int diagonal);
arma::mat coskew_sigma(const arma::vec sigmas);
arma::mat coskewness_block(const arma::rowvec skew);

arma::mat gogarch_covariance(const arma::mat V, const arma::mat A);
arma::mat gogarch_correlation(const arma::mat V, const arma::mat A);
struct gogarch_coskewness_worker : public RcppParallel::Worker {
    const arma::mat& S;
    const arma::mat& A;
    const arma::mat& V;
    const arma::mat& kronA;
    bool standardize;
    arma::cube& result;
    gogarch_coskewness_worker(const arma::mat& S, const arma::mat& A, const arma::mat& V,
                          const arma::mat& kronA, bool standardize, arma::cube& result);
    void operator()(std::size_t begin, std::size_t end);
};
arma::cube gogarch_coskewness(const arma::mat A, const arma::mat S, const arma::mat V, bool standardize);
arma::umat combn(const arma::uvec n, int m);
arma::field<arma::umat> cokurtosis_pairs(int n);
arma::mat cokurtosis_block(const arma::vec s, const arma::vec values);
arma::mat cokurtosis_sigma(const arma::vec sigmas);
struct gogarch_cokurtosis_worker : public RcppParallel::Worker {
    const arma::mat& K;
    const arma::mat& A;
    const arma::mat& V;
    const arma::mat& kronA;
    bool standardize;
    arma::cube& result;
    gogarch_cokurtosis_worker(const arma::mat& K, const arma::mat& A, const arma::mat& V,
                          const arma::mat& kronA, bool standardize, arma::cube& result);
    void operator()(std::size_t begin, std::size_t end);
};
arma::cube gogarch_cokurtosis(const arma::mat A, const arma::mat K, const arma::mat V, bool standardize);

struct gogarch_coskewness_weighted_worker : public RcppParallel::Worker {
    const arma::mat& S;
    const arma::mat& A;
    const arma::mat& kronA;
    const arma::mat& W;
    arma::vec& result;
    gogarch_coskewness_weighted_worker(const arma::mat& S, const arma::mat& A,
                                       const arma::mat& kronA, const arma::mat& W,
                                       arma::vec& result);
    void operator()(std::size_t begin, std::size_t end);
};

arma::vec gogarch_skewness_weighted(const arma::mat A, const arma::mat S, const arma::mat w);


struct gogarch_cokurtosis_weighted_worker : public RcppParallel::Worker {
    const arma::mat& K;
    const arma::mat& V;
    const arma::mat& A;
    const arma::mat& kronA;
    const arma::mat& W;
    arma::vec& result;
    gogarch_cokurtosis_weighted_worker(const arma::mat& K, const arma::mat& V, const arma::mat& A,
                                       const arma::mat& kronA, const arma::mat& W,
                                       arma::vec& result);

    void operator()(std::size_t begin, std::size_t end);
};
arma::vec gogarch_cokurtosis_weighted(const arma::mat A, const arma::mat K, const arma::mat V, const arma::mat w);
arma::vec gogarch_covariance_weighted(const arma::mat V, const arma::mat A, const arma::mat w);

struct gogarch_cokurtosis_weighted_worker_sim : public RcppParallel::Worker {
  // Input data
  const arma::cube& sig;
  const arma::mat& ku;
  const arma::mat& A;
  const arma::mat& kronA;
  const arma::mat& W;
  arma::mat& result;
  int n;
  int nsim;
  gogarch_cokurtosis_weighted_worker_sim(const arma::cube& sig, const arma::mat& ku, 
                                          const arma::mat& A, const arma::mat& kronA, 
                                          const arma::mat& W, arma::mat& result, 
                                          int n, int nsim);
  void operator()(std::size_t begin, std::size_t end);
};

arma::mat gogarch_cokurtosis_weighted_sim(const arma::mat& A, const arma::cube& sig, 
                                           const arma::mat& ku, const arma::mat& weights, int nsim, int n);
#endif // GOGARCH_H
