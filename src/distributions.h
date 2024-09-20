#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
#include <R.h>
#include <RcppArmadillo.h>
#include <RcppBessel.h>
#include <complex>
#include <cmath>

arma::mat rmvt(arma::mat R, arma::mat Z, const double nu);
arma::rowvec rmvt(arma::mat R, arma::rowvec Z, const double nu, const double rc);
arma::mat rmvnorm(arma::mat R, arma::mat Z);
arma::rowvec rmvnorm(arma::mat R, arma::rowvec Z);
arma::mat mpnorm(const arma::mat& x);
arma::mat mpstd(const arma::mat& x, const double shape);
arma::uvec interval(const arma::vec& x, const arma::vec& y);
arma::vec interpolate_window(const arma::vec& x, const arma::vec& y, const arma::vec& z, int w);
arma::cx_vec nigmvcf(const arma::vec& z, const arma::vec& alpha, const arma::vec& beta, const arma::vec& delta, const arma::vec& mu);
arma::vec cfinvnig(const arma::vec& z, double step, const arma::vec& alpha, const arma::vec& beta, const arma::vec& delta, const arma::vec& mu);
std::complex<double> ghypfn(double lambda, double alpha, double beta, double delta, double z);
arma::cx_vec ghypmvcf(const arma::vec& z, const arma::vec& lambda, const arma::vec& alpha,
                      const arma::vec& beta, const arma::vec& delta, const arma::vec& mu);
arma::vec cfinvghyp(const arma::vec& z, double step, const arma::vec& lambda, const arma::vec& alpha, const arma::vec& beta, const arma::vec& delta, const arma::vec& mu);
#endif // DISTRIBUTIONS_H
