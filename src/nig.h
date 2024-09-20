#ifndef NIG_H
#define NIG_H

#include <Rcpp.h>
#include <cmath>
#include <vector>

void heap_sort(int n, std::vector<double>& x, std::vector<int>& order);
double besselk1(double x);
void dnig(std::vector<double>& x, double mu, double delta, double alpha, double beta, std::vector<double>& d);
double fdnig(double x, double mu, double delta, double alpha, double beta);
void intdei(double a, double mu, double delta, double alpha, double beta, double *i, double *err);
void pnig(std::vector<double>& x, double mu, double delta, double alpha, double beta, std::vector<double>& p);
double fpnig(double x, double mu, double delta, double alpha, double beta, double sp);
double zbrent(double x1, double x2, double mu, double delta, double alpha, double beta, double sp);
Rcpp::NumericVector qnig(Rcpp::NumericVector p, double mu, double delta, double alpha, double beta);



#endif // NIG_H
