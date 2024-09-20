// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include "radical.h"

using namespace Rcpp;

struct radicalrot {
    double theta;
    arma::mat rotate;
};

static inline void calculate_rotation(arma::mat &rotate, double theta) {
    rotate(0,0) = cos(theta);
    rotate(0,1) = -1.0 * sin(theta);
    rotate(1,0) = sin(theta);
    rotate(1,1) = cos(theta);
}

inline arma::mat repmat(const arma::mat &a, const int n, const int m) {
    return arma::kron(arma::mat(n, m, arma::fill::ones), a);
}

void calculate_rotation(arma::mat &rotate, double theta);

// Worker for parallel computation
radical_rotate_worker::radical_rotate_worker(const Rcpp::NumericMatrix &X, int m, int k, int nr, int nc, Rcpp::NumericVector entropy)
    : X(X), m(m), k(k), nr(nr), nc(nc), entropy(entropy) {}

void radical_rotate_worker::operator()(std::size_t begin, std::size_t end) {
    arma::mat rotate(2, 2);
    arma::mat rotatep(nr, nc);
    double theta = 0.0;
    arma::rowvec values(nc);
    arma::rowvec marginal_theta(2);
    arma::rowvec tmp;

    for (std::size_t i = begin; i < end; ++i) {
        theta = ((double)i / ((double)k - 1.0) * (0.5 * M_PI) - (0.25 * M_PI));
        rotate.zeros();
        rotatep.zeros();
        calculate_rotation(rotate, theta);
        for (int j = 0; j < nr; ++j) {
            for (int l = 0; l < nc; ++l) {
                rotatep(j, l) = 0.0;
                for (int p = 0; p < 2; ++p) {
                    rotatep(j, l) += rotate(j, p) * X(p, l);
                }
            }
        }
        values = rotatep.row(0);
        values = arma::sort(values);
        tmp = arma::log(values.subvec(m, nc - 1) - values.subvec(0, nc - m - 1));
        marginal_theta[0] = arma::as_scalar(arma::accu(tmp));
        values = rotatep.row(1);
        values = arma::sort(values);
        tmp = arma::log(values.subvec(m, nc - 1) - values.subvec(0, nc - m - 1));
        marginal_theta[1] = arma::as_scalar(arma::accu(tmp));
        entropy[i] = arma::as_scalar(arma::accu(marginal_theta));
    }
};


arma::vec radical_rotate(const arma::mat &X, const int m, const int k, const int nr, const int nc) {
    NumericMatrix X_numeric(X.n_rows, X.n_cols);
    std::copy(X.begin(), X.end(), X_numeric.begin());
    NumericVector entropy(k);
    radical_rotate_worker worker(X_numeric, m, k, nr, nc, entropy);
    RcppParallel::parallelFor(0, k, worker);
    return Rcpp::as<arma::vec>(entropy);
}


double rotation_fun(const arma::mat &X, const double sigma, const int m,
                    const int replications, const int k,
                    const double range, const int d, const int N, bool trace)
{
    arma::mat x_augmented;
    if (replications == 1) {
        x_augmented = X;
    } else {
        arma::mat noise = arma::randn(d, N * replications) * sigma;
        x_augmented = noise + repmat(X, 1, replications);
    }
    const int nr = x_augmented.n_rows;
    const int nc = x_augmented.n_cols;
    arma::vec entropy = radical_rotate(x_augmented, m, k, nr, nc);
    arma::uvec indices = arma::sort_index(entropy);
    int theta_index = indices(0);
    double theta_star = static_cast<double>(theta_index) / (k - 1) * M_PI / 2.0 - M_PI / 4.0;

    if (trace) {
        Rcpp::Rcout << " rotated " << std::round(theta_star / (2 * M_PI) * 360 * 100) / 100 << " degrees." << std::endl;
    }
    return theta_star;
}

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::export(.radical_recursion)]]
Rcpp::List radical_recursion(const int k, const double sigma, const double samples,
                             const int replications, const arma::mat whitening_signal,
                             const arma::mat whitening_matrix, const arma::mat dewhitening_matrix,
                             const arma::mat mixed_signal, const arma::vec mixed_mean, bool trace)
{
    int N = whitening_signal.n_rows;
    int sweeps = N - 1;
    arma::mat total_rotation = arma::eye(N,N);
    arma::mat old_total_rotation = arma::eye(N,N);
    double start_k_float = (double)k / std::pow(1.3, std::ceil(sweeps / 2.0));
    double new_k_float = start_k_float;
    arma::mat x_current = whitening_signal;
    double new_k = 0.0;
    double range = M_PI / 2.0;
    const int spacings = (int)std::floor(std::sqrt(samples));
    arma::rowvec done = arma::ones<arma::rowvec>(samples);
    arma::mat current_sub_space = arma::zeros(2, whitening_signal.n_cols);
    const int d = 2;
    const int nc = whitening_signal.n_cols;
    double theta_star = 0.0;
    for (int sweep_num = 1; sweep_num <= sweeps; ++sweep_num) {
        if (trace) Rcpp::Rcout << "sweep: " << sweep_num << "/" << sweeps << std::endl;
        if (((double)sweep_num > ((double)sweeps / 2.0))) {
            new_k_float = new_k_float * 1.3;
            new_k = std::floor(new_k_float);
        } else {
            new_k_float = start_k_float;
            new_k = std::max(30.0, std::floor(new_k_float));
        }
        for (int i = 0; i < N - 1; ++i) {
            for (int j = i + 1; j < N; ++j) {
                if (trace) Rcpp::Rcout << "unmixing dimensions " << i + 1 << "..." << j + 1 << std::endl;
                current_sub_space.row(0) = x_current.row(i);
                current_sub_space.row(1) = x_current.row(j);
                theta_star = rotation_fun(current_sub_space, sigma, spacings, replications, new_k, range, d, nc, trace);
                arma::mat new_rotation_component = arma::eye<arma::mat>(N, N);
                new_rotation_component.at(i, i) =  cos(theta_star);
                new_rotation_component.at(i, j) = -sin(theta_star);
                new_rotation_component.at(j, i) =  sin(theta_star);
                new_rotation_component.at(j, j) =  cos(theta_star);
                total_rotation = new_rotation_component * total_rotation;
                x_current = total_rotation * whitening_signal;
            }
        }
        old_total_rotation = total_rotation;
    }
    arma::mat W = total_rotation * whitening_matrix;
    arma::mat S = W * mixed_signal + (W * mixed_mean) * done;
    arma::mat U = total_rotation;
    arma::mat A = (arma::inv(U).t() * dewhitening_matrix.t()).t();

    return List::create(
        Named("A") = A.t(),
        Named("W") = W.t(),
        Named("U") = U.t(),
        Named("S") = S.t());
}
