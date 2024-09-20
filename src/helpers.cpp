// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp14)]]
#include "helpers.h"
using namespace Rcpp;

double quadratic_form(const arma::vec& values, const arma::mat& w, int m, arma::uvec& lower_indices, arma::mat& V) {
    V.fill(0.0);
    V.elem(lower_indices) = values;
    V = arma::symmatl(V);
    //std::cout<<w<<std::endl;
    //std::cout<<V<<std::endl;
    double result = arma::as_scalar(w * V * w.t());
    return result;
}

// [[Rcpp::export(.aggregate_mu)]]
arma::mat aggregate_mu(arma::cube mu, arma::mat w)
{
    int s = mu.n_slices;
    int h = mu.slice(0).n_rows;
    arma::mat wmu = arma::zeros(s, h);
    for(int i=0;i<s;i++){
        arma::mat tmp = mu.slice(i);
        for(int j=0;j<h;j++){
            wmu(i,j) = arma::accu(tmp.row(j) % w.row(j));
        }
    }
    return wmu;
}

// [[Rcpp::export(.aggregate_sigma)]]
arma::mat aggregate_sigma(arma::cube sigma, arma::mat w)
{
    int s = sigma.n_slices;
    int m = w.n_cols;
    int h =  sigma.slice(0).n_rows;
    arma::mat wsigma = arma::zeros(s, h);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(m, m));
    arma::mat V(m, m, arma::fill::zeros);  // Pre-allocate matrix V
    for(int i = 0;i<s;++i) {
        for(int j = 0;j<h;j++) {
            arma::vec row_vec = sigma.slice(i).row(j).t();  // Get the row vector
            double result = quadratic_form(row_vec, w.row(j), m, lower_indices, V);
            wsigma(i, j) = result;
        }
    }
    return wsigma;
}

// [[Rcpp::export(.tril2sym)]]
arma::cube tril2sym(arma::mat values, const int m, const bool diag)
{
    const int h = values.n_rows;
    arma::cube S(m, m, h);
    arma::uvec lower_indices;
    if (diag) {
        lower_indices = arma::trimatl_ind(arma::size(m, m));
    } else {
        lower_indices = arma::trimatl_ind(arma::size(m, m), -1);
    }
    arma::mat V(m, m, arma::fill::zeros);
    for(int i=0;i<h;i++){
        V.fill(0.0);
        V.elem(lower_indices) = values.row(i).t();
        if (!diag) {
            V.diag().ones();
        }
        V = arma::symmatl(V);
        S.slice(i) = V;
    }
    return S;
}

// [[Rcpp::export(.sym2tril)]]
arma::mat sym2tril(arma::cube S, const bool diag)
{
    int n_slices = S.n_slices;
    int m = S.slice(0).n_cols;
    arma::uvec lower_indices;
    if (diag) {
        lower_indices = arma::trimatl_ind(arma::size(m,m));
    } else {
        lower_indices = arma::trimatl_ind(arma::size(m,m), -1);
    }
    int n = lower_indices.n_elem;
    arma::mat M = arma::zeros(n_slices, n);
    for(int i=0;i<n_slices;++i) {
        M.row(i) = arma::conv_to<arma::rowvec>::from(S.slice(i).elem(lower_indices));
    }
    return M;
}


// [[Rcpp::export(.generate_constant_covariance)]]
Rcpp::List generate_constant_covariance(const arma::mat& correlation, const arma::mat& sigmas, const arma::mat& residuals) {
    int n = sigmas.n_rows;
    int m = sigmas.n_cols;

    arma::mat whitened_residuals(n, m);
    arma::cube covariance(m, m, n);

    for (int i = 0; i < n; i++) {
        arma::mat covariance_slice = arma::diagmat(sigmas.row(i)) * correlation * arma::diagmat(sigmas.row(i));

        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, covariance_slice);
        arma::mat tmp_z = eigvec * arma::diagmat(1.0/arma::sqrt(eigval)) * eigvec.t();

        whitened_residuals.row(i) = residuals.row(i) * tmp_z;

        covariance.slice(i) = covariance_slice;
    }

    return Rcpp::List::create(
        Rcpp::Named("H") = covariance,
        Rcpp::Named("W") = whitened_residuals
    );
}

// [[Rcpp::export(.generate_dynamic_covariance)]]
Rcpp::List generate_dynamic_covariance(arma::cube correlation, const arma::mat& sigmas, const arma::mat& residuals) {
    int n = sigmas.n_rows;
    int m = sigmas.n_cols;

    arma::mat whitened_residuals(n, m);
    arma::cube covariance(m, m, n);
    for (int i = 0; i < n; i++) {
        arma::mat R = correlation.slice(i);
        arma::mat covariance_slice = arma::diagmat(sigmas.row(i)) * R * arma::diagmat(sigmas.row(i));
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, covariance_slice);
        arma::mat tmp_z = eigvec * arma::diagmat(1.0/arma::sqrt(eigval)) * eigvec.t();
        whitened_residuals.row(i) = residuals.row(i) * tmp_z;
        covariance.slice(i) = covariance_slice;
    }
    return Rcpp::List::create(
        Rcpp::Named("H") = covariance,
        Rcpp::Named("W") = whitened_residuals
    );
}

inline arma::rowvec cor2covin(const arma::rowvec r, const arma::rowvec sigma, const int m, arma::uvec lower_indices) {
    arma::mat C(m, m, arma::fill::zeros);
    C.elem(lower_indices) = r;
    C = arma::symmatl(C);
    C.diag().ones();
    arma::mat s = arma::diagmat(sigma);
    arma::mat V = s * C * s;
    arma::uvec lower_indices_v = arma::trimatl_ind(arma::size(m, m));
    arma::rowvec v = V.elem(lower_indices_v).t();
    return v;
}

cor2cov_worker::cor2cov_worker(const arma::cube& r, const arma::cube& sigma, arma::cube& v, int m, arma::uvec lower_indices)
    : r(r), sigma(sigma), v(v), m(m), lower_indices(lower_indices) {}

void cor2cov_worker::operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
        arma::mat out = arma::zeros(r.slice(i).n_rows, r.slice(i).n_cols + m);
        for (arma::uword j = 0; j < r.slice(i).n_rows; ++j) {
            out.row(j) = cor2covin(r.slice(i).row(j), sigma.slice(i).row(j), m, lower_indices);
        }
        v.slice(i) = out;
    }
}

// [[Rcpp::export(.cor2cov)]]
arma::cube cor2cov(const arma::cube r, const arma::cube sigma, const int m) {
    int slice_n = r.n_slices;
    int mat_n = r.slice(0).n_rows;
    int mat_m = r.slice(0).n_cols;
    int n = mat_m + m;
    arma::cube v(mat_n, n, slice_n, arma::fill::zeros);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(m, m), -1);
    cor2cov_worker worker(r, sigma, v, m, lower_indices);
    parallelFor(0, slice_n, worker);
    return v;
}

cor2cov2_worker::cor2cov2_worker(const arma::rowvec& r, const arma::cube& sigma, arma::cube& v, const int m, arma::uvec lower_indices)
    : r(r), sigma(sigma), v(v), m(m), lower_indices(lower_indices) {}

void cor2cov2_worker::operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
        arma::mat out = arma::zeros(sigma.slice(i).n_rows, r.n_elem + m);
        for (arma::uword j = 0; j < sigma.slice(i).n_rows; ++j) {
            out.row(j) = cor2covin(r, sigma.slice(i).row(j), m, lower_indices);
        }
        v.slice(i) = out;
    }
}

// [[Rcpp::export(.cor2cov2)]]
arma::cube cor2cov2(const arma::rowvec r, const arma::cube sigma, const int m) {
    int slice_n = sigma.n_slices;
    int mat_n = sigma.slice(0).n_rows;
    int mat_m = r.n_elem;
    int n = mat_m + m;
    arma::cube v(mat_n, n, slice_n, arma::fill::zeros);
    arma::uvec lower_indices = arma::trimatl_ind(arma::size(m, m), -1);
    cor2cov2_worker worker(r, sigma, v, m, lower_indices);
    parallelFor(0, slice_n, worker);
    return v;
}
