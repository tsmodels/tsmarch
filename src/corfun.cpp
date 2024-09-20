// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include "corfun.h"
using namespace Rcpp;


// [[Rcpp::export(.matrix_sign)]]
arma::mat matrix_sign(arma::mat x)
{
    int m = x.n_cols;
    int n = x.n_rows;
    arma::mat s = arma::zeros(n, m);
    for (int i = 0; i < m; i++) {
        arma::vec col_i = s.col(i);
        col_i = (-1.0 * arma::sign(x.col(i)) + 1.0) / 2.0;
        col_i(arma::find(col_i == 0.5)).zeros();
        s.col(i) = col_i;
    }
    return s;
}

arma::vec rank_vector(const arma::vec& x) {
    int n = x.n_elem;
    arma::vec ranks(n);
    arma::vec sorted_vec = x;
    arma::uvec sorted_indices = arma::regspace<arma::uvec>(0, n-1);
    // Sort the values and indices
    arma::uvec sorted_order = arma::stable_sort_index(sorted_vec);
    sorted_vec = sorted_vec.elem(sorted_order);
    sorted_indices = sorted_indices.elem(sorted_order);
    // Assign ranks, handling ties by assigning the average rank
    int i = 0;
    while (i < n) {
        int tie_start = i;
        while (i < n - 1 && sorted_vec[i] == sorted_vec[i + 1]) {
            ++i;
        }
        double rank = (tie_start + i + 2) / 2.0;  // average rank for ties
        for (int j = tie_start; j <= i; ++j) {
            ranks[sorted_indices[j]] = rank;
        }
        ++i;
    }

    return ranks;
}

arma::mat spearman_correlation(const arma::mat& data) {
    int n = data.n_rows;
    int m = data.n_cols;
    // Rank the columns
    arma::mat ranked_data(n, m);
    for (int j = 0; j < m; ++j) {
        ranked_data.col(j) = rank_vector(data.col(j));
    }
    // Compute the correlation matrix
    arma::mat correlation_matrix(m, m);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (i == j) {
                correlation_matrix(i, j) = 1.0;
            } else {
                double rank_covariance = arma::as_scalar((ranked_data.col(i) - arma::mean(ranked_data.col(i))).t() *
                                                         (ranked_data.col(j) - arma::mean(ranked_data.col(j)))) / (n - 1);
                double rank_variance_i = arma::as_scalar(arma::sum(arma::square(ranked_data.col(i) - arma::mean(ranked_data.col(i))))) / (n - 1);
                double rank_variance_j = arma::as_scalar(arma::sum(arma::square(ranked_data.col(j) - arma::mean(ranked_data.col(j))))) / (n - 1);
                correlation_matrix(i, j) = rank_covariance / std::sqrt(rank_variance_i * rank_variance_j);
            }
        }
    }

    return correlation_matrix;
}

uint64_t insertion_sort(double* arr, size_t len) {
    size_t maxJ, i;
    uint64_t swapCount = 0;
    if (len < 2) {
        return 0;
    }
    maxJ = len - 1;
    for (i = len - 2; i < len; --i) {
        size_t j = i;
        double val = arr[i];
        for (; j < maxJ && arr[j + 1] < val; ++j) {
            arr[j] = arr[j + 1];
        }
        arr[j] = val;
        swapCount += (j - i);
    }
    return swapCount;
}

uint64_t merge(double* from, double* to, size_t middle, size_t len) {
    size_t bufIndex, leftLen, rightLen;
    uint64_t swaps;
    double* left;
    double* right;
    bufIndex = 0;
    swaps = 0;
    left = from;
    right = from + middle;
    rightLen = len - middle;
    leftLen = middle;
    while (leftLen && rightLen) {
        if (right[0] < left[0]) {
            to[bufIndex] = right[0];
            swaps += leftLen;
            rightLen--;
            right++;
        } else {
            to[bufIndex] = left[0];
            leftLen--;
            left++;
        }
        bufIndex++;
    }
    if (leftLen) {
        std::copy(left, left + leftLen, to + bufIndex);
    } else if (rightLen) {
        std::copy(right, right + rightLen, to + bufIndex);
    }
    return swaps;
}

uint64_t merge_sort(double* x, double* buf, size_t len) {
    uint64_t swaps;
    size_t half;
    if (len < 10) {
        return insertion_sort(x, len);
    }
    swaps = 0;
    if (len < 2) {
        return 0;
    }
    half = len / 2;
    swaps += merge_sort(x, buf, half);
    swaps += merge_sort(x + half, buf + half, len - half);
    swaps += merge(x, buf, half, len);
    std::copy(buf, buf + len, x);
    return swaps;
}

uint64_t tied_pairs(double* data, size_t len) {
    uint64_t Ms = 0, tieCount = 0;
    size_t i;
    for (i = 1; i < len; i++) {
        if (data[i] == data[i - 1]) {
            tieCount++;
        } else if (tieCount) {
            Ms += (tieCount * (tieCount + 1)) / 2;
            tieCount++;
            tieCount = 0;
        }
    }
    if (tieCount) {
        Ms += (tieCount * (tieCount + 1)) / 2;
        tieCount++;
    }
    return Ms;
}

double kendall_tau(arma::vec x, arma::vec y) {
    size_t len = x.n_elem;
    uint64_t m1 = 0, m2 = 0, tieCount, swapCount, nPair;
    int64_t s;
    size_t i;

    nPair = (uint64_t)len * ((uint64_t)len - 1) / 2;
    s = nPair;
    tieCount = 0;

    double* arr1 = x.memptr();
    double* arr2 = y.memptr();

    for (i = 1; i < len; i++) {
        if (arr1[i - 1] == arr1[i]) {
            tieCount++;
        } else if (tieCount > 0) {
            std::sort(arr2 + i - tieCount - 1, arr2 + i);
            m1 += tieCount * (tieCount + 1) / 2;
            s += tied_pairs(arr2 + i - tieCount - 1, tieCount + 1);
            tieCount++;
            tieCount = 0;
        }
    }
    if (tieCount > 0) {
        std::sort(arr2 + i - tieCount - 1, arr2 + i);
        m1 += tieCount * (tieCount + 1) / 2;
        s += tied_pairs(arr2 + i - tieCount - 1, tieCount + 1);
        tieCount++;
    }
    double* buf = new double[len];
    swapCount = merge_sort(arr2, buf, len);
    delete[] buf;

    m2 = tied_pairs(arr2, len);
    s -= (m1 + m2) + 2 * swapCount;

    double denominator1 = nPair - m1;
    double denominator2 = nPair - m2;
    double cor = s / sqrt(denominator1) / sqrt(denominator2);

    return cor;
}

arma::uvec sort_index(const arma::vec& v) {
    size_t n = v.n_elem;
    arma::uvec indices(n);
    // Initialize indices with values from 0 to n-1
    for (size_t i = 0; i < n; ++i) {
        indices(i) = i;
    }
    // Sort indices based on the corresponding values in v
    std::sort(indices.begin(), indices.end(), [&v](size_t i1, size_t i2) {
        return v(i1) < v(i2);
    });
    return indices;
}

arma::mat kendall_correlation(const arma::mat data) {
    int m = data.n_cols;
    int n = data.n_rows;
    arma::mat correlation_matrix(m, m);
    correlation_matrix.eye();
    for (int i = 0; i < m; ++i) {
        arma::vec cur_x = data.col(i);
        arma::uvec ord = sort_index(cur_x);
        arma::vec sorted_x(n);
        for (int k = 0; k < n; ++k) {
            sorted_x(k) = cur_x(ord(k));
        }
        for (int j = i + 1; j < m; ++j) {
            arma::vec cur_y(n);
            for (int k = 0; k < n; ++k) {
                cur_y(k) = data(ord(k), j);
            }
            correlation_matrix(i, j) = correlation_matrix(j, i) = kendall_tau(sorted_x, cur_y);
        }
    }
    return correlation_matrix;
}

arma::mat transform_correlation(const arma::mat correlation_matrix, Rcpp::String method)
{
    arma::mat tcor = arma::zeros(correlation_matrix.n_rows, correlation_matrix.n_cols);
    if (method == "spearman") {
        tcor = arma::sin(M_PI * correlation_matrix/6.0) * 2.0;
    } else if (method == "kendall") {
        tcor = arma::sin(M_PI * correlation_matrix/2.0);
    } else {
        Rcpp::stop("transform_correlation: method not recognized");
    }
    return tcor;
}

// [[Rcpp::export(.make_correlation)]]
arma::mat make_correlation(const arma::mat data, Rcpp::String method)
{
    arma::mat correl = arma::zeros(data.n_cols, data.n_cols);
    if (method == "pearson") {
        correl = arma::cor(data);
    } else if (method == "spearman") {
        correl = spearman_correlation(data);
    } else if (method == "kendall") {
        correl = kendall_correlation(data);
    } else {
        Rcpp::stop("make_correlation: method not recognized");
    }
    return correl;
}

// [[Rcpp::export(.make_psd)]]
arma::mat make_psd(arma::mat x, const double eig_tol, const double conv_tol, const double posd_tol, const int maxit) {
    arma::mat D_S = arma::zeros(x.n_rows, x.n_cols);
    arma::mat X = x;
    int iter = 0;
    bool converged = false;
    double conv = std::numeric_limits<double>::infinity();
    while (iter < maxit && !converged) {
        arma::mat Y = X;
        arma::mat R = Y - D_S;
        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, R);

        arma::uvec p = arma::find(eigval > eig_tol * eigval(0));
        if (p.n_elem == 0) {
            Rcpp::stop("Matrix appears negative semi-definite");
        }

        arma::mat Q_p = eigvec.cols(p);
        arma::vec d_p = eigval(p);

        X = Q_p * arma::diagmat(d_p) * Q_p.t();
        D_S = X - R;
        X = (X + X.t()) / 2;
        X.diag().ones();

        conv = arma::norm(Y - X, "inf") / arma::norm(Y, "inf");
        iter++;
        converged = (conv <= conv_tol);
    }

    if (!converged) {
        Rcpp::warning("'make_posdef' did not converge in " + std::to_string(iter) + " iterations");
    }
    return X;
}

// [[Rcpp::export(.p2P)]]
arma::mat p2P(arma::vec values, const int m)
{
    arma::mat P = arma::zeros(m,m);
    int index = 0;
    for(int j = 0;j < m; j++) {
        for(int i = j + 1; i < m; i++) {
            P(i,j) = values[index];
            P(j,i) = values[index];
            index++;
        }
    }
    P.diag().ones();
    return P;
}

// [[Rcpp::export(.is_psd)]]
bool is_psd(arma::mat x) {
    arma::mat L;
    bool chol_success = arma::chol(L, x, "lower");
    if (chol_success) {
        return true;
    } else {
        return false;
    }
}


// [[Rcpp::export(.array_mean)]]
arma::mat array_mean(arma::cube x)
{
    int s = x.n_slices;
    arma::mat z = x.slice(0);
    for(int i=1;i<s;++i) {
        z = z + x.slice(i);
    }
    z = z/s;
    return z;
}
