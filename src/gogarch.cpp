// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp14)]]
#include "gogarch.h"
using namespace Rcpp;

static inline arma::vec lower_triangular(const arma::mat M, const int diagonal) {
    arma::uvec lower_indices = arma::trimatl_ind(size(M), diagonal);
    arma::vec lower_tri = M(lower_indices);
    return lower_tri;
}

// [[Rcpp::export(.gogarch_covariance)]]
arma::mat gogarch_covariance(const arma::mat V, const arma::mat A) {
    int n = V.n_rows;
    int m = A.n_cols;
    int k = (m * m + m)/2;
    arma::mat covariance = arma::zeros(n, k);
    arma::mat tmp = arma::zeros(m,m);
    const int diagonal = 0;
    for (int i = 0;i<n;++i){
        tmp = A.t() * arma::diagmat(V.row(i)) * A;
        covariance.row(i) = lower_triangular(tmp, diagonal).t();
    }
    return covariance;
}

// [[Rcpp::export(.gogarch_correlation)]]
arma::mat gogarch_correlation(const arma::mat V, const arma::mat A) {
    int n = V.n_rows;
    int m = A.n_cols;
    int k = (m * m - m)/2;
    arma::mat correlation = arma::zeros(n, k);
    arma::mat tmp1 = arma::zeros(m,m);
    arma::mat tmp2 = arma::zeros(m,m);
    arma::mat tmp3 = arma::zeros(m,m);
    const int diagonal = -1;
    for (int i = 0;i<n;++i){
        tmp1 = A.t() * arma::diagmat(V.row(i)) * A;
        tmp2 = arma::diagmat(1/arma::sqrt(arma::diagvec(tmp1)));
        tmp3 = tmp2 * tmp1 * tmp2.t();
        correlation.row(i) =  lower_triangular(tmp3, diagonal).t();
    }
    return correlation;
}

// [[Rcpp::export(.coskewness_sigma)]]
arma::mat coskewness_sigma(const arma::vec sigmas) {
    int n = sigmas.n_elem;
    arma::umat idx1 = arma::sort(arma::vectorise(arma::repmat(arma::regspace<arma::uvec>(0, n - 1), 1, n * n)));
    arma::umat idx2 = arma::repmat(arma::sort(arma::vectorise(arma::repmat(arma::regspace<arma::uvec>(0, n - 1), 1, n))), n, 1);
    arma::umat idx3 = arma::repmat(arma::regspace<arma::uvec>(0, n - 1), n * n, 1);
    arma::umat idx = arma::join_horiz(idx1, idx2, idx3);

    arma::vec results(n * n * n);
    for (unsigned int i = 0; i < idx.n_rows; ++i) {
        results(i) = sigmas(idx(i, 0)) * sigmas(idx(i, 1)) * sigmas(idx(i, 2));
    }
    arma::mat cs = arma::reshape(results, n, n * n);
    return cs;
}



// [[Rcpp::export(.coskewness_block)]]
arma::mat coskewness_block(const arma::rowvec skew) {
    int n = skew.n_cols, i;
    int nn = n*n;
    arma::mat S(n, nn);
    S.zeros();
    for(i=0;i<n;i++) {
        int indx = (i)*nn + (i)*n + i;
        S(indx) = skew(i);
    }
    return S;
}


gogarch_coskewness_worker::gogarch_coskewness_worker(const arma::mat& S, const arma::mat& A,
                                             const arma::mat& V, const arma::mat& kronA,
                                             bool standardize, arma::cube& result)
    : S(S), A(A), V(V), kronA(kronA), standardize(standardize), result(result) {}

// Operator() implementation
void gogarch_coskewness_worker::operator()(std::size_t begin, std::size_t end) {
    arma::mat At = A.t();
    for (size_t i = begin; i < end; ++i) {
        arma::rowvec s_row = S.row(i);
        arma::mat SK = coskewness_block(s_row);
        arma::mat AV = At * SK;
        arma::mat r_slice = AV * kronA;
        if (standardize) {
            arma::vec s = arma::sqrt(arma::diagvec(At * arma::diagmat(V.row(i)) * A));
            arma::mat tmp = coskewness_sigma(s);
            r_slice /= tmp;
        }
        result.slice(i) = r_slice;
    }
}

// [[Rcpp::export(.gogarch_coskewness)]]
arma::cube gogarch_coskewness(const arma::mat A, const arma::mat S, const arma::mat V, bool standardize) {
    int M = A.n_cols;
    int N = A.n_rows;
    int T = S.n_rows;
    arma::cube result(M, M * M, T);
    arma::mat kronA = arma::kron(A, A);
    arma::mat AV = arma::zeros(M, N);
    arma::mat SK = arma::zeros(N, N * N);
    arma::mat At = A.t();
    gogarch_coskewness_worker worker(S, A, V, kronA, standardize, result);
    RcppParallel::parallelFor(0, T, worker);
    return result;

}

// [[Rcpp::export(.combn)]]
arma::umat combn(const arma::uvec n, int m) {
    int len = n.size();
    if (m > len) Rcpp::stop("m cannot be greater than the length of n");
    // Compute the number of combinations using lgamma to prevent overflow
    double log_combinations = std::lgamma(len + 1) - std::lgamma(m + 1) - std::lgamma(len - m + 1);
    int combinations = static_cast<int>(std::round(std::exp(log_combinations)));
    arma::umat results(m, combinations);
    arma::uvec indices(m);
    for (int i = 0; i < m; i++) {
        indices[i] = i;
    }
    int col = 0;
    do {
        for (int i = 0; i < m; i++) {
            results(i, col) = n(indices[i]);
        }
        col++;
        int t = m - 1;
        while (t != -1 && indices[t] == len - m + t) t--;
        if (t == -1) break;
        indices[t]++;
        for (int i = t + 1; i < m; i++) {
            indices[i] = indices[t] + i - t;
        }
    } while (true);

    return results;
}

// [[Rcpp::export(.cokurt_pairs)]]
arma::field<arma::umat> cokurtosis_pairs(int n) {
  // Safely calculate the number of unique pairs using combination formula
  int unique_pairs = (n * (n - 1)) / 2;

  // Initialize the idx matrix to store pairs (6 columns)
  arma::umat idx(unique_pairs, 6);

  // Create indices and generate all combinations of 2 elements
  arma::uvec indices = arma::linspace<arma::uvec>(1, n, n);
  arma::umat prsx = combn(indices, 2);  // Ensure this function works correctly
  arma::umat prs(prsx.n_cols, 3);

  // Fill prs with combinations and differences
  prs.col(0) = prsx.row(1).t();
  prs.col(1) = prsx.row(0).t();
  prs.col(2) = prsx.row(1).t() - prsx.row(0).t();

  // Sort the rows based on the difference between columns 0 and 1
  prs = prs.rows(arma::sort_index(prs.col(0) - prs.col(1)));

  // Initialize the first row of idx safely
  unsigned int ix = n + 2;  // First index initialization
  arma::uvec first = {ix,
                      static_cast<unsigned int>((n-1) * n),
                      static_cast<unsigned int>(n-1),
                      static_cast<unsigned int>((n+1) * (n-1) * (n-1)),
                      static_cast<unsigned int>(n-1),
                      static_cast<unsigned int>((n-1) * n)};

  // Ensure the first indices are within bounds
  if (arma::all(first < std::numeric_limits<unsigned int>::max())) {
    idx.row(0) = arma::cumsum(first).t();
  } else {
    Rcpp::Rcout << "Error: First indices out of bounds" << std::endl;
  }

  // Fill in the rest of the idx matrix for each unique pair
  for (int i = 1; i < unique_pairs; ++i) {
    unsigned int d = prs(i, 2);
    unsigned int i1;

    // Ensure safe index calculation and prevent overflow
    if (d == prs(i - 1, 2)) {
      i1 = idx(i - 1, 0) + (n * n * n) + (n * n) + n + 1;  // Avoid std::pow for large numbers
    } else {
      i1 = ix + n + 1;
      ix = i1;  // Update ix
    }

    // Construct the current indices for this pair
    arma::uvec current = {i1,
                          d * (n-1) * n,
                          d * (n-1),
                          d * (n+1) * (n-1) * (n-1),
                          d * (n-1),
                          d * (n-1) * n};

    // Ensure the current indices are within bounds before storing
    if (arma::all(current < std::numeric_limits<unsigned int>::max())) {
      idx.row(i) = arma::cumsum(current).t();
    } else {
      Rcpp::Rcout << "Error: Index out of bounds for row " << i << std::endl;
    }
  }

  // Store idx and prs in a field to return
  arma::field<arma::umat> result(2);
  result(0) = idx;
  result(1) = prs;

  return result;
}
// arma::field<arma::umat> cokurtosis_pairs(int n) {
//     Rcpp::NumericVector vec = Rcpp::NumericVector::create(n, n-2, 2);
//     Rcpp::NumericVector fact_vec = Rcpp::factorial(vec);
//     int unique_pairs = fact_vec[0] / (fact_vec[1] * fact_vec[2]);
//     arma::umat idx(unique_pairs, 6);
//     arma::uvec indices = arma::linspace<arma::uvec>(1, n, n);
//     arma::umat prsx = combn(indices, 2);
//     arma::umat prs(prsx.n_cols, 3);
//     prs.col(0) = prsx.row(1).t();
//     prs.col(1) = prsx.row(0).t();
//     prs.col(2) = prsx.row(1).t() - prsx.row(0).t();
//     prs = prs.rows(arma::sort_index(prs.col(0) - prs.col(1)));
//
//     unsigned int ix = n + 2;
//     arma::uvec first = {ix, static_cast<unsigned int>((n-1) * n), static_cast<unsigned int>(n-1), static_cast<unsigned int>((n+1) * (n-1) * (n-1)), static_cast<unsigned int>(n-1), static_cast<unsigned int>((n-1) * n)};
//     idx.row(0) = arma::cumsum(first).t();
//
//     if (unique_pairs > 1) {
//         for (int i = 1; i < unique_pairs; ++i) {
//             unsigned int d = prs(i, 2);
//             unsigned int i1;
//             if (d == prs(i-1, 2)) {
//                 i1 = idx(i-1, 0) + std::pow(n, 3) + std::pow(n, 2) + n + 1;
//             } else {
//                 i1 = ix + n + 1;
//                 ix = i1;
//             }
//             arma::uvec current = {i1, d * (n-1) * n, d * (n-1), d * (n+1) * (n-1) * (n-1), d * (n-1), d * (n-1) * n};
//             idx.row(i) = arma::cumsum(current).t();
//         }
//     }
//
//     arma::field<arma::umat> result(2);
//     result(0) = idx;
//     result(1) = prs;
//
//     return result;
// }


// [[Rcpp::export(.cokurt_index)]]
arma::mat cokurtosis_block(const arma::vec s, const arma::vec values)
{
    int n = s.n_elem;
    arma::field<arma::umat> x = cokurtosis_pairs(n);
    arma::uvec indices = arma::regspace<arma::uvec>(1, n);
    int n2 = n * n;
    int n3 = n * n * n;
    int n4 = n3 * n;
    arma::uvec ix = (indices - 1) * n3 + (indices - 1) * n2 + (indices - 1) * n + indices - 1;
    arma::vec z(n4, arma::fill::zeros);
    z(ix) = values;
    arma::umat indx = x(0);
    arma::umat pairs = x(1);
    int k = indx.n_rows;
    for(int i = 0;i < k; ++i) {
        arma::uvec row = indx.row(i).t();
        double value = s(pairs(i,0) - 1) * s(pairs(i,1) - 1);
        z.elem(row - 1) = arma::vec(row.n_elem, arma::fill::value(value));
    }
    arma::mat zmat = arma::reshape(z, n, n * n * n);
    return zmat;
}

// [[Rcpp::export(.cokurtosis_sigma)]]
arma::mat cokurtosis_sigma(const arma::vec sigmas) {
    int n = sigmas.n_elem;
    arma::umat idx1 = arma::sort(arma::vectorise(arma::repmat(arma::regspace<arma::uvec>(0, n - 1), 1, n * n * n)));
    arma::umat idx2 = arma::repmat(arma::sort(arma::vectorise(arma::repmat(arma::regspace<arma::uvec>(0, n - 1), 1, n * n))), n, 1);
    arma::umat idx3 = arma::repmat(arma::sort(arma::vectorise(arma::repmat(arma::regspace<arma::uvec>(0, n - 1), 1, n))), n * n, 1);
    arma::umat idx4 = arma::repmat(arma::regspace<arma::uvec>(0, n - 1), n * n * n, 1);
    arma::umat idx = arma::join_horiz(idx1, idx2, idx3, idx4);

    arma::vec results(n * n * n * n);
    for (unsigned int i = 0; i < idx.n_rows; ++i) {
        results(i) = sigmas(idx(i, 0)) * sigmas(idx(i, 1)) * sigmas(idx(i, 2)) * sigmas(idx(i, 3));
    }
    arma::mat ks = arma::reshape(results, n, n * n * n);
    return ks;
}


gogarch_cokurtosis_worker::gogarch_cokurtosis_worker(const arma::mat& K, const arma::mat& A, const arma::mat& V,
                                                     const arma::mat& kronA, bool standardize, arma::cube& result)
    : K(K), A(A), V(V), kronA(kronA), standardize(standardize), result(result) {}

void gogarch_cokurtosis_worker::operator()(std::size_t begin, std::size_t end) {
    arma::mat At = A.t();
    arma::mat V_t = V.t();
    arma::mat K_t = K.t();
    for (size_t i = begin; i < end; ++i) {
        arma::vec v_col = V_t.col(i);
        arma::vec k_col = K_t.col(i);
        arma::mat KU = cokurtosis_block(v_col, k_col);
        arma::mat AV = At * KU;
        arma::mat r_slice = AV * kronA;
        if (standardize) {
            arma::vec s = arma::sqrt(arma::diagvec(At * arma::diagmat(V.row(i)) * A));
            arma::mat tmp = cokurtosis_sigma(s);
            r_slice /= tmp;
        }
        result.slice(i) = r_slice;
    }
}


// [[Rcpp::export(.gogarch_cokurtosis)]]
arma::cube gogarch_cokurtosis(const arma::mat A, const arma::mat K, const arma::mat V, bool standardize) {
    int M = A.n_cols;
    int N = A.n_rows;
    int T = K.n_rows;
    arma::cube result(M, M * M * M, T);
    arma::mat kronA = arma::kron(A, A);
    kronA = arma::kron(kronA, A);
    arma::mat AV = arma::zeros(M, N);
    arma::mat KU = arma::zeros(N, N * N * N);
    gogarch_cokurtosis_worker worker(K, A, V, kronA, standardize, result);
    RcppParallel::parallelFor(0, T, worker);
    return result;
}


gogarch_coskewness_weighted_worker::gogarch_coskewness_weighted_worker(
    const arma::mat& S, const arma::mat& A, const arma::mat& kronA,
    const arma::mat& W, arma::vec& result)
    : S(S), A(A), kronA(kronA), W(W), result(result) {}

void gogarch_coskewness_weighted_worker::operator()(std::size_t begin, std::size_t end) {
    arma::mat At = A.t();

    for (size_t i = begin; i < end; ++i) {
        arma::rowvec s_row = S.row(i);
        arma::mat SK = coskewness_block(s_row);
        arma::mat AV = At * SK;
        arma::mat r_slice = AV * kronA;
        arma::rowvec w_i = W.row(i);
        arma::mat w_kron_i = arma::kron(w_i.t(), w_i.t());
        result(i) = arma::as_scalar(w_i * r_slice * w_kron_i);
    }
}

// [[Rcpp::export(.gogarch_skewness_weighted)]]
arma::vec gogarch_skewness_weighted(const arma::mat A, const arma::mat S,
                                    const arma::mat w) {
    int T = S.n_rows;
    arma::vec result(T);
    arma::mat kronA = arma::kron(A, A);
    gogarch_coskewness_weighted_worker worker(S, A, kronA, w, result);
    RcppParallel::parallelFor(0, T, worker);
    return result;
}


gogarch_cokurtosis_weighted_worker::gogarch_cokurtosis_weighted_worker(
    const arma::mat& K, const arma::mat& V, const arma::mat& A, const arma::mat& kronA,
    const arma::mat& W, arma::vec& result)
    : K(K), V(V), A(A), kronA(kronA), W(W), result(result) {}

void gogarch_cokurtosis_weighted_worker::operator()(std::size_t begin, std::size_t end) {
    arma::mat At = A.t();
    arma::mat V_t = V.t();
    arma::mat K_t = K.t();
    for (size_t i = begin; i < end; ++i) {
        arma::vec v_col = V_t.col(i);
        arma::vec k_col = K_t.col(i);
        arma::mat KU = cokurtosis_block(v_col, k_col);
        arma::mat AV = At * KU;
        arma::mat r_slice = AV * kronA;
        arma::rowvec w_i = W.row(i);
        arma::vec w_it = w_i.t();
        arma::mat tmpw = arma::kron(w_it, w_it);
        arma::mat w_kron_i = arma::kron(w_it, tmpw);
        arma::mat w_kron_j = arma::kron(w_it, w_kron_i);
        arma::vec rvec = arma::vectorise(r_slice);
        result(i) = arma::as_scalar(rvec.t() * w_kron_j);
    }
}

// [[Rcpp::export(.gogarch_kurtosis_weighted)]]
arma::vec gogarch_cokurtosis_weighted(const arma::mat A, const arma::mat K,
                                      const arma::mat V, const arma::mat w) {
    int T = K.n_rows;
    arma::vec result(T);
    arma::mat tmpA = arma::kron(A, A);
    arma::mat kronA = arma::kron(tmpA, A);
    gogarch_cokurtosis_weighted_worker worker(K, V, A, kronA, w, result);
    RcppParallel::parallelFor(0, T, worker);
    return result;
}

// [[Rcpp::export(.gogarch_covariance_weighted)]]
arma::vec gogarch_covariance_weighted(const arma::mat V, const arma::mat A, const arma::mat w) {
    int n = V.n_rows;
    int m = A.n_cols;
    arma::vec wcovariance = arma::zeros(n);
    arma::mat tmp = arma::zeros(m,m);
    for (int i = 0;i<n;++i){
        tmp = A.t() * arma::diagmat(V.row(i)) * A;
        wcovariance(i) = arma::as_scalar(w.row(i) * tmp * w.row(i).t());
    }
    return wcovariance;
}


gogarch_cokurtosis_weighted_worker_sim::gogarch_cokurtosis_weighted_worker_sim(
  const arma::cube& sig, const arma::mat& ku, const arma::mat& A,
  const arma::mat& kronA, const arma::mat& W, arma::mat& result,
  int n, int nsim)
  : sig(sig), ku(ku), A(A), kronA(kronA), W(W), result(result), n(n), nsim(nsim) {}

void gogarch_cokurtosis_weighted_worker_sim::operator()(std::size_t begin, std::size_t end) {
  arma::mat At = A.t();
  for (std::size_t idx = begin; idx < end; ++idx) {
    int sim_idx = idx / n;
    int col_idx = idx % n;
    arma::mat sig_mat = sig.slice(sim_idx);
    arma::mat V = arma::pow(sig_mat, 2);
    arma::mat V_t = V.t();
    arma::mat K = ku % arma::pow(sig_mat, 4);
    arma::mat K_t = K.t();
    arma::vec v_col = V_t.col(col_idx);
    arma::vec k_col = K_t.col(col_idx);
    arma::mat KU = cokurtosis_block(v_col, k_col);
    arma::mat AV = At * KU;
    arma::mat r_slice = AV * kronA;
    arma::rowvec w_i = W.row(col_idx);
    arma::vec w_it = w_i.t();
    arma::mat tmpw = arma::kron(w_it, w_it);
    arma::mat w_kron_i = arma::kron(w_it, tmpw);
    arma::mat w_kron_j = arma::kron(w_it, w_kron_i);
    arma::vec rvec = arma::vectorise(r_slice);
    result(sim_idx, col_idx) = arma::as_scalar(rvec.t() * w_kron_j);;
  }
}

// [[Rcpp::export(.gogarch_kurtosis_weighted_sim)]]
arma::mat gogarch_cokurtosis_weighted_sim(const arma::mat& A, const arma::cube& sig,
                                           const arma::mat& ku, const arma::mat& weights, int nsim, int n) {
  arma::mat tmpA = arma::kron(A, A);
  arma::mat kronA = arma::kron(tmpA, A);
  arma::mat result(nsim, n);
  gogarch_cokurtosis_weighted_worker_sim worker(sig, ku, A, kronA, weights, result, n, nsim);
  RcppParallel::parallelFor(0, n * nsim, worker);
  return result;
}
