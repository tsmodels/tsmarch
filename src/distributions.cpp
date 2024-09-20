// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppBessel)]]
// [[Rcpp::plugins(cpp14)]]
#include "distributions.h"
using namespace Rcpp;
using namespace std::complex_literals;

// [[Rcpp::export(.rmvnorm)]]
arma::mat rmvnorm(arma::mat R, arma::mat Z){
    Rcpp::RNGScope scope;
    int m = R.n_rows;
    int n = Z.n_rows;
    arma::vec eigval(m);
    arma::mat eigvec(m, m);
    arma::mat temp(m, m);
    arma::eig_sym(eigval, eigvec, R);
    arma::mat ans = arma::zeros(n, m);
    temp = (eigvec * arma::diagmat(arma::sqrt(eigval)) * eigvec.t());
    for (int i = 0;i<n;++i){
        ans.row(i) = Z.row(i) * temp;
    }
    return ans;
}

// [[Rcpp::export(.rmvt)]]
arma::mat rmvt(arma::mat R, arma::mat Z, const double nu){
    Rcpp::RNGScope scope;
    int m = R.n_rows;
    int n = Z.n_rows;
    arma::mat R_adj = ((nu - 2.0)/nu) * R;
    arma::vec eigval(m);
    arma::mat eigvec(m, m);
    arma::mat temp(m, m);
    arma::eig_sym(eigval, eigvec, R_adj);
    temp = (eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::inv(eigvec));
    arma::mat ans = arma::zeros(n,m);
    for (int i = 0;i<n;++i) {
        double rc = Rf_rchisq(nu);
        double v = sqrt(nu/rc);
        ans.row(i) = v * (Z.row(i) * temp);
    }
    return ans;
}

arma::rowvec rmvt(arma::mat R, arma::rowvec Z, const double nu, const double rc){
    Rcpp::RNGScope scope;
    int m = R.n_rows;
    arma::mat R_adj = ((nu - 2.0)/nu) * R;
    arma::vec eigval(m);
    arma::mat eigvec(m, m);
    arma::mat temp(m, m);
    arma::eig_sym(eigval, eigvec, R_adj);
    temp = (eigvec * arma::diagmat(arma::sqrt(eigval)) * arma::inv(eigvec));
    double v = sqrt(nu/rc);
    arma::rowvec ans = v * (Z * temp);
    return ans;
}

arma::rowvec rmvnorm(arma::mat R, arma::rowvec Z){
    Rcpp::RNGScope scope;
    int m = R.n_rows;
    arma::vec eigval(m);
    arma::mat eigvec(m, m);
    arma::mat temp(m, m);
    arma::eig_sym(eigval, eigvec, R);
    temp = (eigvec * arma::diagmat(arma::sqrt(eigval)) * eigvec.t());
    arma::rowvec ans = Z * temp;
    return ans;
}

arma::mat mpnorm(const arma::mat& x) {
    Rcpp::NumericVector vec(x.begin(), x.end());
    Rcpp::NumericVector pnorm_vec = Rcpp::pnorm(vec, 0.0, 1.0);
    arma::mat result(pnorm_vec.begin(), x.n_rows, x.n_cols);
    return result;
}

arma::mat mpstd(const arma::mat& x, const double shape) {
    double s = sqrt(shape / (shape - 2.0));
    Rcpp::NumericVector vec(x.begin(), x.end());
    vec = vec * s;
    Rcpp::NumericVector pt_vec = Rcpp::pt(vec, shape);
    arma::mat result(pt_vec.begin(), x.n_rows, x.n_cols);
    return result;
}

arma::uvec interval(const arma::vec& x, const arma::vec& y) {
    arma::uvec i1 = arma::find(x < y.min());
    arma::uvec i2 = arma::find(x > y.max());

    arma::uvec i = arma::regspace<arma::uvec>(i1.min(), i2.min());

    int j = i(0);
    while (x(j) > y(0)) {
        j--;
    }

    int k = i.max();
    while (x(k) < y.max()) {
        k++;
    }

    arma::uvec result;
    if (i(0) - 1 != 0) {
        result = arma::regspace<arma::uvec>(j, i(0) - 1);
    }
    result = arma::join_vert(result, i);
    if (i.max() + 1 < k) {
        result = arma::join_vert(result, arma::regspace<arma::uvec>(i.max() + 1, k));
    }

    result = arma::unique(result);
    result = arma::sort(result);

    return result;
}

// [[Rcpp::export(.interpolate_window)]]
arma::vec interpolate_window(const arma::vec& x, const arma::vec& y, const arma::vec& z, int w = -1) {
    if (w == -1) {
        w = z.n_elem;
    }
    int m = z.n_elem;
    int n = std::floor(m / w);
    arma::vec r;

    for (int i = 0; i < n; i++) {
        arma::vec dz = z.subvec(i * w, (i + 1) * w - 1);
        arma::uvec k = interval(x, dz);
        arma::vec yz(dz.n_elem);
        arma::interp1(x.elem(k), y.elem(k), dz, yz);
        r = arma::join_vert(r, yz);
    }
    if (m % w != 0) {
        int e = m % w;
        arma::vec dz = z.subvec(n * w, n * w + e - 1);
        arma::uvec k = interval(x, dz);
        arma::vec yz(dz.n_elem);
        arma::interp1(x.elem(k), y.elem(k), dz, yz);
        r = arma::join_vert(r, yz);
    }

    return r;
}

// [[Rcpp::export(.nigmvcf)]]
arma::cx_vec nigmvcf(const arma::vec& z, const arma::vec& alpha, const arma::vec& beta, const arma::vec& delta, const arma::vec& mu) {
    int N = z.n_elem;
    int m = mu.n_elem;

    arma::cx_vec x1 = 1i * z * sum(mu);
    arma::cx_mat zx(N, m);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < m; ++j) {
            zx(i, j) = delta[j] * (std::sqrt(alpha[j] * alpha[j] - beta[j] * beta[j]) -
                std::sqrt(alpha[j] * alpha[j] - (beta[j] + 1i * z[i]) * (beta[j] + 1i * z[i])));
        }
    }

    arma::cx_vec x2 = sum(zx, 1);

    return exp(x1 + x2);
}

// [[Rcpp::export(.cfinvnig)]]
arma::vec cfinvnig(const arma::vec& z, double step, const arma::vec& alpha, const arma::vec& beta, const arma::vec& delta, const arma::vec& mu) {
    int pmax = 18;
    int p = 14;
    double maxz = round(arma::max(arma::abs(z))) + 5;
    while ((maxz / step + 1) > std::pow(2, (p - 1))) {
        p++;
    }
    if (p > pmax) p = pmax;
    if ((maxz / step + 1) > std::pow(2, (p - 1))) {
        step = (maxz + 1) * (1 + step / 10) / std::pow(2, (p - 1));
    }
    arma::vec zs = arma::sort(z);
    int n = std::pow(2, p);
    arma::vec x = arma::linspace(0, n - 1, n) * step - (n * step / 2);
    double s = 1 / (step * n);
    arma::vec tt = 2 * M_PI * s * (arma::linspace(0, n - 1, n) - n / 2);
    arma::vec sgn = arma::ones<arma::vec>(n);
    for (int i = 1; i < n; i += 2) {
        sgn[i] = -1;
    }
    arma::cx_vec cf = nigmvcf(tt, alpha, beta, delta, mu);
    arma::cx_vec phi = sgn % cf;
    phi[n / 2] = sgn[n / 2];
    arma::vec p_result = s * abs(arma::fft(phi));
    arma::vec pdf = interpolate_window(x, p_result, zs, -1);
    return pdf;
}

// [[Rcpp::export(.ghypfn)]]
std::complex<double> ghypfn(double lambda, double alpha, double beta, double delta, double z) {
    std::complex<double> i(0.0, 1.0);
    std::complex<double> beta_plus_iz = beta + i * z;
    std::complex<double> x1 = (lambda / 2.0) *
        (std::log(alpha * alpha - beta * beta) -
        std::log(alpha * alpha - std::pow(beta_plus_iz, 2)));
    std::complex<double> bessel_arg1 = delta * std::sqrt(alpha * alpha - std::pow(beta_plus_iz, 2));
    SEXP bessel_k1_sexp = RcppBessel::bessel_k(Rcpp::wrap(bessel_arg1), lambda);
    std::complex<double> bessel_k1 = Rcpp::as<std::complex<double>>(bessel_k1_sexp);
    double bessel_arg2 = delta * std::sqrt(alpha * alpha - beta * beta);
    SEXP bessel_k2_sexp = RcppBessel::bessel_k(Rcpp::wrap(bessel_arg2), lambda);
    std::complex<double> bessel_k2 = Rcpp::as<std::complex<double>>(bessel_k2_sexp);
    std::complex<double> x2 = std::log(bessel_k1) - std::log(bessel_k2);
    return x1 + x2;
}


// [[Rcpp::export(.ghypmvcf)]]
arma::cx_vec ghypmvcf(const arma::vec& z, const arma::vec& lambda, const arma::vec& alpha,
                      const arma::vec& beta, const arma::vec& delta, const arma::vec& mu) {
    int N = z.n_elem;
    int m = mu.n_elem;
    if(lambda.n_elem != m || alpha.n_elem != m || beta.n_elem != m || delta.n_elem != m) {
        stop("Vectors lambda, alpha, beta, delta, and mu must have the same length.");
    }
    arma::cx_double i(0.0, 1.0);
    double sum_mu = arma::accu(mu);
    arma::cx_vec result = i * z * sum_mu;

    arma::cx_mat zx(N, m, arma::fill::zeros);
    for (int j = 0; j < m; ++j) {
        for (int k = 0; k < N; ++k) {
            zx(k, j) = ghypfn(lambda[j], alpha[j], beta[j], delta[j], z[k]);
        }
    }
    arma::cx_vec x2 = sum(zx, 1);
    arma::cx_vec ans(N);
    for (int k = 0; k < N; ++k) {
        ans[k] = arma::as_scalar(std::exp(result[k] + x2[k]));
    }
    return ans;
}

// [[Rcpp::export(.cfinvghyp)]]
arma::vec cfinvghyp(const arma::vec& z, double step, const arma::vec& lambda, const arma::vec& alpha, const arma::vec& beta, const arma::vec& delta, const arma::vec& mu) {
    int pmax = 18;
    int p = 14;
    double maxz = round(arma::max(arma::abs(z))) + 5;
    while ((maxz / step + 1) > std::pow(2, (p - 1))) {
        p++;
    }
    if (p > pmax) p = pmax;
    if ((maxz / step + 1) > std::pow(2, (p - 1))) {
        step = (maxz + 1) * (1 + step / 10) / std::pow(2, (p - 1));
    }
    arma::vec zs = arma::sort(z);
    int n = std::pow(2, p);
    arma::vec x = arma::linspace(0, n - 1, n) * step - (n * step / 2);
    double s = 1 / (step * n);
    arma::vec tt = 2 * M_PI * s * (arma::linspace(0, n - 1, n) - n / 2);
    arma::vec sgn = arma::ones<arma::vec>(n);
    for (int i = 1; i < n; i += 2) {
        sgn[i] = -1;
    }
    arma::cx_vec cf = ghypmvcf(tt, lambda, alpha, beta, delta, mu);
    arma::cx_vec phi = sgn % cf;
    phi[n / 2] = sgn[n / 2];
    arma::vec p_result = s * abs(arma::fft(phi));
    arma::vec pdf = interpolate_window(x, p_result, zs, -1);
    return pdf;
}
