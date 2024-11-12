// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]
#include "dcc.h"
using namespace Rcpp;

// [[Rcpp::export(.dcc_constant_normal)]]
Rcpp::List dcc_constant_normal(const arma::mat Z, const arma::mat S)
{
    int timesteps = Z.n_rows;
    int m = Z.n_cols;
    // not used
    Rcpp::String distribution = "mvn";
    Rcpp::String method = "pearson";
    arma::mat R = make_correlation(Z, method);
    bool ispd = is_psd(R);
    if (!ispd) {
        R = make_psd(R);
    }
    arma::mat r_inverse = arma::inv_sympd(R);
    double part1 = arma::log_det_sympd(R);
    double const_part = m * log(2 * M_PI);
    double part2 = 0.0;
    double garchll = 0.0;
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);

    for (int i = 0;i<timesteps;i++) {
        part2 = arma::as_scalar(Z.row(i) * Z.row(i).t());
        garchll = 0.5 * (const_part + 2.0 * log(arma::prod(S.row(i))) + part2);
        dccll_vec(i) =  0.5 * (arma::as_scalar(Z.row(i) * (r_inverse  * Z.row(i).t())) + part1 - part2);
        ll_vec(i) = dccll_vec(i) + garchll;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("R") = R , _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec,  _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}

// [[Rcpp::export(.dcc_constant_student)]]
Rcpp::List dcc_constant_student(const arma::mat Z, const arma::mat S, const double shape)
{
    int timesteps = Z.n_rows;
    int m = Z.n_cols;
    // not used
    Rcpp::String distribution = "mvt";
    Rcpp::String method = "kendall";
    arma::mat R = make_correlation(Z, method);
    R = transform_correlation(R, method);
    bool ispd = is_psd(R);
    if (!ispd) {
        R = make_psd(R);
    }
    double const_term = -1.0 * lgamma(0.5 * (m + shape)) + lgamma(0.5 * shape) + 0.5 * m * log(M_PI * (shape - 2.0));
    arma::mat r_inverse = arma::inv_sympd(R);
    double part1 = arma::log_det_sympd(R);
    double part2 = 0.0;
    double part3 = 0.0;
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);

    for (int i = 0;i<timesteps;i++) {
        part2 = arma::as_scalar(Z.row(i) * (r_inverse  * Z.row(i).t()));
        part3 = log(arma::prod(S.row(i)));
        dccll_vec(i) = const_term + 0.5 * part1 + 0.5 * (shape + m) * log(1.0 + (1.0/(shape - 2.0)) * part2);
        ll_vec(i) = dccll_vec(i) + part3;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("R") = R ,  _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec,  _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}

// [[Rcpp::export(.dcc_dynamic_normal)]]
Rcpp::List dcc_dynamic_normal(const arma::vec alpha, const arma::vec gamma, const arma::vec beta,
                                 const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder)
{
    const int alpha_order = dccorder[0];
    const int gamma_order = dccorder[1];
    const int beta_order = dccorder[2];
    const int maxpq = std::max(alpha_order, beta_order);
    const int m = z.n_cols;
    double sum_alpha = arma::accu(alpha);
    double sum_gamma = arma::accu(gamma);
    double sum_beta = arma::accu(beta);
    const arma::mat zero_matrix = arma::zeros(maxpq, m);
    Rcpp::String distribution = "mvn";
    arma::mat Qbar = arma::cov(z);
    arma::mat AsyZinit = arma::zeros(z.n_rows, m);
    arma::mat Nbar = arma::zeros(m, m);
    if (gamma_order > 0) {
        AsyZinit = matrix_sign(z) % z;
        Nbar = arma::cov(AsyZinit);
    }
    const int timesteps = z.n_rows + maxpq;
    const arma::mat Z = arma::join_cols(zero_matrix, z);
    const arma::mat S = arma::join_cols(zero_matrix, s);
    arma::mat AsyZ = arma::join_cols(zero_matrix, AsyZinit);
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);
    arma::mat Omega = Qbar * (1.0 - sum_alpha - sum_beta);
    if (gamma_order > 0) {
        Omega = Omega - sum_gamma * Nbar;
    }
    arma::cube C(m,m,timesteps);
    arma::cube R(m,m,timesteps);
    int i,j;
    for(i = 0;i<maxpq;i++){
        C.slice(i) = Qbar;
    }
    arma::mat Q = arma::zeros(m,m);

    double const_part = m * log(2 * M_PI);
    double part1 = 0.0;
    double part2 = 0.0;
    double garch_ll = 0.0;

    for (i = maxpq; i < timesteps; ++i) {
        Q = Omega;
        if (alpha_order > 0) {
            for(j = 0;j<alpha_order;++j) {
                Q = Q + alpha(j) * (Z.row(i - (j + 1)).t() * Z.row(i - (j + 1)));
            }
        }
        if (gamma_order > 0) {
            for(j = 0;j<gamma_order;++j) {
                Q = Q + gamma(j) * (AsyZ.row(i - (j + 1)).t() * AsyZ.row(i - (j + 1)));
            }
        }
        if (beta_order > 0) {
            for(j = 0;j<beta_order;++j) {
                Q = Q + beta(j) * C.slice(i - (j + 1));
            }
        }
        C.slice(i) = Q;
        arma::vec tempx = arma::sqrt(Q.diag());
        arma::mat tempy = tempx * tempx.t();
        arma::mat Rt = Q/tempy;
        R.slice(i) = Rt;
        arma::mat id_matrix = arma::eye(m, m);
        arma::mat Rinv = arma::inv(Rt);
        part1 = arma::log_det_sympd(Rt);
        part2 = arma::as_scalar(Z.row(i) * Z.row(i).t());
        garch_ll = 0.5 * (const_part + 2.0 * log(arma::prod(S.row(i))) + part2);
        dccll_vec(i) =  0.5 * (arma::as_scalar(Z.row(i) * (Rinv  * Z.row(i).t())) + part1 - part2);
        ll_vec(i) = dccll_vec(i) + garch_ll;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("Qbar") = Qbar, _("Nbar") = Nbar, _("R") = R , _("Q") = C,  _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec, _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}


// [[Rcpp::export(.dcc_dynamic_student)]]
Rcpp::List dcc_dynamic_student(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, double shape, const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder)
{
    const int alpha_order = dccorder[0];
    const int gamma_order = dccorder[1];
    const int beta_order = dccorder[2];
    const int maxpq = std::max(alpha_order, beta_order);
    const int m = z.n_cols;
    double sum_alpha = arma::accu(alpha);
    double sum_gamma = arma::accu(gamma);
    double sum_beta = arma::accu(beta);
    arma::mat zero_matrix = arma::zeros(maxpq, m);
    Rcpp::String distribution = "mvt";
    arma::mat Qbar = arma::cov(z);
    arma::mat AsyZinit = arma::zeros(z.n_rows, m);
    arma::mat Nbar = arma::zeros(m, m);
    if (gamma_order > 0) {
        AsyZinit = matrix_sign(z) % z;
        Nbar = arma::cov(AsyZinit);
    }
    const int timesteps = z.n_rows + maxpq;
    arma::mat Z = arma::join_cols(zero_matrix, z);
    arma::mat S = arma::join_cols(zero_matrix, s);
    arma::mat AsyZ = arma::join_cols(zero_matrix, AsyZinit);
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);
    arma::mat Omega = Qbar * (1.0 - sum_alpha - sum_beta);
    if (gamma_order > 0) {
        Omega = Omega - sum_gamma * Nbar;
    }
    arma::cube C(m,m,timesteps);
    arma::cube R(m,m,timesteps);
    double const_term = -1.0 * lgamma(0.5 * (m + shape)) + lgamma(0.5 * shape) + 0.5 * m * log(M_PI * (shape - 2.0));
    double part1 = 0.0;
    double part2 = 0.0;
    double part3 = 0.0;
    int i,j;
    for(i = 0;i<maxpq;i++){
        C.slice(i) = Qbar;
    }
    arma::mat Q = arma::zeros(m,m);
    for (i = maxpq; i < timesteps; ++i) {
        Q = Omega;
        if (alpha_order > 0) {
            for(j = 0;j<alpha_order;++j) {
                Q = Q + alpha(j) * (Z.row(i - (j + 1)).t() * Z.row(i - (j + 1)));
            }
        }
        if (gamma_order > 0) {
            for(j = 0;j<gamma_order;++j) {
                Q = Q + gamma(j) * (AsyZ.row(i - (j + 1)).t() * AsyZ.row(i - (j + 1)));
            }
        }
        if (beta_order > 0) {
            for(j = 0;j<beta_order;++j) {
                Q = Q + beta(j) * C.slice(i - (j + 1));
            }
        }
        C.slice(i) = Q;
        arma::vec tempx = arma::sqrt(Q.diag());
        arma::mat tempy = tempx * tempx.t();
        arma::mat Rt = Q/tempy;
        R.slice(i) = Rt;
        arma::mat Rinv = arma::inv(Rt);
        part1 = arma::log_det_sympd(Rt);
        part2 = arma::as_scalar(Z.row(i) * (Rinv  * Z.row(i).t()));
        part3 = log(arma::prod(S.row(i)));
        dccll_vec(i) = const_term + 0.5 * part1 + 0.5 * (shape + m) * log(1.0 + (1.0/(shape - 2.0)) * part2);
        ll_vec(i) = dccll_vec(i) + part3;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("Qbar") = Qbar, _("Nbar") = Nbar, _("R") = R , _("Q") = C,  _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec,  _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}

// [[Rcpp::export(.adcc_constraint)]]
double adcc_constraint(const arma::vec alpha, const arma::vec gamma, const arma::vec beta, double shape, const arma::mat z, Rcpp::IntegerVector dccorder)
{
    double sum_alpha = arma::accu(alpha);
    double sum_gamma = arma::accu(gamma);
    double sum_beta = arma::accu(beta);
    arma::mat Qbar = arma::cov(z);
    arma::mat AsyZ = matrix_sign(z) % z;
    arma::mat Nbar = arma::cov(AsyZ);
    arma::mat Qbar2 = arma::inv_sympd(arma::sqrtmat_sympd(Qbar));
    arma::mat tmp = Qbar2 * Nbar * Qbar2;
    arma::vec eigenv = arma::eig_sym(tmp);
    double delta = eigenv.max();
    double out = sum_alpha + sum_beta + delta * sum_gamma;
    return out;
}


// [[Rcpp::export(.dcc_constant_normal_filter)]]
Rcpp::List dcc_constant_normal_filter(const arma::mat Z, const arma::mat S, const int n_update)
{
    int timesteps = Z.n_rows;
    int m = Z.n_cols;
    // not used
    Rcpp::String distribution = "mvn";
    Rcpp::String method = "pearson";
    arma::mat R = make_correlation(Z.head_rows(n_update), method);
    bool ispd = is_psd(R);
    if (!ispd) {
        R = make_psd(R);
    }
    arma::mat r_inverse = arma::inv_sympd(R);
    double part1 = arma::log_det_sympd(R);
    double const_part = m * log(2 * M_PI);
    double part2 = 0.0;
    double garchll = 0.0;
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);

    for (int i = 0;i<timesteps;i++) {
        part2 = arma::as_scalar(Z.row(i) * Z.row(i).t());
        garchll = 0.5 * (const_part + 2.0 * log(arma::prod(S.row(i))) + part2);
        dccll_vec(i) =  0.5 * (arma::as_scalar(Z.row(i) * (r_inverse  * Z.row(i).t())) + part1 - part2);
        ll_vec(i) = dccll_vec(i) + garchll;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("R") = R ,   _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec,  _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}

// [[Rcpp::export(.dcc_constant_student_filter)]]
Rcpp::List dcc_constant_student_filter(double shape, const arma::mat Z, arma::mat S, const int n_update)
{
    int timesteps = Z.n_rows;
    int m = Z.n_cols;
    // not used
    Rcpp::String distribution = "mvt";
    Rcpp::String method = "kendall";
    arma::mat R = make_correlation(Z.head_rows(n_update), method);
    R = transform_correlation(R, method);
    bool ispd = is_psd(R);
    if (!ispd) {
        R = make_psd(R);
    }
    double const_term = lgamma(0.5 * (m + shape)) - lgamma(0.5 * shape) - 0.5 * m * log(M_PI * (shape - 2.0));
    arma::mat r_inverse = arma::inv_sympd(R);
    double part1 = arma::log_det_sympd(R);
    double part2 = 0.0;
    double part3 = 0.0;
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);

    for (int i = 0;i<timesteps;i++) {
        part2 = arma::as_scalar(Z.row(i) * (r_inverse  * Z.row(i).t()));
        part3 = log(arma::prod(S.row(i)));
        dccll_vec(i) = const_term - 0.5 * part1 - 0.5 * (shape + m) * log(1.0 + (1.0/(shape - 2.0)) * part2);
        ll_vec(i) = dccll_vec(i) - part3;
        ll_vec(i) *= -1.0;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = -1.0 * arma::accu(dccll_vec);
    List L = List::create(Named("R") = R ,  _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec,  _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}

// [[Rcpp::export(.dcc_dynamic_normal_filter)]]
Rcpp::List dcc_dynamic_normal_filter(const arma::vec alpha, const arma::vec gamma, const arma::vec beta,
                                     const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder, const int n_update)
{
    const int alpha_order = dccorder[0];
    const int gamma_order = dccorder[1];
    const int beta_order = dccorder[2];
    const int maxpq = std::max(alpha_order, beta_order);
    const int m = z.n_cols;
    double sum_alpha = arma::accu(alpha);
    double sum_gamma = arma::accu(gamma);
    double sum_beta = arma::accu(beta);
    arma::mat zero_matrix = arma::zeros(maxpq, m);
    arma::mat Qbar = arma::cov(z.head_rows(n_update));
    arma::mat AsyZinit = arma::zeros(z.n_rows, m);
    arma::mat Nbar = arma::zeros(m, m);
    if (gamma_order > 0) {
        AsyZinit = matrix_sign(z) % z;
        Nbar = arma::cov(AsyZinit.head_rows(n_update));
    }
    const int timesteps = z.n_rows + maxpq;
    arma::mat Z = arma::join_cols(zero_matrix, z);
    const arma::mat S = arma::join_cols(zero_matrix, s);
    arma::mat AsyZ = arma::join_cols(zero_matrix, AsyZinit);
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);
    arma::mat Omega = Qbar * (1.0 - sum_alpha - sum_beta);
    if (gamma_order > 0) {
        Omega = Omega - sum_gamma * Nbar;
    }
    arma::cube C(m,m,timesteps);
    arma::cube R(m,m,timesteps);
    int i,j;
    for(i = 0;i<maxpq;i++){
        C.slice(i) = Qbar;
    }
    arma::mat Q = arma::zeros(m,m);

    double const_part = m * log(2 * M_PI);
    double part1 = 0.0;
    double part2 = 0.0;
    double garch_ll = 0.0;

    for (i = maxpq; i < timesteps; ++i) {
        Q = Omega;
        if (alpha_order > 0) {
            for(j = 0;j<alpha_order;++j) {
                Q = Q + alpha(j) * (Z.row(i - (j + 1)).t() * Z.row(i - (j + 1)));
            }
        }
        if (gamma_order > 0) {
            for(j = 0;j<gamma_order;++j) {
                Q = Q + gamma(j) * (AsyZ.row(i - (j + 1)).t() * AsyZ.row(i - (j + 1)));
            }
        }
        if (beta_order > 0) {
            for(j = 0;j<beta_order;++j) {
                Q = Q + beta(j) * C.slice(i - (j + 1));
            }
        }
        C.slice(i) = Q;
        arma::vec tempx = arma::sqrt(Q.diag());
        arma::mat tempy = tempx * tempx.t();
        arma::mat Rt = Q/tempy;
        R.slice(i) = Rt;
        arma::mat id_matrix = arma::eye(m, m);
        arma::mat Rinv = arma::inv(Rt);
        part1 = arma::log_det_sympd(Rt);
        part2 = arma::as_scalar(Z.row(i) * Z.row(i).t());
        garch_ll = 0.5 * (const_part + 2.0 * log(arma::prod(S.row(i))) + part2);
        dccll_vec(i) =  0.5 * (arma::as_scalar(Z.row(i) * (Rinv  * Z.row(i).t())) + part1 - part2);
        ll_vec(i) = dccll_vec(i) + garch_ll;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("Qbar") = Qbar, _("Nbar") = Nbar, _("R") = R , _("Q") = C,  _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec, _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}


// [[Rcpp::export(.dcc_dynamic_student_filter)]]
Rcpp::List dcc_dynamic_student_filter(const arma::vec alpha, const arma::vec gamma,
                                         const arma::vec beta, double shape,
                                         const arma::mat z, const arma::mat s, Rcpp::IntegerVector dccorder,
                                         const int n_update)
{
    const int alpha_order = dccorder[0];
    const int gamma_order = dccorder[1];
    const int beta_order = dccorder[2];
    const int maxpq = std::max(alpha_order, beta_order);
    const int m = z.n_cols;
    double sum_alpha = arma::accu(alpha);
    double sum_gamma = arma::accu(gamma);
    double sum_beta = arma::accu(beta);
    arma::mat zero_matrix = arma::zeros(maxpq, m);
    Rcpp::String distribution = "mvt";
    arma::mat Qbar = arma::cov(z.head_rows(n_update));
    arma::mat AsyZinit = arma::zeros(z.n_rows, m);
    arma::mat Nbar = arma::zeros(m, m);
    if (gamma_order > 0) {
        AsyZinit = matrix_sign(z) % z;
        Nbar = arma::cov(AsyZinit.head_rows(n_update));
    }
    const int timesteps = z.n_rows + maxpq;
    arma::mat Z = arma::join_cols(zero_matrix, z);
    arma::mat S = arma::join_cols(zero_matrix, s);
    arma::mat AsyZ = arma::join_cols(zero_matrix, AsyZinit);
    arma::vec ll_vec = arma::zeros(timesteps);
    arma::vec dccll_vec = arma::zeros(timesteps);
    arma::mat Omega = Qbar * (1.0 - sum_alpha - sum_beta);
    if (gamma_order > 0) {
        Omega = Omega - sum_gamma * Nbar;
    }
    arma::cube C(m,m,timesteps);
    arma::cube R(m,m,timesteps);
    double const_term = -1.0 * lgamma(0.5 * (m + shape)) + lgamma(0.5 * shape) + 0.5 * m * log(M_PI * (shape - 2.0));
    double part1 = 0.0;
    double part2 = 0.0;
    double part3 = 0.0;
    int i,j;
    for(i = 0;i<maxpq;i++){
        C.slice(i) = Qbar;
    }
    arma::mat Q = arma::zeros(m,m);
    for (i = maxpq; i < timesteps; ++i) {
        Q = Omega;
        if (alpha_order > 0) {
            for(j = 0;j<alpha_order;++j) {
                Q = Q + alpha(j) * (Z.row(i - (j + 1)).t() * Z.row(i - (j + 1)));
            }
        }
        if (gamma_order > 0) {
            for(j = 0;j<gamma_order;++j) {
                Q = Q + gamma(j) * (AsyZ.row(i - (j + 1)).t() * AsyZ.row(i - (j + 1)));
            }
        }
        if (beta_order > 0) {
            for(j = 0;j<beta_order;++j) {
                Q = Q + beta(j) * C.slice(i - (j + 1));
            }
        }
        C.slice(i) = Q;
        arma::vec tempx = arma::sqrt(Q.diag());
        arma::mat tempy = tempx * tempx.t();
        arma::mat Rt = Q/tempy;
        R.slice(i) = Rt;
        arma::mat Rinv = arma::inv(Rt);
        part1 = arma::log_det_sympd(Rt);
        part2 = arma::as_scalar(Z.row(i) * (Rinv  * Z.row(i).t()));
        part3 = log(arma::prod(S.row(i)));
        dccll_vec(i) = const_term + 0.5 * part1 + 0.5 * (shape + m) * log(1.0 + (1.0/(shape - 2.0)) * part2);
        ll_vec(i) = dccll_vec(i) + part3;
    }
    double nll = arma::accu(ll_vec);
    double dccnll = arma::accu(dccll_vec);
    List L = List::create(Named("Qbar") = Qbar, _("Nbar") = Nbar, _("R") = R , _("Q") = C,  _["dccll_vec"] = dccll_vec, _["ll_vec"] = ll_vec,  _["dcc_nll"] = dccnll, _["nll"] = nll);
    return L;
}


// [[Rcpp::export(.dcc_dynamic_simulate)]]
Rcpp::List dcc_dynamic_simulate(const arma::vec alpha, const arma::vec gamma,
                                   const arma::vec beta, double shape,
                                   const arma::mat Qbar, arma::mat Nbar,
                                   const arma::cube Qinit,
                                   const arma::mat Zinit,
                                   const arma::mat std_noise,
                                   const int timesteps,
                                   const int burn,
                                   Rcpp::IntegerVector dccorder,
                                   Rcpp::String distribution)
{
    const int alpha_order = dccorder[0];
    const int gamma_order = dccorder[1];
    const int beta_order = dccorder[2];
    double sum_alpha = arma::accu(alpha);
    double sum_gamma = arma::accu(gamma);
    double sum_beta = arma::accu(beta);
    const int maxpq = std::max(alpha_order, beta_order);
    const int m = Qbar.n_cols;
    const int exc = burn + maxpq;
    arma::mat AsyZ = arma::zeros(timesteps + maxpq, m);
    arma::mat AsyZinit = matrix_sign(Zinit) % Zinit;
    arma::mat Z = arma::zeros(timesteps + maxpq, m);
    arma::mat Omega = Qbar * (1.0 - sum_alpha - sum_beta);
    if (gamma_order > 0) {
        Omega = Omega - sum_gamma * Nbar;
    }
    arma::cube C(m, m, timesteps + maxpq);
    arma::cube R(m, m, timesteps + maxpq);
    arma::vec chisqrv = as<arma::vec>(wrap(Rcpp::rchisq(timesteps + maxpq, shape)));
    int i,j;
    for(int i = 0;i<maxpq;++i){
        Z.row(i) = Zinit.row(i);
        AsyZ.row(i) = AsyZinit.row(i);
        C.slice(i) = Qinit.slice(i);
    }
    arma::mat Q = arma::zeros(m,m);
    for(i = maxpq;i<(timesteps + maxpq);++i){
        Q = Omega;
        if (alpha_order > 0) {
            for(j = 0;j<alpha_order;++j) {
                Q = Q + alpha(j) * (Z.row(i - (j + 1)).t() * Z.row(i - (j + 1)));
            }
        }
        if (gamma_order > 0) {
            for(j = 0;j<gamma_order;++j) {
                Q = Q + gamma(j) * (AsyZ.row(i - (j + 1)).t() * AsyZ.row(i - (j + 1)));
            }
        }
        if (beta_order > 0) {
            for(j = 0;j<beta_order;++j) {
                Q = Q + beta(j) * C.slice(i - (j + 1));
            }
        }
        C.slice(i) = Q;
        arma::vec tempx = arma::sqrt(Q.diag());
        arma::mat tempy = tempx * tempx.t();
        arma::mat Rt = Q/tempy;
        Rt = arma::symmatu(Rt);
        R.slice(i) = Rt;
        arma::rowvec ztmp = std_noise.row(i);
        if (distribution == "mvt") {
            Z.row(i) = rmvt(Rt, ztmp, shape, (double) chisqrv(i));
        } else if (distribution == "mvn") {
            Z.row(i) = rmvnorm(Rt, ztmp);
        } else {
            Rf_error("dccsim: unknown distribution");
        }
        AsyZ.row(i) = matrix_sign(Z.row(i)) % Z.row(i);
    }
    if (exc > 0) {
        R.shed_slices(0, exc - 1);
        Z.shed_rows(0, exc - 1);
    }
    arma::mat RR = sym2tril(R, false);
    List L = List::create(Named("R") = RR, _("Z") = Z,  _("Q") = C, _("chisqrv") = chisqrv);
    return L;
}


// [[Rcpp::export(.dcc_constant_simulate)]]
Rcpp::List dcc_constant_simulate(const double shape, const arma::mat R, const arma::mat std_noise, const int timesteps, Rcpp::String distribution)
{
    int m = R.n_cols;
    arma::mat Z = arma::zeros(timesteps, m);
    arma::vec chisqrv = as<arma::vec>(wrap(Rcpp::rchisq(timesteps, shape)));
    for(int i = 0;i<timesteps;++i){
        arma::rowvec ztmp = std_noise.row(i);
        if (distribution == "mvt") {
            Z.row(i) = rmvt(R, ztmp, shape, (double) chisqrv(i));
        } else if (distribution == "mvn") {
            Z.row(i) = rmvnorm(R, ztmp);
        } else {
            Rf_error("cgarchsim: unknown distribution");
        }
    }
    List L = List::create(Named("R") = R, _("Z") = Z, _("chisqrv") = chisqrv);
    return L;
}
