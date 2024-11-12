.dmvnorm <- function(x, mu, sigma, log = FALSE) {
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mu)) {
        mean <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(1, ncol(x), ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mu) != NROW(sigma)) {
        stop("mu and sigma have non-conforming size")
    }
    distval <- mahalanobis(x, center = mu, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    value <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
    if (!log) value <- exp(value)
    return(value)
}

.dmvt <- function(x, mu, sigma, shape = 1, log = FALSE)
{
    if (shape < 2.00001) stop("\nshape must be greater than 2.0")
    if (is.vector(x)) {
        x <- matrix(x, ncol = length(x))
    }
    if (missing(mu)) {
        mu <- rep(0, length = ncol(x))
    }
    if (missing(sigma)) {
        sigma <- diag(1, ncol(x), ncol(x))
    }
    if (NCOL(x) != NCOL(sigma)) {
        stop("x and sigma have non-conforming size")
    }
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mu) != NROW(sigma)) {
        stop("mu and sigma have non-conforming size")
    }
    m <- NCOL(sigma)
    distval <- mahalanobis(x, center = mu, cov = sigma)
    logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
    logretval <- lgamma((m + shape)/2) - (lgamma(shape/2) + 0.5 * (logdet + m * logb(pi * shape))) - 0.5 * (shape + m) * logb(1 + distval/shape)
    if (log) retval = logretval else retval = exp(logretval)
    return(retval)
}


rmvnorm <- function(n = 1, mu = rep(0, ncol(sigma)), sigma)
{
    R <- cov2cor(sigma)
    m <- NCOL(R)
    Z <- matrix(rnorm(n * m), ncol = m, nrow = n)
    s <- sqrt(diag(sigma))
    out <- .rmvnorm(R, Z)
    out <- sweep(out, 2, s, FUN = "*")
    out <- sweep(out, 2, mu, FUN = "+")
    return(out)
}

dcopula_gauss <- function(u, s, log = FALSE) {
    m <- dim(s)[2]
    z <- apply(u, 2, FUN = function(x) qnorm(x))
    value <- .dmvnorm(z, mu = rep(0, m), sigma = s, log = TRUE) - apply(dnorm(z, log = TRUE), 1,  "sum")
    if (!log) value <- exp(value)
    return(value)
}

rcopula_gauss <- function(n, u, s)
{
    m <- dim(s)[2]
    mu <- rep(0, m)
    value <- pnorm(rmvnorm(n, mu = mu, sigma = s))
    return(value)
}


# .interpolate <- function(x, y, z) {
#     ans = approx(x, y = y, xout = z, method = "linear", n = length(z),
#                  yleft = min(z), yright = max(z), rule = 1, f = 0, ties = mean)$y
#     return(ans)
# }
#
#
# cfinv <- function(z, f, step = 0.01, ...)
# {
#     pmax <- 18
#     p <- 14
#     maxz <- round(max(abs(z))) + 5
#     while ((maxz/step + 1) > 2^(p - 1)) {
#         p <- p + 1
#     }
#     if (p > pmax) p <- pmax
#     if ((maxz/step + 1) > 2^(p - 1)) {
#         step <- (maxz + 1) * (1 + step/10)/(2^(p - 1))
#     }
#     zs <- sort(z)
#     # run the fast fourier transform
#     n <- 2^p
#     x <- seq(0, n - 1, by = 1) * step - (n * step/2)
#     s <- 1/(step * n)
#     tt <- 2 * pi * s * (seq(0, n - 1, by = 1) - n/2)
#     sgn <- rep(1, n)
#     ds <- seq(2, n, by = 2)
#     sgn[ds] <- -1 * rep(1, n/2)
#     cf <- f(tt, ...)
#     phi <- sgn * cf
#     phi[n/2 + 1] <- sgn[n/2 + 1]
#     p <- s * abs(fft(phi))
#     pdf <- .interpolate_window(x, p, zs, w = -1)
#     return(pdf)
# }
#
# # characteristic function of nig with independent margins
# nigmvcf <- function(z, alpha, beta, delta, mu)
# {
#     N <- length(z)
#     m <- length(mu)
#     x1 <- 1i * z * sum(mu)
#     zx <- matrix(0, ncol = m, nrow = N)
#     zx <- apply(cbind(delta, alpha, beta), 1, FUN = function(x) x[1] * (sqrt(x[2]^2 - x[3]^2) - sqrt(x[2]^2 - (x[3] + 1i * z)^2)))
#     x2 <- apply(t(zx), 2, FUN = function(x) sum(x))
#     ans <- x1 + x2
#     return(exp(ans))
# }
#
# .ghypfn <- function(lambda, alpha, beta, delta, z)
# {
#     x1 <- (lambda/2) * (log(alpha^2 - beta^2) - log(alpha^2 - (beta + 1i * z)^2))
#     x2 <- log(bessel_k(delta * sqrt(alpha^2 - (beta + 1i * z)^2), abs(lambda))) - log(bessel_k(delta * sqrt(alpha^2 - beta^2), abs(lambda)))
#     return(x1 + x2)
# }
#
# ghypmvcf <- function(z, lambda, alpha, beta, delta, mu)
# {
#     N <- length(z)
#     m <- length(mu)
#     x0 <- 1i * z * sum(mu)
#     zx <- matrix(0, ncol = m, nrow = N)
#     zx <- apply(cbind(lambda, alpha, beta, delta), 1, FUN = function(x) .ghypfn(x[1], x[2], x[3], x[4], z))
#     x <- apply(t(zx), 2, FUN = function(x) sum(x))
#     ans <- exp(x0 + x)
#     return(ans)
# }


# .dfft <- function(object, index = 1, type = "predict")
# {
#     if (type %in% c("predict","backtest")) {
#         if (object$h == 1) index <- index + 1
#     }
#     switch(object$distribution,
#            mvnorm = .dfft_mvnorm(object, index),
#            manig  = .dfft_numerical(object, index),
#            magh   = .dfft_numerical(object, index))
# }



.csimpsum <- function(fx) {
    l <- length(fx)
    l2 <- l %/% 2
    if (l %% 2 == 0) {
        fx <- c(fx[1:l2], (fx[l2] + fx[l2 + 1])/2, fx[(l2 + 1):l])
        l <- l + 1
    }
    f_even <- fx[seq(l) %% 2 == 0]
    f_odd  <- fx[seq(l) %% 2 == 1]
    fs <- 2 * cumsum(f_odd) - f_odd - f_odd[1]
    fsm <- 4 * cumsum(f_even)
    ff <- c(0, (fs[2:(l2 + 1)] + fsm)/3)
    return(ff)
}

.cdf_fun <- function(x, y, yleft, yright){
    stepfun(x = x, y = c(yleft, y))
}

.icdf_fun <- function(x, y, y_left, y_right){
    y_left <- y_left[1]
    y_right <- y_right[1]
    yl <- if (is.finite(y_left)) y_left  else y[1]
    yr <- if (is.finite(y_right)) y_right else y[length(y)]
    f1 <- approxfun(x = x, y = y, method = "linear", yleft = yl, yright = yr, ties = "ordered")
    f <- function(x) {
        y1 <- f1(x)
        y1[.is_equal(x,0)] <- y_left
        y1[.is_equal(x,1)] <- y_right
        return(y1)
    }
    return(f)
}

.create_pdf <- function(x, dx, h = NULL, standM = "sum", mu = NULL) {
    if (is.null(mu)) mu <- 0
    dx <- (dx >= .Machine$double.eps) * dx
    if (length(dx) < length(x)) dx <- c(0, dx)
    if (is.null(h)) h <- 1
    dx1 <- dx / h
    ## density
    df1 <- approxfun(x = x + mu, y = dx1, yleft = 0, yright = 0, ties = "ordered")
    if (standM == "sum") {
        stand <- sum(dx)
    } else {
        stand <- try(integrate(df1, -Inf, Inf)$value, TRUE)
        if (is(stand,"try-error")) {
            warning("'integrate' function threw an error. Results may be inaccurate.")
            stand <- sum(df1(x)) * h * (x[2] - x[1])
        }
    }
    dfun <- function(x, log = FALSE) {
        if (log) {
            d0 <- log(df1(x)) - log(stand)
        } else {
            d0 <- df1(x) / stand
        }
        return(d0)
    }
    return(dfun)
}

.create_cdf <- function(x, dx, h = NULL, mu = NULL) {
    if (is.null(mu)) mu <- 0
    if (is.null(h)) h <- 0
    x_upper <- x_lower <- x
    l <- length(x)
    if ((l %% 2 == 0)) {
        l2 <- l/2
        x_lower <- c(x[1:l2], (x[l2] + x[l2 + 1])/2, x[(l2 + 1):l])
        x_upper <- c(x[1:l2], (x[l2] + x[l2 + 1])/2, x[(l2 + 1):l])
        l <- l + 1
    }
    x_lower <- x_lower[seq(l) %% 2 == 1]
    x_upper <- x_upper[seq(l) %% 2 == 1]
    p_lower <- .csimpsum(dx)
    nm_lower <- max(p_lower + mu)
    p1_lower <- approxfun(x = x_lower + mu, y = p_lower, yleft = 0, yright = nm_lower, ties = "ordered")
    nm_upper <- p1_lower(max(x + mu))
    p_upper <- rev(.csimpsum(rev(dx)))
    ## continuity correction by h/2
    p1_upper <- approxfun(x = x_upper + mu, y = p_upper, yright = 0, yleft = nm_upper, ties = "ordered")
    pfun <- function(q, lower_tail = TRUE, log_p = FALSE){
        if (lower_tail) {
            p0 <- p1_lower(q)
            p0 <- if (log_p) log(p0) - log(nm_lower) else p0/nm_lower
        } else {
            p0 <- p1_upper(q)
            p0 <- if (log_p) log(p0) - log(nm_upper) else p0/nm_upper
        }
        return(p0)
    }
    return(pfun)
}

.create_icdf <- function(x, dx, h, mu = NULL)
{
    if (is.null(mu)) mu <- 0
    cdf_f <- .create_cdf(x, dx, h, mu = mu)
    y_left <- min(x + mu)
    y_right <- max(x + mu)
    px_lower <- cdf_f(x)
    px_upper <- cdf_f(x, lower_tail = FALSE)
    ix_lower <- .is_equal_01(px_lower)
    xx_lower <- px_lower[!ix_lower]
    yy_lower <- x[!ix_lower]
    q_lower <- .icdf_fun(x = xx_lower, y = yy_lower, y_left = y_left, y_right = y_right)
    x_rev <- rev(x)
    px_upper <- rev(px_upper)
    ix_upper <- .is_equal_01(px_upper)
    xx_upper <- px_upper[!ix_upper]
    yy_upper <- x_rev[!ix_upper]
    q_upper <- .icdf_fun(x = xx_upper, y = yy_upper, y_left = y_right, y_right = y_left)
    p <- NULL
    qfun <- function(p, lower_tail = TRUE, log_p = FALSE) {
        if (log_p) p <- exp(p)
        if (any((p < -.Machine$double.eps) | (p > 1 + .Machine$double.eps))) warning("q method produced NaN's")
        i01 <- (-.Machine$double.eps <= p) & (p <= 1 + .Machine$double.eps)
        p01 <- p[i01]
        q0  <- p * 0
        q0[!i01] <- NaN
        if (lower_tail) {
            q0[i01] <- q_lower(p01)
        } else {
            q0[i01] <- q_upper(p01)
        }
        return(as.numeric(q0))
    }
    return(qfun)
}
