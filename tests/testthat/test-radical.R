test_that("RADICAL [DC] Check",{
    X <- coredata(y)
    n_comp <- 5
    pca_method <- "ML"
    ic <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, seed = 100, trace = FALSE)
    mu <- ic$mu
    X_demeaned <- scale(X, center = mu, scale = FALSE)
    C <- (t(X_demeaned) %*% X_demeaned)/dim(X_demeaned)[1]
    ed <- eigen(C)
    D <- diag(ed$values)
    E <- ed$vectors
    D <- D[1:n_comp, 1:n_comp]
    E <- E[, 1:n_comp]
    K <- solve(sqrt(D)) %*% t(E)
    L <- E %*% sqrt(D)
    Z <- X_demeaned %*% t(K)
    W <- t(K) %*% ic$U
    S <- X_demeaned %*% W + matrix(1, nrow = nrow(X_demeaned)) %*% matrix(mu, nrow = 1) %*% W
    expect_equal(S, ic$S)
})

test_that("RADICAL [DC] Custom Input",{
    X <- coredata(y)
    n_comp <- 5
    pca_method <- "ML"
    mu <- colMeans(X)
    X_demeaned <- scale(X, center = mu, scale = FALSE)
    C <- (t(X_demeaned) %*% X_demeaned)/dim(X_demeaned)[1]
    ed <- eigen(C)
    D <- diag(ed$values)
    E <- ed$vectors
    D <- D[1:n_comp, 1:n_comp]
    E <- E[, 1:n_comp]
    K <- solve(sqrt(D)) %*% t(E)
    L <- E %*% sqrt(D)
    Z <- X_demeaned %*% t(K)
    ic_default <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, seed = 100, trace = FALSE)
    ic_custom_1 <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, E = E, D = D, seed = 100, trace = FALSE)
    expect_equal(ic_default$S, ic_custom_1$S)
    ic_custom_2 <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, Z = Z, K = K, L = L, seed = 100, trace = FALSE)
    expect_equal(ic_default$S, ic_custom_2$S)
})

test_that("RADICAL Custom Input",{
    X <- coredata(y)[,1:8]
    n_comp <- NCOL(X)
    pca_method <- "ML"
    mu <- colMeans(X)
    X_demeaned <- scale(X, center = mu, scale = FALSE)
    C <- (t(X_demeaned) %*% X_demeaned)/dim(X_demeaned)[1]
    ed <- eigen(C)
    D <- diag(ed$values)
    E <- ed$vectors
    D <- D[1:n_comp, 1:n_comp]
    E <- E[, 1:n_comp]
    K <- solve(sqrt(D)) %*% t(E)
    L <- E %*% sqrt(D)
    Z <- X_demeaned %*% t(K)
    ic_default <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, seed = 100, trace = FALSE)
    ic_custom_1 <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, E = E, D = D, seed = 100, trace = FALSE)
    expect_equal(ic_default$S, ic_custom_1$S)
    ic_custom_2 <- radical(X, components = n_comp, demean = TRUE, pca_cov = pca_method, Z = Z, K = K, L = L, seed = 100, trace = FALSE)
    expect_equal(ic_default$S, ic_custom_2$S)
})
