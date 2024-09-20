.adjust_breaks <- function(x, min_val, max_val)
{
    dens <- x$density
    brks <- x$breaks
    min_brk <- min(brks)
    max_brk <- max(brks)
    if (min_brk < min_val) {
        lower_index <- max(which(brks < min_val))
        dens <- dens[-(1:lower_index)]
        brks <-  brks[-(1:lower_index)]
    } else if (min_brk > min_val) {
        d <- median(diff(x$breaks))
        augment <- seq(min_val, min_brk, by = d)
        brks <- c(augment, brks[-1])
        dens <- c(rep(0, length(augment) - 1), dens)
    }

    if (max_brk > max_val) {
        n <- length(brks)
        upper_index <- min(which(x$breaks > max_val))
        dens <- dens[-(upper_index:n)]
        brks <-  brks[-(upper_index:n)]
    } else if (max_brk < max_val) {
        d <- median(diff(x$breaks))
        augment <- seq(max_brk, max_val, by = d)
        brks <- c(brks, augment[-1])
        dens <- c(dens, rep(0, length(augment) - 1))
    }
    return(list(density = dens, breaks = brks))
}


.constant_bivariate_correlation_plot <- function(x, pair = c(1,2), xlim = c(-5, 5), ...)
{
    parameter <- NULL
    min_val <- xlim[1]
    max_val <- xlim[2]
    cnames <- names(x$spec$univariate)
    # Set the desired number of breaks
    R <- tscor(x)[(pair),(pair)]
    if (inherits(x, "cgarch.estimate")) {
        distribution <- x$spec$copula
    } else {
        distribution <- x$spec$distribution
    }
    if (distribution == "gaussian") {
        f <- function(x, y) .dmvnorm(cbind(x,y), mu = rep(0, 2), sigma = R)
    } else {
        shape <- x$parmatrix[parameter == "shape"]$value
        f <- function(x, y) .dmvt(cbind(x,y), mu = rep(0, 2), sigma = R, shape = shape, log = FALSE)
    }
    vx <- seq(min_val, max_val, length = 401)
    vy <- seq(min_val, max_val, length = 401)
    z <- outer(vx,vy,f)
    u_distribution <- c(x$spec$univariate[[pair[1]]]$spec$distribution, x$spec$univariate[[pair[2]]]$spec$distribution)
    skew <- c(x$spec$univariate[[pair[1]]]$parmatrix[parameter == "skew"]$value,
              x$spec$univariate[[pair[2]]]$parmatrix[parameter == "skew"]$value)
    shape <- c(x$spec$univariate[[pair[1]]]$parmatrix[parameter == "shape"]$value,
               x$spec$univariate[[pair[2]]]$parmatrix[parameter == "shape"]$value)
    lambda <- c(x$spec$univariate[[pair[1]]]$parmatrix[parameter == "lambda"]$value,
                x$spec$univariate[[pair[2]]]$parmatrix[parameter == "lambda"]$value)
    d1 <- function(y) ddist(u_distribution[1], y, 0, 1, skew = skew[1], shape = shape[1], lambda = lambda[1])
    d2 <- function(y) ddist(u_distribution[2], y, 0, 1, skew = skew[2], shape = shape[2], lambda = lambda[2])
    X <- cbind(as.numeric(residuals(x$spec$univariate[[pair[1]]], standardize = TRUE)),
               as.numeric(residuals(x$spec$univariate[[pair[2]]], standardize = TRUE)))
    # calculate breakpoints
    x_hist <- hist(X[,1], plot = FALSE, breaks = "fd")
    x_d <- .adjust_breaks(x_hist, min_val, max_val)
    y_hist <- hist(X[,2], plot = FALSE, breaks = "fd")
    y_d <- .adjust_breaks(y_hist, min_val, max_val)
    top <- max(c(x_d$density, y_d$density, dnorm(0)))
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE), c(3,1), c(1,3), TRUE)
    par(mar = c(3, 3, 0, 0))
    image(vx, vy, z, col = rev(heat.colors(101)))
    contour(vx, vy, z, col = "blue", add = TRUE)
    points(X,cex = .2)
    title(ylab = cnames[pair[2]], xlab = cnames[pair[1]], line = -1)
    par(mar = c(0,3,1,1))
    barplot(x_d$density, axes = FALSE, ylim = c(0, top), space = 0, col = "steelblue")
    lines((density(X[,1])$x - x_d$breaks[1])/diff(x_d$breaks)[1], d1(density(X[,1])$x),col = "forestgreen")
    par(mar = c(3,0,1,1))
    barplot(y_d$density, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE, col = "steelblue")
    lines(d2(density(X[,2])$x),(density(X[,2])$x - y_d$breaks[1])/diff(y_d$breaks)[1], col = "forestgreen")
    return(invisible(x))
}


## NEED TO FIX
.dynamic_pairwise_correlation_plot <- function(x, index_format = "%Y", col = hcl.colors(12, "YlOrRd", rev = TRUE), ...)
{
    tf <- function(m) m[, ncol(m):1, drop = FALSE]
    R <- x$R
    R <- R[,order(colMeans(R), decreasing = TRUE)]
    d_index <- format(x$spec$target$index, index_format)
    col_key <- data.table("color" = col, "value" = seq(from = min(z), to = max(z), along.with = 1:length(col)))
    par.orig <- par(no.readonly = TRUE)
    on.exit(par(par.orig))
    par(fig = c(0,.9,0,1), mar = c(2,2,2,0))
    image(z = tf(R), col = col, axes = FALSE)
    at_index <- 1:NROW(R)/NROW(R)
    axis(1, at = at_index, labels = d_index, tick = F)
    title(main = "(Pairwise) Time Varying Correlation Heatmap")
    par(fig = c(.9,1,.3,.7), mar = c(1,1,1,2.5), new = T)
    image(t(matrix(col_key$value, ncol = 1)), col = col_key$color, axes = FALSE)
    axis(side = 4, lwd = 0, las = 2, line = -.75, cex.axis = 0.8)
    mtext(text = "", adj = 0, line = 1, cex = 0.8)
    return(invisible(R))
}


