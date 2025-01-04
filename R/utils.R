##' Gaussian kernel
##' h: bandwidth. u: index variable.
ker <- function(h, u) {
  p <- (1 / h) * exp(-0.5 * (u / h)^2)
  return(p)
}

##' Kernel estimated mean when the index variable equals newu
##' x: observed predictors of the training sample.
##' u: observed index variable of the training sample.
##' h: bandwidth.
mean_ker <- function(x, u, newu, h) {
  weights <- ker(h, u - newu)
  mean_x <- as.vector(t(weights) %*% x / sum(weights))
  return(mean_x)
}

##' Kernel estimated covariance when the index variable equals newu
##' x: observed predictors of the training sample.
##' u: observed index variable of the training sample.
##' h_cov: bandwidth for variance estimation.
##' h_mean: bandwidth for mean estimation
cov_ker <- function(x, u, newu, h_cov, h_mean = NULL) {
  if (is.null(h_mean)) h_mean <- h_cov
  mean_x <- mean_ker(x, u, newu, h_mean)
  x_centered <- t(t(x) - mean_x)
  weights <- ker(h_cov, u - newu)
  cov_x <- t(x_centered) %*% diag(weights, length(weights)) %*%
    x_centered / sum(weights)
  return(cov_x)
}

##' Leave one out cross validation for choosing the best bandwidth in
##' the kernel estimation of the mean.
##' h_mean: a range of bandwidths to be chosen from
##' x: observed predictors of the training sample.
##' u: observed index variable of the training sample.
cv_h_mean <- function(h_mean, x, u) {
  err <- rep(0, length(h_mean))
  for (i in seq_along(h_mean)) {
    for (j in seq_along(u)) {
      err[i] <- err[i] +
        sum((x[j, ] - mean_ker(x[-j, ], u[-j], u[j], h_mean[i]))^2)
    }
  }
  return(h_mean[which.min(err)])
}

##' Leave one out cross validation for choosing the best bandwidth in
##' the kernel estimation of the covariance.
##' h_cov: a range of bandwidths to be chosen from
##' u: observed index variable of the training sample.
##' x: observed predictors of the training sample.
##' h_mean: the bandwidth for mean estimation.
cv_h_cov <- function(h_cov, x, u, h_mean) {
  err <- rep(0, length(h_cov))
  for (i in seq_along(h_cov)) {
    for (j in seq_along(u)) {
      cov <- cov_ker(x[-j, ], u[-j], u[j], h_cov[i], h_mean)
      mean <- mean_ker(x, u, u[j], h_mean)
      err[i] <- err[i] +
        norm((x[j, ] - mean) %*% t(x[j, ] - mean) - cov, "F")^2
    }
  }
  return(h_cov[which.min(err)])
}

##' Moore–Penrose inverse
ginv <- function(x, tol = sqrt(.Machine$double.eps)) {
  svd <- svd(x)
  d <- svd$d
  d <- d * (d > tol)
  d[which(d == 0)] <- Inf
  u <- svd$u
  v <- svd$v
  return(v %*% diag(1 / d, length(d)) %*% t(u))
}
