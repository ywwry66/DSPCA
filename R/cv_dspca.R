new_cv_dspca <- function(x) {
  stopifnot(is.list(x))
  stopifnot(identical(names(x),
                      c("x1", "x2", "u1", "u2", "h1_mean", "h2_mean",
                        "h1_cov", "h2_cov", "rho_opt", "k_opt", "lvs")))
  structure(x, class = "cv_dspca")
}

##' Runs K-fold cross-validation for DSPCA.
##'
##' Runs K-fold cross-validation to select the optimal values of
##' `h_mean`, `h_cov`, `rho` and `k`. Returns an object of class
##' `"cv_dspca"` that contains the optimal tuning parameters,
##' including the bandwiths for Nadaraya-Watson estimator, the weight
##' in the total covariance, and the number of eigenvectors used for
##' projection. Call the `predict` function with the returned object
##' to predict the class label of new observations.
##'
##' @title Cross-validation for DSPCA
##' @param x the input matrix. Each column consists of observations
##'   from an input variable.
##' @param u a numeric vector, consisting of observations from the
##'   index variable.
##' @param y a two-level factor, consisting of observations from the
##'   response variable. Levels represent class labels.
##' @param n_folds number of folds used for cross-validation.
##' @param type the classification method used in conjunction with the
##'   dynamic supervised principal component analysis framework.
##'   Supported methods are `type = "lda"` for linear discriminant
##'   analysis and `type = "qda"` for quadratic discriminant analysis.
##' @param h_mean a numeric vector of kernel bandwidths for the
##'   Nadaraya-Watson estimator of the means, to be chosen by
##'   leave-one-out cross-validation.
##' @param h_cov a numeric vector of kernel bandwidths for the
##'   Nadaraya-Watson estimator of the covariances, to be chosen by
##'   leave-one-out cross-validation.
##' @param rho a numeric vector of weight parameters for constructing
##'   the total covariance, to be chosen by cross-validation.
##' @param k the number of eigenvectors used for projection.
##' @return An object of S3 class `"cv_dspca"`, which is a list with
##'   the following components:
##'
##' \item{x1}{the submatrix of `x` with observations from the first
##' class.}
##'
##' \item{x2}{the submatrix of `x` with observations from the second
##' class.}
##'
##' \item{u1}{the subvector of `u` with observations from the first
##' class.}
##'
##' \item{u2}{the subvector of `u` with observations from the second
##' class.}
##'
##' \item{h1_mean}{the optimal kernel bandwidth for the
##' Nadaraya-Watson estimator of the first mean chosen by
##' leave-one-out cross-validation.}
##'
##' \item{h2_mean}{the optimal kernel bandwidth for the
##' Nadaraya-Watson estimator of the second mean chosen by
##' leave-one-out cross-validation.}
##'
##' \item{h1_cov}{the optimal kernel bandwidth for the Nadaraya-Watson
##' estimator of the first covariance chosen by leave-one-out
##' cross-validation.}
##'
##' \item{h2_cov}{the optimal kernel bandwidth for the Nadaraya-Watson
##' estimator of the second covariance chosen by leave-one-out
##' cross-validation.}
##'
##' \item{rho_opt}{the optimal weight parameter for constructing the
##' total covariance chosen by cross-validation.}
##'
##' \item{k_opt}{the optimal number of eigenvectors used for
##' projection chosen by cross-validation.}
##'
##' \item{lvs}{the two levels of the factor `y`.}
##' @author Wenbo Ouyang, Ruiyang Wu, Ning Hao, Hao Helen Zhang
##' @export
cv_dspca <- function(x, u, y, n_folds = 5, type = "lda",
                     h_mean = (1:15) / 10, h_cov = (1:15) / 10,
                     rho = exp(-1:6), k = 1:5) {
  x <- as.matrix(x)
  lvs <- levels(y)
  n1_ind <- which(y == lvs[1])
  n2_ind <- which(y == lvs[2])
  n1 <- length(n1_ind)
  n2 <- length(n2_ind)
  x1 <- x[n1_ind, ]
  x2 <- x[n2_ind, ]
  u1 <- u[n1_ind]
  u2 <- u[n2_ind]

  h1_mean <- cv_h_mean(h_mean, x1, u1)
  h2_mean <- cv_h_mean(h_mean, x2, u2)
  h1_cov <- cv_h_cov(h_cov, x1, u1, h1_mean)
  h2_cov <- cv_h_cov(h_cov, x2, u2, h2_mean)

  folds1 <- split(sample(1:n1), rep(1:n_folds, length = n1))
  folds2 <- split(sample(1:n2), rep(1:n_folds, length = n2))

  err <- matrix(0, nrow = length(rho), ncol = length(k))
  for (i in seq_len(n_folds)) {
    x1_cv <- x1[-folds1[[i]], ]
    x2_cv <- x2[-folds2[[i]], ]
    u1_cv <- u1[-folds1[[i]]]
    u2_cv <- u2[-folds2[[i]]]

    for (j in folds1[[i]]) {
      dspca_obj <- dspca_internal1(x1_cv, x2_cv, u1_cv, u2_cv, u1[j], rho,
                                   max(k), h1_mean, h2_mean, h1_cov, h2_cov)
      for (l in seq_along(rho)) {
        for (r in seq_along(k)) {
          obj <- c(dspca_obj$ob1,
                   list(rot = dspca_obj$rot[[l]][, 1:k[r]], type = type))
          err[l, r] <- err[l, r] +
            !predict_dspca_internal(new_dspca(obj), x1[j, ])
        }
      }
    }
    for (j in folds2[[i]]) {
      dspca_obj <- dspca_internal1(x1_cv, x2_cv, u1_cv, u2_cv, u2[j], rho,
                                   max(k), h1_mean, h2_mean, h1_cov, h2_cov)
      for (l in seq_along(rho)) {
        for (r in seq_along(k)) {
          obj <- c(dspca_obj$ob1,
                   list(rot = dspca_obj$rot[[l]][, 1:k[r]], type = type))
          err[l, r] <- err[l, r] +
            predict_dspca_internal(new_dspca(obj), x2[j, ])
        }
      }
    }
  }
  opt <- arrayInd(which.min(err), dim(err))
  rho_opt <- rho[opt[1]]
  k_opt <- k[opt[2]]
  object <- list(x1 = x1, x2 = x2, u1 = u1, u2 = u2,
                 h1_mean = h1_mean, h2_mean = h2_mean,
                 h1_cov = h1_cov, h2_cov = h2_cov,
                 rho_opt = rho_opt, k_opt = k_opt, lvs = lvs)
  return(new_cv_dspca(object))
}

##' This function is used to predict the class labels for new
##' observations
##'
##' `predict.cv_dspca` outputs the predicted class labels for new
##' observations. The new values of the input variables are passed in
##' `newx`. The new values of the index variable are passed in `newu`.
##'
##' @title DSPCA Prediction
##' @param object fitted object from the `cv_dspca` function.
##' @param newx a matrix of new values of `x`. Each row is a new
##'   observation whose class label is to be precdicted.
##' @param newu a vector of new values of the index variable `u`.
##' @param type the classification method used in conjunction with the
##'   dynamic supervised principal component analysis framework.
##'   Supported methods are `type = "lda"` for linear discriminant
##'   analysis and `type = "qda"` for quadratic discriminant analysis.
##' @return Predicted class labels for `newx` and `newu`.
##' @author Wenbo Ouyang, Ruiyang Wu, Ning Hao, Hao Helen Zhang
##' @export
predict.cv_dspca <- function(object, newx, newu, type = "lda") {
  newx <- as.matrix(newx)
  n_new <- nrow(newx)
  sn <- numeric(n_new)
  for (i in seq_len(n_new)) {
    dspca_obj <- dspca_internal(object$x1, object$x2, object$u1, object$u2,
                                newu[i], object$rho_opt, object$k_opt,
                                object$h1_mean, object$h2_mean,
                                object$h1_cov, object$h2_cov, type = type)
    sn[i] <- predict_dspca_internal(dspca_obj, newx[i, ])
  }
  return(object$lvs[2 - sn])
}
