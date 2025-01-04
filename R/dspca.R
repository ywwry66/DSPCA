new_dspca <- function(x) {
  stopifnot(is.list(x))
  stopifnot(identical(names(x),
                      c("newu", "n1", "n2", "mean1", "mean2",
                        "cov1", "cov2", "cov", "rot", "type")))
  structure(x, class = "dspca")
}

dspca_internal1 <- function(x1, x2, u1, u2, newu, rho, k,
                            h1_mean, h2_mean, h1_cov, h2_cov) {
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  mean1 <- mean_ker(x1, u1, newu, h1_mean)
  mean2 <- mean_ker(x2, u2, newu, h2_mean)
  cov1 <- cov_ker(x1, u1, newu, h1_cov, h1_mean)
  cov2 <- cov_ker(x2, u2, newu, h2_cov, h2_mean)
  cov <- (1 / (n1 + n2)) * (n1 * cov1 + n2 * cov2)
  rot <- rotmat(x1, x2, u1, u2, mean1, mean2, h1_cov, h2_cov, rho, k, newu)

  object <- list(ob1 = list(newu = newu, n1 = n1, n2 = n2,
                            mean1 = mean1, mean2 = mean2,
                            cov1 = cov1, cov2 = cov2, cov = cov),
                 rot = rot)
  return(object)
}

dspca_internal <- function(x1, x2, u1, u2, newu, rho, k,
                           h1_mean, h2_mean, h1_cov, h2_cov, type = "lda") {
  object <- dspca_internal1(x1, x2, u1, u2, newu, rho, k,
                            h1_mean, h2_mean, h1_cov, h2_cov)
  return(new_dspca(c(object$ob1, list(rot = object$rot[[1]], type = type))))
}

dspca <- function(x, u, y, newu, rho, k,
                  h1_mean, h2_mean, h1_cov, h2_cov, type = "lda") {
  x <- as.matrix(x)
  lvs <- levels(y)
  n1_ind <- which(y == lvs[1])
  n2_ind <- which(y == lvs[2])
  x1 <- x[n1_ind, ]
  x2 <- x[n2_ind, ]
  u1 <- u[n1_ind]
  u2 <- u[n2_ind]

  return(dspca_internal(x1, x2, u1, u2, newu, rho, k,
                        h1_mean, h2_mean, h1_cov, h2_cov, type))
}

predict_dspca_internal <- function(object, newx) {
  n1 <- object$n1
  n2 <- object$n2
  mean1 <- object$mean1
  mean2 <- object$mean2
  cov1 <- object$cov1
  cov2 <- object$cov2
  cov <- object$cov
  rot <- object$rot
  type <- object$type

  if (type == "lda") {
    ser <- t(rot) %*% cov %*% rot
    delr <- t(rot) %*% (mean1 - mean2)
    dscrm <- t(ginv(ser) %*% (delr)) %*% (t(rot) %*% newx) + log(n1 / n2) -
      (1 / 2) * t((t(rot) %*% (mean1 + mean2))) %*% ginv(ser) %*% delr
    dscrm <- if (all(abs(Im(dscrm)) < 1e-6)) as.double(dscrm)
    return(dscrm >= 0)
  } else if (type == "qda") {
    m1 <- t(rot) %*% mean1
    m2 <- t(rot) %*% mean2
    v1 <- t(rot) %*% cov1 %*% rot
    v2 <- t(rot) %*% cov2 %*% rot
    beta0 <- (-1 / 2) * (log(det(v1) / det(v2)) + t(m1) %*% ginv(v1) %*% m1 -
                           t(m2) %*% ginv(v2) %*% m2) - log(n2 / n1)
    beta <- ginv(v1) %*% m1 - ginv(v2) %*% m2
    ome <- (-1 / 2) * (ginv(v1) - ginv(v2))

    return(beta0 + t(beta) %*% (t(rot) %*% newx) + t(t(rot) %*% newx) %*%
             ome %*% (t(rot) %*% newx) >= 0)
  }
}

## predict.dspca <- function(object, newx, newu) {
##   if (object$newu != newu)
##     stop("Please use the same newu that is used to train the model!")

##   return(predict_dspca_internal(object, newx))
## }

## dspca <- function(newx, newu, x1, x2, u1, u2, rho, k,
##                   h1_mean, h2_mean, h1_cov, h2_cov,
##                   pre_trained = NULL, type = "lda") {
##   newx <- as.numeric(newx)
##   if (is.null(pre_trained)) {
##     n1 <- nrow(x1)
##     n2 <- nrow(x2)
##     me1 <- mean_ker(x1, u1, newu, h1_mean)
##     me2 <- mean_ker(x2, u2, newu, h2_mean)
##     se1 <- cov_ker(x1, u1, newu, h1_cov, h1_mean)
##     se2 <- cov_ker(x2, u2, newu, h2_cov, h2_mean)
##     se <- (1 / (n1 + n2)) * (n1 * se1 + n2 * se2)
##     rot <- rotmat(x1, x2, u1, u2, me1, me2, h1_cov, h2_cov, rho, k, newu)[[1]]
##   } else {
##     n1 <- pre_trained$n1
##     n2 <- pre_trained$n2
##     me1 <- pre_trained$me1
##     me2 <- pre_trained$me2
##     se1 <- pre_trained$se1
##     se2 <- pre_trained$se2
##     se <- pre_trained$se
##     rot <- pre_trained$rot
##   }
##   Ur <- rot

##   if (type == "lda") {
##     ser <- t(Ur) %*% se %*% Ur
##     delr <- t(Ur) %*% (me1 - me2)
##     dscrm <- (t(ginv(ser) %*% (delr)) %*% (t(Ur) %*% newx) + log(n1 / n2) -
##                 (1 / 2) * t((t(Ur) %*% (me1 + me2))) %*% ginv(ser) %*% delr)
##     dscrm <- if (all(abs(Im(dscrm)) < 1e-6)) as.double(dscrm)
##     return(dscrm >= 0)
##   } else if (type == "qda") {
##     Ume1 <- t(Ur) %*% me1
##     Ume2 <- t(Ur) %*% me2
##     Use1 <- t(Ur) %*% se1 %*% Ur
##     Use2 <- t(Ur) %*% se2 %*% Ur
##     beta0 <-
##       (-1 / 2) * (log(det(Use1) / det(Use2)) + t(Ume1) %*% ginv(Use1) %*% Ume1 -
##                   t(Ume2) %*% ginv(Use2) %*% Ume2) - log(n2 / n1)
##     beta <- ginv(Use1) %*% Ume1 - ginv(Use2) %*% Ume2
##     ome <- (-1 / 2) * (ginv(Use1) - ginv(Use2))

##     return(beta0 + t(beta) %*% (t(Ur) %*% newx) + t(t(Ur) %*% newx) %*%
##              ome %*% (t(Ur) %*% newx) >= 0)
##   }
## }
