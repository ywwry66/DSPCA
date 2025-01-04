rotmat <- function(x1, x2, u1, u2, mux1, mux2, hvarx1, hvarx2,
                   rho, k = min(dim(x1), dim(x2)), un) {
  rotation_matrix <- list()
  x1cen <- t(t(x1) - as.numeric(mux1))
  sigxn11 <- ker(hvarx1, u1 - un)
  x2cen <- t(t(x2) - as.numeric(mux2))
  sigyn11 <- ker(hvarx2, u2 - un)
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  p <- dim(x1)[2]
  XX <- rbind(sqrt(n1 / (n1 + n2)) * t(t(x1cen) %*% sqrt(diag(c(sigxn11)))) /
                sqrt(sum(sigxn11)),
              sqrt(n2 / (n1 + n2)) * t(t(x2cen) %*% sqrt(diag(c(sigyn11)))) /
                sqrt(sum(sigyn11)))
  delta <- mux1 - mux2
  M <- length(rho)
  if (p < n1 + n2) {
    for (m in 1:M) {
      W <- crossprod(XX)
      B <- crossprod(t(delta))
      Sigma_rho <- W + rho[m] * B
      current_rotation_matrix <-
        (eigen(Sigma_rho, symmetric = TRUE)$vectors)[, 1:k]
      rotation_matrix <- c(rotation_matrix, list(current_rotation_matrix))
    }
  } else {
    for (m in 1:M) {
      NewX <- rbind(XX, t(sqrt(rho[m]) * delta))
      small_Sigma_rho <- tcrossprod(NewX)
      temp <- (eigen(small_Sigma_rho, symmetric = TRUE)$vectors)[, 1:k]
      current_rotation_matrix <- scale(t(NewX) %*% temp,
                                       center = FALSE) / sqrt(p - 1)
      rotation_matrix <- c(rotation_matrix, list(current_rotation_matrix))
    }
  }
  return(rotation_matrix)
}
