simulate_data <- function(KK, TT, w, alpha, beta, a, b, inits = NULL) {
  set.seed(42)

  H <- rep(list(matrix(0, nrow = KK , ncol = KK)), times = TT)
  D <- rep(list(matrix(0, nrow = KK , ncol = KK)), times = TT)
  R <- rep(list(matrix(0, nrow = KK , ncol = KK)), times = TT)
  Q <- rep(list(matrix(0, nrow = KK , ncol = KK)), times = TT)

  y <- matrix(0, nrow = KK, ncol = TT)
  u <- matrix(0, nrow = KK, ncol = TT)

  R_init <- matrix(rnorm(KK^2), nrow = KK, ncol = KK)
  R_init <- t(R_init) %*% R_init
  Q_init <- matrix(rnorm(KK^2), nrow = KK, ncol = KK)
  Q_init <- t(Q_init) %*% Q_init
  if (is.null(inits)) {
    # H0 <- diag(1, KK)
    D0 <- diag(1, KK)
  } else {
    # H0 <- inits[[2]]
    D0 <- inits[[3]]
  }
  H0 <- D0 %*% R_init %*% D0
  y0 <- MASS::mvrnorm(1, rep(0, times = KK), Sigma = H0)
  u0 <- solve(D0) %*% y0
  Q0 <- matrix(rnorm(KK^2), nrow = KK, ncol = KK)
  Q0 <- t(Q0) %*% Q0
  #
  # browser()
  #
  diag(D[[1]]) <- w + alpha * y0^2 + beta * diag(H0)

  Q[[1]] <- (1 - a - b) * R_init + a * crossprod(u0) + b * Q0
  H[[1]] <- D[[1]] %*% R %*% D[[1]]
  y[, 1] <- MASS::mvrnorm(1, rep(0, times = KK), Sigma = H[[1]])

  for (t in 2:TT) {
    diag(D[[t]]) <- w + alpha * y[, t - 1]^2 + beta * diag(H[[t - 1]])
    Q[[t]] <- (1 - a - b) * R_init + a * crossprod(u[, t - 1]) + b * Q[[t - 1]]
    H[[t]] <- D[[t]] %*% R %*% D[[t]]
    y[, t] <- MASS::mvrnorm(1, rep(0, times = KK), Sigma = H[[t]])
  }
  return(list(y, ))
}