KK <- 2
TT <- 3

w     <- rep(0.5, times = KK)
alpha <- rep(0.1, times = KK)
beta  <- rep(0.1, times = KK)
a <- 0.35
b <- 0.45

bayesDCC::simulate_data(KK, TT, w, alpha, beta, a, b, inits = NULL)
