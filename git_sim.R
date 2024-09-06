# R code to simulate from bivariate copula MPCMP with (0, 0) inflation
# Load in look-up table
load("lam_grid_poly_10K.RData")

# Simulation
y1predA <- matrix(NA, nrow(mu1), niter)
y2predA <- matrix(NA, nrow(mu2), niter)
date()
for (i in 1:niter){
  if(i%%1000==0) {
    cat(paste0("iteration: ", i, "\n"))
  }
  log_lam1 <- log(lam_grid_poly_10K[cbind(round(mu1[, i]*1000), 
                                          round(nu1big[, i]*1000))])
  log_lam1_mat <- matrix(log_lam1, 11, n, T)
  unique <- 0:10
  lfac <- lfactorial(unique)
  G_mat1 <- matrix(rowSums(exp(tcrossprod(log_lam1, unique) - 
                                 tcrossprod(nu1big[, i], lfac))), 11, n, T)
  x <- 0:10
  lfac_mat <- matrix(lfactorial(0:10), 11, 1)
  pmf1 <- exp(x*log_lam1_mat - lfac_mat%*%nu1big[, i] - log(G_mat1))
  # Repeat for y2
  log_lam2 <- log(lam_grid_poly_10K[cbind(round(mu2[, i]*1000), 
                                          round(nu2big[, i]*1000))])
  log_lam2_mat <- matrix(log_lam2, 11, n, T)
  G_mat2 <- matrix(rowSums(exp(tcrossprod(log_lam2, unique) - 
                                 tcrossprod(nu2big[, i], lfac))), 11, n, T)
  pmf2 <- exp(x*log_lam2_mat - lfac_mat%*%nu2big[, i] - log(G_mat2))
  cdf_y1 <- apply(pmf1, 2, cumsum)
  cdf_y2 <- apply(pmf2, 2, cumsum)
  y1predA[ ,i] <- rowSums(t(cdf_y1) < u1[, i])
  y2predA[ ,i] <- rowSums(t(cdf_y2) < u2[, i])
  # Now inflate specific pairs
  pi_00 <- rbinom(n, 1, pi[i, 1]) == 1
  y1predA[pi_00, i] <- 0
  y2predA[pi_00, i] <- 0
  pi_02 <- rbinom(n, 1, pi[i, 2]) == 1
  y1predA[pi_02, i] <- 0
  y2predA[pi_02, i] <- 2
}
