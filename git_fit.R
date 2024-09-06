# R code to fit bivariate copula (0, 0) inflated MPCMP model 
fit_Big5_mean_centre <- function(dat, niters = 1000, thin = 1, 
                                 burnin = 1000, block_beta = TRUE){
  niters <- niters + burnin
  print("Getting initial values")
  # Fit initial model for starting values
  YC <- dat$YC
  Home_NoFans <- dat$Home * dat$NoFans
  # Indicators for (0, 0) and (0, 2) outcomes
  zzI <- dat$YC[dat$Home == 1] == 0 & dat$YC[dat$Home == 0] == 0
  zzI2 <- dat$YC[dat$Home == 1] == 0 & dat$YC[dat$Home == 0] == 2
  # Poisson fit assuming independence
  fit.original <- glm(YC ~ -1 + League + Home + NoFans + Home_NoFans, family = "poisson", data = dat, 
                      x = TRUE)
  beta.hat <- fit.original$coefficients[1:8]
  p <- length(beta.hat)
  # Referee effects
  nref <- nlevels(dat$Referee)
  ntheta <- nref
  theta.hat <- rep(0, nref)
  # No. of estimated referees in each league - used in sig_theta_sq update
  niref <- colSums(with(dat, table(Referee, League)) != 0)
  # Referee indices needed for update too
  refind <- rep(1:5, niref)
  # Team effects
  nteam <- nlevels(dat$Team)
  ngamma <- nteam
  gamma.hat <- rep(0, nteam)
  # No. of estimated teams in each league - used in sig_gamma_sq update
  niteam <- colSums(with(Big5_18_22_Jun_24, table(Team, League)) != 0)
  # Team indices needed for update too
  teamind <- rep(1:5, niteam)
  # Fitted dispersion - set as overall for each league to initialise
  nu.hat <- rep(1, 5)
  
  # Set-up two datasets and extract responses
  x.matrix <- fit.original$x
  dat$Home_NoFans <- dat$Home * dat$NoFans
  dat1 <- dat[dat$Home == 1, ]
  dat2 <- dat[dat$Home == 0, ]
  n <- nrow(dat1)
  nleague <- as.numeric(table(dat1$League))
  y1 <- dat1$YC
  y2 <- dat2$YC
  # Lagged counts
  dat1d <- dat1
  dat1d$YC <- dat1d$YC - 1
  y1d <- dat1d$YC
  dat2d <- dat2
  dat2d$YC <- dat2d$YC - 1
  y2d <- dat2d$YC 
  print("Design matrices found")
  Sigma_beta <- 0.2*fit.original$variance_beta[1:p, 1:p]
  Sigma_theta <- 0.02*fit.original$variance_beta[(p + 1):(p + 1 + ntheta - 1), 
                                                 (p + 1):(p + 1 + ntheta - 1)]
  Sigma_gamma <- 0.20*fit.original$variance_beta[(p + 1 + ntheta):(p + 1 + ntheta + ngamma - 1), 
                                                 (p + 1 + ntheta):(p + 1 + ntheta + ngamma - 1)]
  
  # Random walk proposal values
  rw_nu <- 0.10
  rw_kappa <- 0.35
  
  # Prior means and SDs
  beta_m <- 0
  beta_s <- 0.25
  theta_m <- 0
  gamma_m <- 0
  nu_m <- 0 # Equidispersion (on log-scale)
  nu_s <- 1
  # Bivariate specific paras
  # Normal prior for kappa
  kappa_m <- 0.5
  kappa_s <- 1
  # Hyperparameters for sigma_theta_sq
  a <- 3
  b <- 0.2
  # Beta hyperparameters for inflation terms
  pi_za <- 1
  pi_zb <- 49
  
  # Initial values 
  beta0 <- beta.hat
  theta0 <- theta.hat
  gamma0 <- gamma.hat
  # Initial values for nu (log-scale) for each league
  nu10 <- nu20 <- log(nu.hat)
  # Initialise linear predictors and set up design matrices
  x.matrix11 <- x.matrix[dat$Home == 1, 1:p] 
  mu10 <- with(dat1, 
               exp(x.matrix11 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
  x.matrix21 <- x.matrix[dat$Home == 0, 1:p] 
  x.matrix22 <- Refs[dat$Home == 0, ]  
  x.matrix23 <- Teams[dat$Home == 0, ] 
  mu20 <- with(dat2, 
               exp(x.matrix21 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
  if (max(mu10) >= 10){mu10 <- pmin(9.9, mu10)}
  if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
  if (max(mu20) >= 10){mu20 <- pmin(9.9, mu20)}
  if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
  # Initialise copula dependence parameters
  kappa0 <- rep(2, 5)
  # Initialise theta variance hyperparameters
  sig_theta_sq <- rinvgamma(5, a, b)
  sig_theta_sq_long <- rep(sig_theta_sq, niref)
  # Initialise gamma variance hyperparameters
  sig_gamma_sq <- rinvgamma(5, a, b)
  sig_gamma_sq_long <- rep(sig_gamma_sq, niteam)
  # Initialise inflation parameters
  pi_z0 <- c(rbeta(1, pi_za, pi_zb), rbeta(1, pi_za, pi_zb))
  
  K <- niters
  betas <- matrix(NA, nrow = K, ncol = p)
  thetas <- matrix(NA, nrow = K, ncol = ntheta)
  gammas <- matrix(NA, nrow = K, ncol = ngamma)
  nus <- matrix(NA, nrow = K, ncol = 10)
  kappa <- matrix(NA, nrow = K, ncol = 5)
  var_theta <- matrix(NA, nrow = K, ncol = 5) # One for each league
  var_gamma <- matrix(NA, nrow = K, ncol = 5) # One for each league
  pis <- matrix(NA, nrow = K, ncol = 2)
  ll <- matrix(NA, nrow = K, ncol = 1)
  C <- matrix(NA, nrow = K, ncol = n)
  colnames(betas) <-  paste("x", 0:(p-1), sep = "")
  colnames(thetas) <-  paste("ref", 1:length(theta.hat), sep = "")
  colnames(gammas) <-  paste("team", 1:length(gamma.hat), sep = "")
  colnames(nus) <-  c(paste("nu_1", 1:5, sep = ""),
                      paste("nu_2", 1:5, sep = ""))
  colnames(kappa) <- paste("kap", 1:5, sep = "")
  colnames(var_theta) <- paste("var_th", 1:5, sep = "")
  colnames(var_gamma) <- paste("var_ga", 1:5, sep = "")
  colnames(pis) <- c("pi_00", "pi_02")
  iter_times <- rep(0, K)
  
  # Data preparation
  # Summaries for home yellow data
  summ1 <- summarydata(dat1, col = "YC")
  lfac1 <- summ1$lfac
  lfac_y1 <- summ1$lfac_y
  lfac_y1_max <- lfactorial(0:max(y1))
  N <- summ1$N
  unique1 <- summ1$unique
  # Summaries for away yellow data
  summ2 <- summarydata(dat2, col = "YC")
  lfac2 <- summ2$lfac
  lfac_y2 <- summ2$lfac_y
  lfac_y2_max <- lfactorial(0:max(y2))
  unique2 <- summ2$unique
  # Summaries for lagged home yellow data
  summ1d <- summarydata(dat1d, col = "YC")
  lfac_y1d <- summ1d$lfac_y
  lfac_y1d_max <- lfactorial(0:max(y1d))
  # Summaries for lagged away yellow data
  summ2d <- summarydata(dat2d, col = "YC")
  lfac_y2d <- summ2d$lfac_y
  lfac_y2d_max <- lfactorial(0:max(y2d))
  
  # Set up copula indicators
  M1 <- matrix(y1, n, max(y1) + 1)
  M1max <- matrix(0:max(y1), n, max(y1) + 1, byrow = T)
  M1ind <- M1max <= M1
  M1ind_d <- M1max[, -(max(y1) + 1)] < M1[, -(max(y1) + 1)]
  M2 <- matrix(y2, n, max(y2) + 1)
  M2max <- matrix(0:max(y2), n, max(y2) + 1, byrow = T)
  M2ind <- M2max <= M2
  M2ind_d <- M2max[, -(max(y2) + 1)] < M2[, -(max(y2) + 1)]
  
  ptm <- proc.time()
  print("Starting MCMC algorithm")
  for(k in 1:K){
    #print(k)
    start <- proc.time()[3]
    if(k%%1000==0) {
      cat(paste0("iteration: ", k, "\n"))
    }
    #### beta update ####
    if (block_beta) {beta1 <- rmvnorm(1, beta0, Sigma_beta)[1, ]}
    else {beta1 <- rmvnorm(1, beta0, 
                           c(rep(0.001, 8))*diag(p))[1, ]}
    mu11 <- with(dat1, 
                 exp(x.matrix11 %*% beta1 + theta0[RefID] + gamma0[TeamID]))
    mu21 <- with(dat2, 
                 exp(x.matrix21 %*% beta1 + theta0[RefID] + gamma0[TeamID]))
    nu10_lp <- rep(exp(nu10), nleague)
    nu20_lp <- rep(exp(nu20), nleague)
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    if (max(nu10_lp >= 9.9)){nu10_lp <- pmin(9.9, nu10_lp)}
    if (min(nu10_lp) < 0.001){nu10_lp <- pmax(0.001, nu10_lp)}
    if (max(nu20_lp >= 9.9)){nu20_lp <- pmin(9.9, nu20_lp)}
    if (min(nu20_lp) < 0.001){nu20_lp <- pmax(0.001, nu20_lp)}
    # Find lambda
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000), 
                                               round(nu10_lp*1000))]) 
    log_lambda11 <- log(lam_grid_poly_10K[cbind(round(mu11*1000), 
                                                round(nu10_lp*1000))])
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000), 
                                               round(nu20_lp*1000))]) 
    log_lambda21 <- log(lam_grid_poly_10K[cbind(round(mu21*1000), 
                                                round(nu20_lp*1000))])
    # Calculate normalising constant
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) - 
                        tcrossprod(nu10_lp, lfac1)))
    G11 <- rowSums(exp(tcrossprod(log_lambda11, unique1) - 
                         tcrossprod(nu10_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) - 
                        tcrossprod(nu20_lp, lfac2))) 
    G21 <- rowSums(exp(tcrossprod(log_lambda21, unique2) - 
                         tcrossprod(nu20_lp, lfac2))) 
    # Some terms for copula calculation
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) - 
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) - 
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) - 
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) - 
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    Gy11 <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1)) - 
                          tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G11
    Gy21 <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2)) - 
                          tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G21
    Gy11d <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1d)) - 
                           tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G11
    Gy21d <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2d)) - 
                           tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G21
    # Define copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z0[1]
    C11[zzI2] <- C11[zzI2] + pi_z0[2]
    C21 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C31 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C41 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    x1c <- x.matrix11
    x2c <- x.matrix21
    denrat <- colSums(x1c*log(C11 - C21 - C31 + C41) + 
                        x2c*log(C11 - C21 - C31 + C41)) - 
      colSums(x1c*log(C1 - C2 - C3 + C4) + 
                x2c*log(C1 - C2 - C3 + C4))
    denrat <- denrat - 0.5*((beta1 - beta_m)^2 - (beta0 - beta_m)^2)/beta_s^2
    laccept <- pmin(0, denrat)
    accept <- (log(runif(p)) < laccept)
    beta0[accept] <- beta1[accept]
    mu10 <- with(dat1, 
                 exp(x.matrix11 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
    mu20 <- with(dat2, 
                 exp(x.matrix21 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    betas[k,] <- beta0
    
    #### kappa update ####
    kappa1 <- rnorm(5, kappa0, rw_kappa)
    # lambda & G terms don't change but need latest update
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                               round(nu10_lp*1000))])
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                               round(nu20_lp*1000))])
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) -
                        tcrossprod(nu10_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) -
                        tcrossprod(nu20_lp, lfac2)))
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) -
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) -
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) -
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) -
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    # Copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa1[LeagueID]*Gy1) - 1)*
                                                         (exp(-kappa1[LeagueID]*Gy2) - 1)/
                                                         (exp(-kappa1[LeagueID]) - 1))/kappa1[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z0[1]
    C11[zzI2] <- C11[zzI2] + pi_z0[2]
    C21 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa1[LeagueID]*Gy1d) - 1)*
                                                         (exp(-kappa1[LeagueID]*Gy2) - 1)/
                                                         (exp(-kappa1[LeagueID]) - 1))/kappa1[LeagueID])
    C31 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa1[LeagueID]*Gy1) - 1)*
                                                         (exp(-kappa1[LeagueID]*Gy2d) - 1)/
                                                         (exp(-kappa1[LeagueID]) - 1))/kappa1[LeagueID])
    C41 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa1[LeagueID]*Gy1d) - 1)*
                                                         (exp(-kappa1[LeagueID]*Gy2d) - 1)/
                                                         (exp(-kappa1[LeagueID]) - 1))/kappa1[LeagueID])
    # Either dat1 or dat2 below as league id fixed
    denrat <- with(dat1, tapply(log(C11 - C21 - C31 + C41)  -
                                  log(C1 - C2 - C3 + C4), LeagueID, sum))
    denrat <- denrat - 0.5*((kappa1 - kappa_m)^2 -
                              (kappa0 - kappa_m)^2)/kappa_s^2
    laccept <- pmin(0, denrat)
    accept <- (log(runif(5)) < laccept)
    kappa0[accept] <- kappa1[accept]
    #kappa0 <- rep(0.001, 5)
    kappa[k, ] <- kappa0
    
    #### theta update (referees) ####
    theta1  <- rmvnorm(1, theta0, 0.10*diag(ntheta))[1, ]
    mu11 <- with(dat1, 
                 exp(x.matrix11 %*% beta0 + theta1[RefID] + gamma0[TeamID]))
    mu21 <- with(dat2, 
                 exp(x.matrix21 %*% beta0 + theta1[RefID] + gamma0[TeamID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    # Find lambda for proposed values
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                               round(nu10_lp*1000))])
    log_lambda11 <- log(lam_grid_poly_10K[cbind(round(mu11*1000),
                                                round(nu10_lp*1000))])
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                               round(nu20_lp*1000))])
    log_lambda21 <- log(lam_grid_poly_10K[cbind(round(mu21*1000),
                                                round(nu20_lp*1000))])
    # Calculate normalising constants
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) -
                        tcrossprod(nu10_lp, lfac1)))
    G11 <- rowSums(exp(tcrossprod(log_lambda11, unique1) -
                         tcrossprod(nu10_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) -
                        tcrossprod(nu20_lp, lfac2)))
    G21 <- rowSums(exp(tcrossprod(log_lambda21, unique2) -
                         tcrossprod(nu20_lp, lfac2)))
    # Some terms for copula calculation
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) -
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) -
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) -
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) -
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    Gy11 <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1)) -
                          tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G11
    Gy21 <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2)) -
                          tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G21
    Gy11d <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1d)) -
                           tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G11
    Gy21d <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2d)) -
                           tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G21
    # Define copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z0[1]
    C11[zzI2] <- C11[zzI2] + pi_z0[2]
    C21 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C31 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C41 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    # Doesn't matter to sum using dat1 or dat2 as ref ids are fixed
    denrat <- with(dat1, tapply(log(C11 - C21 - C31 + C41)  -
                                  log(C1 - C2 - C3 + C4), RefID, sum))
    denrat <- denrat - 0.5*((theta1 - theta_m)^2 - (theta0 - theta_m)^2)/sig_theta_sq_long
    laccept <- pmin(0, denrat)
    accept <- (log(runif(ntheta)) < laccept)
    theta0[accept] <- theta1[accept]
    # MEAN-CENTERING HERE
    theta0[1:niref[1]] <- theta0[1:niref[1]] - mean(theta0[1:niref[1]])
    theta0[(niref[1]+1):niref[2]] <- theta0[(niref[1]+1):niref[2]] - 
      mean(theta0[(niref[1]+1):niref[2]])
    theta0[(niref[2]+1):niref[3]] <- theta0[(niref[2]+1):niref[3]] - 
      mean(theta0[(niref[2]+1):niref[3]])
    theta0[(niref[3]+1):niref[4]] <- theta0[(niref[3]+1):niref[4]] - 
      mean(theta0[(niref[3]+1):niref[4]])
    theta0[(niref[4]+1):niref[5]] <- theta0[(niref[4]+1):niref[5]] - 
      mean(theta0[(niref[4]+1):niref[5]])
    # Update linear predictors
    mu10 <- with(dat1, 
                 exp(x.matrix11 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
    mu20 <- with(dat2, 
                 exp(x.matrix21 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    # Store update
    thetas[k,] <- theta0
    
    #### gamma update (teams) ####
    gamma1  <- rmvnorm(1, gamma0, 0.02*diag(ngamma))[1, ]
    mu11 <- with(dat1, 
                 exp(x.matrix11 %*% beta0 + theta0[RefID] + gamma1[TeamID]))
    mu21 <- with(dat2, 
                 exp(x.matrix21 %*% beta0 + theta0[RefID] + gamma1[TeamID]))
    if (max(mu11) >= 9.9){mu11 <- pmin(9.8, mu11)}
    if (min(mu11) < 0.001){mu11 <- pmax(0.001, mu11)}
    if (max(mu21) >= 9.9){mu21 <- pmin(9.81, mu21)}
    if (min(mu21) < 0.001){mu21 <- pmax(0.001, mu21)}
    # Find lambda for proposed values
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                               round(nu10_lp*1000))])
    log_lambda11 <- log(lam_grid_poly_10K[cbind(round(mu11*1000),
                                                round(nu10_lp*1000))])
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                               round(nu20_lp*1000))])
    log_lambda21 <- log(lam_grid_poly_10K[cbind(round(mu21*1000),
                                                round(nu20_lp*1000))])
    # Calculate normalising constants
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) -
                        tcrossprod(nu10_lp, lfac1)))
    G11 <- rowSums(exp(tcrossprod(log_lambda11, unique1) -
                         tcrossprod(nu10_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) -
                        tcrossprod(nu20_lp, lfac2)))
    G21 <- rowSums(exp(tcrossprod(log_lambda21, unique2) -
                         tcrossprod(nu20_lp, lfac2)))
    # Some terms for copula calculation
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) -
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) -
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) -
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) -
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    Gy11 <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1)) -
                          tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G11
    Gy21 <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2)) -
                          tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G21
    Gy11d <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1d)) -
                           tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G11
    Gy21d <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2d)) -
                           tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G21
    # Define copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z0[1]
    C11[zzI2] <- C11[zzI2] + pi_z0[2]
    C21 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C31 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C41 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    denrat <- with(dat1, tapply(log(C11 - C21 - C31 + C41)  -
                                  log(C1 - C2 - C3 + C4), TeamID, sum))
    denrat <- denrat + with(dat2, tapply(log(C11 - C21 - C31 + C41)  -
                                           log(C1 - C2 - C3 + C4), TeamID, sum))
    denrat <- denrat - 0.5*((gamma1 - gamma_m)^2 - (gamma0 - gamma_m)^2)/sig_gamma_sq_long
    laccept <- pmin(0, denrat)
    accept <- (log(runif(ngamma)) < laccept)
    gamma0[accept] <- gamma1[accept]
    # MEAN-CENTERING HERE
    gamma0[1:niteam[1]] <- gamma0[1:niteam[1]] - mean(gamma0[1:niteam[1]])
    gamma0[(niteam[1]+1):niteam[2]] <- gamma0[(niteam[1]+1):niteam[2]] - 
      mean(gamma0[(niteam[1]+1):niteam[2]])
    gamma0[(niteam[2]+1):niteam[3]] <- gamma0[(niteam[2]+1):niteam[3]] - 
      mean(gamma0[(niteam[2]+1):niteam[3]])
    gamma0[(niteam[3]+1):niteam[4]] <- gamma0[(niteam[3]+1):niteam[4]] - 
      mean(gamma0[(niteam[3]+1):niteam[4]])
    gamma0[(niteam[4]+1):niteam[5]] <- gamma0[(niteam[4]+1):niteam[5]] - 
      mean(gamma0[(niteam[4]+1):niteam[5]])
    # Update linear predictors
    mu10 <- with(dat1, 
                 exp(x.matrix11 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
    mu20 <- with(dat2, 
                 exp(x.matrix21 %*% beta0 + theta0[RefID] + gamma0[TeamID]))
    if (max(mu10) >= 9.9){mu10 <- pmin(9.8, mu10)}
    if (min(mu10) < 0.001){mu10 <- pmax(0.001, mu10)}
    if (max(mu20) >= 9.9){mu20 <- pmin(9.81, mu20)}
    if (min(mu20) < 0.001){mu20 <- pmax(0.001, mu20)}
    # Store update
    gammas[k,] <- gamma0
    
    #### nu update (one for each league and H/A within league) ####
    # First nu1
    nu11 <- rnorm(5, nu10, rw_nu)
    nu11_lp <- rep(exp(nu11), nleague)
    if (max(nu11_lp) >= 9.9){nu11_lp <- pmin(9.9, nu11_lp)}
    if (min(nu11_lp) < 0.001){nu11_lp <- pmax(0.001, nu11_lp)}
    nu21_lp <- nu20_lp
    # Find lambda
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                               round(nu10_lp*1000))])
    log_lambda11 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                                round(nu11_lp*1000))])
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                               round(nu20_lp*1000))])
    log_lambda21 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                                round(nu21_lp*1000))])
    # Calculate normalising constant
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) -
                        tcrossprod(nu10_lp, lfac1)))
    G11 <- rowSums(exp(tcrossprod(log_lambda11, unique1) -
                         tcrossprod(nu11_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) -
                        tcrossprod(nu20_lp, lfac2)))
    G21 <- rowSums(exp(tcrossprod(log_lambda21, unique2) -
                         tcrossprod(nu21_lp, lfac2)))
    # Some terms for copula calculation
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) -
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) -
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) -
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) -
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    Gy11 <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1)) -
                          tcrossprod(nu11_lp, lfac_y1_max))*M1ind)/G11
    Gy21 <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2)) -
                          tcrossprod(nu21_lp, lfac_y2_max))*M2ind)/G21
    Gy11d <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1d)) -
                           tcrossprod(nu11_lp, lfac_y1d_max))*M1ind_d)/G11
    Gy21d <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2d)) -
                           tcrossprod(nu21_lp, lfac_y2d_max))*M2ind_d)/G21
    # Define copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z0[1]
    C11[zzI2] <- C11[zzI2] + pi_z0[2]
    C21 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C31 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C41 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    denrat <- with(dat1, tapply(log(C11 - C21 - C31 + C41)  -
                                  log(C1 - C2 - C3 + C4), League, sum))
    denrat <- denrat - 0.5*((nu11 - nu_m)^2 - (nu10 - nu_m)^2)/nu_s^2
    laccept <- pmin(0, denrat)
    accept <- (log(runif(5)) < laccept)
    nu10[accept] <- nu11[accept]
    # Now nu2
    nu21 <- rnorm(5, nu20, rw_nu)
    nu21_lp <- rep(exp(nu21), nleague)
    nu11_lp <- nu10_lp
    if (max(nu21_lp) >= 9.9){nu21_lp <- pmin(9.9, nu21_lp)}
    if (min(nu21_lp) < 0.001){nu21_lp <- pmax(0.001, nu21_lp)}
    # Find lambda
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                               round(nu10_lp*1000))])
    log_lambda11 <- log(lam_grid_poly_10K[cbind(round(mu10*1000),
                                                round(nu11_lp*1000))])
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                               round(nu20_lp*1000))])
    log_lambda21 <- log(lam_grid_poly_10K[cbind(round(mu20*1000),
                                                round(nu21_lp*1000))])
    # Calculate normalising constant
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) -
                        tcrossprod(nu10_lp, lfac1)))
    G11 <- rowSums(exp(tcrossprod(log_lambda11, unique1) -
                         tcrossprod(nu11_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) -
                        tcrossprod(nu20_lp, lfac2)))
    G21 <- rowSums(exp(tcrossprod(log_lambda21, unique2) -
                         tcrossprod(nu21_lp, lfac2)))
    # Some terms for copula calculation
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) -
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) -
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) -
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) -
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    Gy11 <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1)) -
                          tcrossprod(nu11_lp, lfac_y1_max))*M1ind)/G11
    Gy21 <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2)) -
                          tcrossprod(nu21_lp, lfac_y2_max))*M2ind)/G21
    Gy11d <- rowSums(exp(tcrossprod(log_lambda11, 0:max(y1d)) -
                           tcrossprod(nu11_lp, lfac_y1d_max))*M1ind_d)/G11
    Gy21d <- rowSums(exp(tcrossprod(log_lambda21, 0:max(y2d)) -
                           tcrossprod(nu21_lp, lfac_y2d_max))*M2ind_d)/G21
    # Define copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z0[1]
    C11[zzI2] <- C11[zzI2] + pi_z0[2]
    C21 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C31 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C41 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy11d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy21d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    denrat <- with(dat1, tapply(log(C11 - C21 - C31 + C41)  -
                                  log(C1 - C2 - C3 + C4), League, sum))
    denrat <- denrat - 0.5*((nu21 - nu_m)^2 - (nu20 - nu_m)^2)/nu_s^2
    laccept <- pmin(0, denrat)
    accept <- (log(runif(5)) < laccept)
    nu20[accept] <- nu21[accept]
    # Store update
    nus[k,] <- c(nu10, nu20)
    
    #### Referee variance update ####
    sig_theta_sq <- rinvgamma(5, a + niref, 
                              b + 0.5*tapply(theta0^2, refind, sum))
    #print(sig_theta_sq)
    sig_theta_sq_long <- rep(sig_theta_sq, niref)
    var_theta[k, ] <- sig_theta_sq
    
    #### Team variance update ####
    sig_gamma_sq <- rinvgamma(5, a + niteam,
                              b + 0.5*tapply(gamma0^2, teamind, sum))
    sig_gamma_sq_long <- rep(sig_gamma_sq, niteam)
    var_gamma[k, ] <- sig_gamma_sq
    
    #### pi_z updates for (0, 0) and (0, 2) inflation ####
    # Update on logit scale
    pi_z1_logit <- rnorm(2, log(pi_z0/(1 - pi_z0)), 0.45)
    pi_z1 <- exp(pi_z1_logit)/(1 + exp(pi_z1_logit))
    # lambda & G terms don't change but need latest update
    nu10_lp <- rep(exp(nu10), nleague)
    nu20_lp <- rep(exp(nu20), nleague)
    log_lambda1 <- log(lam_grid_poly_10K[cbind(round(mu10*1000), 
                                               round(nu10_lp*1000))]) 
    log_lambda2 <- log(lam_grid_poly_10K[cbind(round(mu20*1000), 
                                               round(nu20_lp*1000))]) 
    G1 <- rowSums(exp(tcrossprod(log_lambda1, unique1) - 
                        tcrossprod(nu10_lp, lfac1)))
    G2 <- rowSums(exp(tcrossprod(log_lambda2, unique2) - 
                        tcrossprod(nu20_lp, lfac2))) 
    Gy1 <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1)) - 
                         tcrossprod(nu10_lp, lfac_y1_max))*M1ind)/G1
    Gy2 <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2)) - 
                         tcrossprod(nu20_lp, lfac_y2_max))*M2ind)/G2
    Gy1d <- rowSums(exp(tcrossprod(log_lambda1, 0:max(y1d)) - 
                          tcrossprod(nu10_lp, lfac_y1d_max))*M1ind_d)/G1
    Gy2d <- rowSums(exp(tcrossprod(log_lambda2, 0:max(y2d)) - 
                          tcrossprod(nu20_lp, lfac_y2d_max))*M2ind_d)/G2
    # Copula terms
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11 <- (1 - pi_z1[1] - pi_z1[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C11[zzI] <- C11[zzI] + pi_z1[1]
    C11[zzI2] <- C11[zzI2] + pi_z1[2]
    C21 <- (1 - pi_z1[1] - pi_z1[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C31 <- (1 - pi_z1[1] - pi_z1[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C41 <- (1 - pi_z1[1] - pi_z1[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                         (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                         (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    denrat <- sum(log(C11 - C21 - C31 + C41))  - 
      sum(log(C1 - C2 - C3 + C4))
    denrat <- denrat + (pi_za - 1)*log(pi_z1/pi_z0) + 
      (pi_zb - 1)*log((1 - pi_z1)/(1 - pi_z0)) + 
      log(1 - pi_z1) - log(1 - pi_z0) + log(pi_z1/pi_z0)
    laccept <- pmin(0, denrat)
    accept <- (log(runif(2)) < laccept)
    pi_z0[accept] <- pi_z1[accept]
    pis[k, ] <- pi_z0
    
    #### Log-likelihood ####
    # Need all of the G & copula terms for this too
    # G terms same as in pi update as doesn't depend on pi
    C1 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C1[zzI] <- C1[zzI] + pi_z0[1]
    C1[zzI2] <- C1[zzI2] + pi_z0[2]
    C2 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C3 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    C4 <- (1 - pi_z0[1] - pi_z0[2])*with(dat1, -log1p((exp(-kappa0[LeagueID]*Gy1d) - 1)*
                                                        (exp(-kappa0[LeagueID]*Gy2d) - 1)/
                                                        (exp(-kappa0[LeagueID]) - 1))/kappa0[LeagueID])
    
    ll[k] <- sum(log(C1 - C2 - C3 + C4))
    # Store terms for WAIC calculation
    C[k, ] <- C1 - C2 - C3 + C4
    
    # Monitor iteration times
    iter_times[k] <- proc.time()[3] - start
  }
  paras <- cbind(betas, thetas, gammas, nus, kappa, 
                 var_theta, var_gamma, pis)
  if (burnin > 1){
    paras <- paras[-(1:burnin), ]
  }
  if (thin > 1){
    paras <- paras[seq(1, niters - burnin, thin), ]
  }
  cpu_time <- proc.time() - ptm
  acc_rates <- rep(NA, ncol(paras))
  if (niters > 3){
    acc_rates <- 100*colSums(apply(paras, 2, diff) != 0)/(K - 1)
  }
  list(cpu_time = cpu_time[3], paras = paras, iter_times = iter_times, 
       acc_rates = acc_rates, X11 = x.matrix11, X12 = x.matrix12,
       X13 = x.matrix13, X21 = x.matrix21, X22 = x.matrix22,
       X23 = x.matrix23, loglike = ll, Cmat = C)
}