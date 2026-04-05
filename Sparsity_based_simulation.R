set.seed(1)

#---------------------------
# Global simulation settings 
#---------------------------

# Dimensions considered in the comparison 
p_list <- c(5, 50, 100)

# Signal-strength grid used for each dimension 
r_grid_fun <- function(p) seq(0, 2.5 * sqrt(p), length.out = 6)

# Sparsity levels considered for each dimension 
k_grid_fun <- function(p) {
  raw_k <- c(1, 2, 5, 10, round(0.1 * p), round(0.2 * p), round(0.5 * p), p)
  sort(unique(pmax(1L, pmin(p, raw_k))))
}

# Monte Carlo sizes and MCMC settings 
Nmc_radial <- 5000
Nmc_hs    <- 500
n_iter_hs <- 3000
burn_hs   <- 1000
thin_hs   <- 2

#--------------
# Basic helpers
#--------------

# Sample from an inverse-gamma distribution via reciprocal gamma sampling
rinv_gamma <- function(n, shape, rate) {
  1 / rgamma(n, shape = shape, rate = rate)
}

# Generate Nmc observations from N_p(theta, I_p)
draw_Y <- function(theta, Nmc) {
  p <- length(theta)
  Y <- matrix(rnorm(Nmc * p), nrow = Nmc, ncol = p)
  sweep(Y, 2, theta, "+")
}

# Construct a k-sparse mean vector with Euclidean norm r
theta_k_sparse <- function(r, p, k) {
  stopifnot(k >= 1, k <= p)
  c(rep(r / sqrt(k), k), rep(0, p - k))
}

#--------------------------
# BetaPrime shrinkage rules
#--------------------------

# Radial shrinkage factor for the BetaPrime minimax rule as a function of s = ||y||^2
shrink_betaprime_minimax <- function(s, p) {
  beta <- p / 2 - 2
  if (beta <= 0) stop("Need p >= 5 so that beta = p/2 - 2 > 0.")
  alpha <- beta + p / 2
  lambda <- s / 2
  out <- numeric(length(s))
  
  # Limit at s = 0 
  i0 <- which(lambda == 0)
  if (length(i0) > 0) {
    out[i0] <- 1 - alpha / (alpha + 1)
  }
  
  # General case s > 0
  i1 <- which(lambda > 0)
  if (length(i1) > 0) {
    lam <- lambda[i1]
    P_a   <- pgamma(lam, shape = alpha,     scale = 1, lower.tail = TRUE)
    P_ap1 <- pgamma(lam, shape = alpha + 1, scale = 1, lower.tail = TRUE)
    Et <- alpha * P_ap1 / (lam * pmax(P_a, .Machine$double.eps))
    out[i1] <- 1 - Et
  }
  pmin(pmax(out, 0), 1)
}

#--------------------------------------
# Risk estimation for radial estimators 
#--------------------------------------

# Compute empirical risk from a matrix of observations for a radial shrinkage rule
risk_from_Y_radial <- function(Y, theta, shrink_fun) {
  s <- rowSums(Y^2)
  a <- shrink_fun(s)
  delta <- Y * a
  err <- sweep(delta, 2, theta, "-")
  mean(rowSums(err^2))
}

# Monte Carlo risk estimate at a fixed theta for a radial estimator 
risk_at_theta_radial <- function(theta, Nmc, shrink_fun) {
  Y <- draw_Y(theta, Nmc)
  risk_from_Y_radial(Y, theta, shrink_fun)
}

#--------------------------------------------
# Horseshoe posterior mean via Gibbs sampling 
#--------------------------------------------

# Compute the posterior mean under the horseshoe prior using a Gibbs sampler 
# with Rao-Blackwellization and latent inverse-gamma updates 

hs_posterior_mean <- function(y,
                              n_iter = 4000,
                              burn   = 1000,
                              thin   = 2,
                              init_tau2 = 1,
                              init_lambda2 = NULL) {
  p <- length(y)
  if (is.null(init_lambda2)) init_lambda2 <- rep(1, p)
  beta <- rep(0, p)
  lambda2 <- init_lambda2
  nu <- rep(1, p)
  tau2 <- init_tau2
  xi <- 1
  keep <- 0L
  post_mean_sum <- rep(0, p)
  for (iter in 1:n_iter) {
    
    # Gaussian full conditional for beta 
    post_var_beta  <- (tau2 * lambda2) / (1 + tau2 * lambda2)
    post_mean_beta <- post_var_beta * y
    beta <- rnorm(p, mean = post_mean_beta, sd = sqrt(post_var_beta))
    
    # Local scale updates 
    rate_lambda <- 1 / nu + beta^2 / (2 * tau2)
    lambda2 <- rinv_gamma(p, shape = 1, rate = rate_lambda)
    nu <- rinv_gamma(p, shape = 1, rate = 1 + 1 / lambda2)
    
    # Global scale updates
    rate_tau <- 1 / xi + sum(beta^2 / (2 * lambda2))
    tau2 <- rinv_gamma(1, shape = (p + 1) / 2, rate = rate_tau)
    xi <- rinv_gamma(1, shape = 1, rate = 1 + 1 / tau2)
    
    # Rao-Blackwellized posterior mean estimate
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      keep <- keep + 1L
      rb_mean <- (tau2 * lambda2 / (1 + tau2 * lambda2)) * y
      post_mean_sum <- post_mean_sum + rb_mean
    }
  }
  if (keep == 0L) {
    stop("No posterior draws kept. Increase n_iter or reduce burn/thin.")
  }
  list(
    mean    = post_mean_sum / keep,
    tau2    = tau2,
    lambda2 = lambda2
  )
}

#-----------------------------------------------
# Risk estimate for the horseshoe posterior mean 
#-----------------------------------------------

# Compute empirical risk from a matrix of observations by fitting the horseshoe 
# posterior mean separately to each row of Y 
risk_from_Y_hs <- function(Y, theta,
                           n_iter = 4000,
                           burn   = 1000,
                           thin   = 2) {
  Nmc <- nrow(Y)
  p <- length(theta)
  errs <- numeric(Nmc)
  
  # Warm starts are used to reduce computation across repeated fits
  last_tau2 <- 1
  last_lambda2 <- rep(1, p)
  for (m in 1:Nmc) {
    fit <- hs_posterior_mean(
      y = Y[m, ],
      n_iter = n_iter,
      burn = burn,
      thin = thin,
      init_tau2 = last_tau2,
      init_lambda2 = last_lambda2
    )
    delta_y <- fit$mean
    last_tau2 <- fit$tau2
    last_lambda2 <- fit$lambda2
    errs[m] <- sum((delta_y - theta)^2)
  }
  mean(errs)
}

# Monte Carlo risk estimate at a fixed theta for the horseshoe estimator
risk_at_theta_hs <- function(theta,
                             Nmc    = 200,
                             n_iter = 4000,
                             burn   = 1000,
                             thin   = 2) {
  Y <- draw_Y(theta, Nmc)
  risk_from_Y_hs(Y, theta, n_iter = n_iter, burn = burn, thin = thin)
}
#---------------------------------
# Main simulator for one dimension 
#---------------------------------
# For a given p: 
# 1. compute the BetaPrime risk curve (common across sparsity levels), 
# 2. compute horseshoe risk curves for each sparsity level, 
# 3. return all quantities needed for plotting 

run_one_dimension <- function(p,
                              Nmc_radial = 5000,
                              Nmc_hs = 25,
                              n_iter_hs = 1200,
                              burn_hs = 400,
                              thin_hs = 2) {
  cat("\n========================================\n")
  cat("Running comparison for p =", p, "\n")
  cat("========================================\n")
  r_grid <- r_grid_fun(p)
  k_grid <- k_grid_fun(p)
  shrink_bp <- function(s) shrink_betaprime_minimax(s, p)
  
  # Risk of the MLE Y under squared error loss 
  risk_mle <- rep(p, length(r_grid))
  
  # BetaPrime rule is radial, so one curve suffices for all sparsity levels
  risk_bp <- numeric(length(r_grid))
  cat("Computing BetaPrime risk curve (common to all k) ...\n")
  for (i in seq_along(r_grid)) {
    r <- r_grid[i]
    theta_ref <- theta_k_sparse(r, p, k = 1)  
    Y_radial <- draw_Y(theta_ref, Nmc_radial)
    risk_bp[i] <- risk_from_Y_radial(Y_radial, theta_ref, shrink_bp)
    cat("  r =", round(r, 3), "done (BetaPrime)\n")
  }
  
  # Horseshoe risk is evaluated separately for each sparsity level 
  risk_hs <- vector("list", length(k_grid))
  names(risk_hs) <- paste0("k=", k_grid)
  for (j in seq_along(k_grid)) {
    k <- k_grid[j]
    cat("Computing horseshoe curve for k =", k, "...\n")
    curve_j <- numeric(length(r_grid))
    for (i in seq_along(r_grid)) {
      r <- r_grid[i]
      theta <- theta_k_sparse(r, p, k)
      curve_j[i] <- risk_at_theta_hs(
        theta,
        Nmc    = Nmc_hs,
        n_iter = n_iter_hs,
        burn   = burn_hs,
        thin   = thin_hs
      )
      cat("  k =", k, ", r =", round(r, 3), "done (HS)\n")
    }
    risk_hs[[j]] <- curve_j
  }
  list(
    p = p,
    r_grid = r_grid,
    k_grid = k_grid,
    risk_mle = risk_mle,
    risk_bp = risk_bp,
    risk_hs = risk_hs
  )
}

#--------------------------------------
# Run the simulation for all dimensions 
#--------------------------------------

results <- vector("list", length(p_list))
names(results) <- paste0("p=", p_list)
for (idx in seq_along(p_list)) {
  p <- p_list[idx]
  results[[idx]] <- run_one_dimension(
    p = p,
    Nmc_radial = Nmc_radial,
    Nmc_hs = Nmc_hs,
    n_iter_hs = n_iter_hs,
    burn_hs = burn_hs,
    thin_hs = thin_hs
  )
}

#-------------
# Plot styling 
#-------------

curve_cols <- c(
  "black",    
  "#D55E00",   
  "#0072B2",   
  "#009E73"    
)
curve_lty <- c(
  1,  
  1,
  3   
)
curve_pch <- c(
  1,  
  16, 
  17  
)

# Add a curve with point markers for readability 
add_styled_curve <- function(x, y, col, lty, pch, mark_every = 1, lwd = 2.5) {
  lines(x, y, col = col, lty = lty, lwd = lwd)
  idx <- seq(1, length(x), by = mark_every)
  points(x[idx], y[idx], col = col, pch = pch, cex = 0.8)
}

# Generate a color palette for the collection of horseshoe curves 
hs_palette <- function(n) {
  if (n == 1) return(curve_cols[2])
  grDevices::colorRampPalette(c(curve_cols[2], curve_cols[4]))(n)
}

#--------------------------------
# Plot results for each dimension
#--------------------------------

plot_results <- function(results, open_new_device = FALSE) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  for (res in results) {
    if (open_new_device) grDevices::dev.new()
    hs_cols <- hs_palette(length(res$k_grid))
    all_vals <- c(res$risk_mle, res$risk_bp, unlist(res$risk_hs))
    plot(
      res$r_grid, res$risk_mle,
      type = "l",
      lwd = 2.5,
      lty = curve_lty[1],
      col = curve_cols[1],
      xlab = expression(paste("\u2016", theta, "\u2016")),
      ylab = "Estimated risk",
      ylim = range(all_vals),
      main = paste0("Risk comparison at p = ", res$p)
    )
    
    # MLE benchmark 
    idx_mle <- seq(1, length(res$r_grid), by = 1)
    points(
      res$r_grid[idx_mle], res$risk_mle[idx_mle],
      col = curve_cols[1],
      pch = curve_pch[1],
      cex = 0.8
    )
    
    # BetaPrime curve
    add_styled_curve(
      res$r_grid, res$risk_bp,
      col = curve_cols[3],
      lty = curve_lty[3],
      pch = curve_pch[3],
      mark_every = 1
    )
    
    # Horseshoe curves for each sparsity level
    for (j in seq_along(res$k_grid)) {
      add_styled_curve(
        res$r_grid, res$risk_hs[[j]],
        col = hs_cols[j],
        lty = curve_lty[2],
        pch = curve_pch[2],
        mark_every = 1
      )
    }
    legend(
      "topleft",
      legend = c(
        "MLE (= p)",
        "BNN with BetaPrime hyperprior",
        paste0("Horseshoe, k = ", res$k_grid)
      ),
      lty = c(curve_lty[1], curve_lty[3], rep(curve_lty[2], length(res$k_grid))),
      pch = c(curve_pch[1], curve_pch[3], rep(curve_pch[2], length(res$k_grid))),
      col = c(curve_cols[1], curve_cols[3], hs_cols),
      lwd = 2.5,
      pt.cex = 0.8,
      cex = 0.75,
      bty = "n"
    )
  }
}
# Produce the risk plots
plot_results(results)