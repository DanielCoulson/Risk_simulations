set.seed(1)
#--------------------------------
# Basic numerical helpers
#--------------------------------

# Stable log-sum-exp for normalizing log-weights
logsumexp <- function(z) {
  m <- max(z)
  m + log(sum(exp(z - m)))
}

# 1D interpolation used to evaluate shrinkage on a precomputed grid
interp1 <- function(xgrid, ygrid, x) {
  approx(xgrid, ygrid, x, rule = 2)$y
}
#----------------------------------------------------------
# Sampling the latent mixture index k for each hidden layer
#----------------------------------------------------------

# Sample k from the binomial-type mixture on {1,...,n}
rk_layer <- function(n, m = 1L) {
  ks <- 1:n
  w <- choose(n, ks)
  w <- w / sum(w)   
  sample(ks, size = m, replace = TRUE, prob = w)
}

# Vectorized version of rk_layer() when the layer widths vary across samples
rk_layer_vec <- function(n_vec) {
  out <- integer(length(n_vec))
  if (!length(n_vec)) return(out)
  n_pos <- n_vec[n_vec > 0]
  if (!length(n_pos)) return(out)
  for (n in sort(unique(n_pos))) {
    idx <- which(n_vec == n)
    ks <- 1:n
    w <- choose(n, ks) 
    w <- w / sum(w)
    out[idx] <- sample(ks, size = length(idx), replace = TRUE, prob = w)
  }
  out
}

#----------------------------------------------
# Sampling from the latent scale distribution V 
#----------------------------------------------

# Draw Monte Carlo samples from the fixed-scale BNN mixing distribution
rV_bnn_fixed <- function(M, d, widths, sigmas) {
  stopifnot(length(sigmas) == d)
  stopifnot(length(widths) == d - 1)
  const <- (2^(d - 1)) * prod(sigmas^2)
  Tprod <- rep(1, M)
  for (ell in 1:(d - 1)) {
    k <- rk_layer(widths[ell], m = M)
    T <- rgamma(M, shape = k / 2, rate = 1)
    Tprod <- Tprod * T
  }
  const * Tprod
}

# Draw Monte Carlo samples from the dropout BNN mixing distribution 
rV_bnn_dropout <- function(M, d, widths, sigmas, keep_probs, inverted = TRUE) {
  stopifnot(length(sigmas) == d)
  stopifnot(length(widths) == d - 1)
  if (length(keep_probs) == 1L) keep_probs <- rep(keep_probs, d - 1)
  stopifnot(length(keep_probs) == d - 1)
  if (any(keep_probs <= 0 | keep_probs > 1)) {
    stop("keep_probs must lie in (0,1].")
  }
  const <- (2^(d - 1)) * prod(sigmas^2)
  if (inverted) {
    const <- const * prod(1 / keep_probs)
  }
  Tprod <- rep(1, M)
  dead <- rep(FALSE, M)
  for (ell in 1:(d - 1)) {
    n_active <- rbinom(M, size = widths[ell], prob = keep_probs[ell])
    dead <- dead | (n_active == 0L)
    alive <- !dead
    if (any(alive)) {
      k <- rk_layer_vec(n_active[alive])
      T <- rgamma(sum(alive), shape = k / 2, rate = 1)
      Tprod[alive] <- Tprod[alive] * T
    }
  }
  V <- const * Tprod
  V[dead] <- 0
  V
}

#---------------------------------------------------------
# Closed form shrinkage for the Beta-prime hyperprior rule
#---------------------------------------------------------

# Evaluate the BetaPrime minimax shrinkage factor as a function of s = ||y||^2
shrink_betaprime_minimax <- function(s, p) {
  beta <- p / 2 - 2
  if (beta <= 0) stop("Need p >= 5 so that beta = p/2 - 2 > 0.")
  alpha <- beta + p / 2
  lambda <- s / 2
  out <- numeric(length(s))
  i0 <- which(lambda == 0)
  if (length(i0) > 0) out[i0] <- 1 - alpha / (alpha + 1)
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

#--------------------------------------------------------
# Precompute radial shrinkage from Monte Carlo draws of V 
#--------------------------------------------------------

# For each s on grid, approximate E[V/(V+1)| ||Y||^2 = s]
precompute_shrink_from_V <- function(Vdraw, s_grid, p) {
  a <- numeric(length(s_grid))
  logV1 <- log(Vdraw + 1)
  for (j in seq_along(s_grid)) {
    s <- s_grid[j]
    logw <- -(p / 2) * logV1 - s / (2 * (Vdraw + 1))
    logw <- logw - logsumexp(logw)
    w <- exp(logw)
    a[j] <- sum(w * (Vdraw / (Vdraw + 1)))
  }
  pmin(pmax(a, 0), 1)
}

#------------------------------------------------------------
# Monte Carlo risk estimation for radial shrinkage estimators 
#------------------------------------------------------------

#Estimate risk at a fixed scale norm r by randomizing the signal direction
risk_at_norm_shrink <- function(r, p, Nmc, shrink_fun) {
  u <- rnorm(p)
  u <- u / sqrt(sum(u^2))
  theta <- r * u
  Y <- matrix(rnorm(Nmc * p), nrow = Nmc, ncol = p)
  Y <- sweep(Y, 2, theta, "+")
  s <- rowSums(Y^2)
  a <- shrink_fun(s)
  delta <- Y * a
  mean(rowSums((delta - matrix(theta, nrow = Nmc, ncol = p, byrow = TRUE))^2))
}

#--------------------
# Simulation settings
#--------------------

p_vals <- c(5, 50, 100)
d <- 3
widths <- c(20, 20)
sigmas <- c(1.0, 1.0, 1.0)
keep_probs <- c(0.8, 0.8)

#Signal strength grid for the risk curves 
r_grid <- seq(0, 500, by = 1)

# Number of Monte Carlo draws used to approximate the V mixture distribution 
M_V <- 200000

# Dimension dependent Monte Carlo settings for risk estimation
mc_settings_for_p <- function(p) {
  if (p <= 5) {
    list(Nmc = 50000, Kdir = 10)
  } else if (p <= 50) {
    list(Nmc = 15000, Kdir = 6)
  } else {
    list(Nmc = 8000, Kdir = 4)
  }
}

#-------------
# Plot styling 
#-------------
curve_labels <- c(
  "Y",
  "BNN with fixed scales",
  "BNN with BetaPrime hyper prior",
  "BNN with dropout"
)
curve_cols <- c(
  "black",
  "#D55E00",
  "#0072B2",
  "#009E73"
)
curve_lty <- c(1, 2, 3, 4)
curve_pch <- c(1, 16, 17, 15)

# Add a line with periodic point markers for readability
add_styled_curve <- function(x, y, col, lty, pch, mark_every = 8) {
  lines(x, y, col = col, lty = lty, lwd = 2.5)
  idx <- seq(1, length(x), by = mark_every)
  points(x[idx], y[idx], col = col, pch = pch, cex = 0.8)
}

#----------------------------------------------------------------
# Precompute shared Monte Carlo samplers from the BNN mixing laws
#----------------------------------------------------------------

V_bnn_shared <- rV_bnn_fixed(
  M = M_V,
  d = d,
  widths = widths,
  sigmas = sigmas
)
V_bnn_do_shared <- rV_bnn_dropout(
  M = M_V,
  d = d,
  widths = widths,
  sigmas = sigmas,
  keep_probs = keep_probs,
  inverted = TRUE
)
#------------------------------------------
# Main routine for one dimension p: 
# 1. precompute shrinkage curves, 
# 2. estimate risk over the r-grid, 
# 3. return all objects needed for plotting 
#------------------------------------------

run_one_p <- function(p, r_grid, V_bnn, V_bnn_do) {
  mc <- mc_settings_for_p(p)
  
  # Grid for s=||y||^2 used in shrinkage interpolation 
  s_max <- (max(r_grid) + 6 * sqrt(p))^2
  s_grid <- seq(0, s_max, length.out = 2500)
  
  # Fixed-scale BNN shrinkage 
  a_bnn_grid <- precompute_shrink_from_V(V_bnn, s_grid, p)
  shrink_bnn_fixed <- function(s) interp1(s_grid, a_bnn_grid, s)
  
  # Dropout BNN shrinkage 
  a_bnn_do_grid <- precompute_shrink_from_V(V_bnn_do, s_grid, p)
  shrink_bnn_dropout <- function(s) interp1(s_grid, a_bnn_do_grid, s)
  
  # Beta-prime shrinkage 
  a_mm_grid <- shrink_betaprime_minimax(s_grid, p)
  shrink_minimax <- function(s) shrink_betaprime_minimax(s, p)
  
  # Risk of the MLE, Y under squared error loss
  risk_mle <- rep(p, length(r_grid))
  
  # Estimate risks by Monte Carlo, averaging over random signal directions
  risk_bnn <- sapply(
    r_grid,
    function(r) mean(replicate(
      mc$Kdir,
      risk_at_norm_shrink(r, p, mc$Nmc, shrink_bnn_fixed)
    ))
  )
  risk_mm <- sapply(
    r_grid,
    function(r) mean(replicate(
      mc$Kdir,
      risk_at_norm_shrink(r, p, mc$Nmc, shrink_minimax)
    ))
  )
  risk_bnn_do <- sapply(
    r_grid,
    function(r) mean(replicate(
      mc$Kdir,
      risk_at_norm_shrink(r, p, mc$Nmc, shrink_bnn_dropout)
    ))
  )
  list(
    p = p,
    s_grid = s_grid,
    r_grid = r_grid,
    a_bnn_grid = a_bnn_grid,
    a_mm_grid = a_mm_grid,
    a_bnn_do_grid = a_bnn_do_grid,
    risk_mle = risk_mle,
    risk_bnn = risk_bnn,
    risk_mm = risk_mm,
    risk_bnn_do = risk_bnn_do
  )
}

#----------------------------------------------------------------------
# Run the simulation for each p and produce the corresponding risk plot
#----------------------------------------------------------------------

results <- vector("list", length(p_vals))
names(results) <- paste0("p_", p_vals)
for (i in seq_along(p_vals)) {
  p <- p_vals[i]
  message("Running p = ", p, " ...")
  res <- run_one_p(
    p = p,
    r_grid = r_grid,
    V_bnn = V_bnn_shared,
    V_bnn_do = V_bnn_do_shared
  )
  results[[i]] <- res
  risk_all <- c(
    res$risk_mle,
    res$risk_bnn,
    res$risk_mm,
    res$risk_bnn_do
  )
  plot(
    res$r_grid, res$risk_mle,
    type = "l",
    lwd = 2.5,
    col = curve_cols[1],
    lty = curve_lty[1],
    xlab = expression(paste("\u2016", theta, "\u2016")),
    ylab = "Estimated risk",
    main = paste("Risk comparison (p =", p, ")"),
    ylim = range(risk_all)
  )
  idx <- seq(1, length(res$r_grid), by = 8)
  points(
    res$r_grid[idx], res$risk_mle[idx],
    col = curve_cols[1],
    pch = curve_pch[1],
    cex = 0.8
  )
  add_styled_curve(res$r_grid, res$risk_bnn,    curve_cols[2], curve_lty[2], curve_pch[2])
  add_styled_curve(res$r_grid, res$risk_mm,     curve_cols[3], curve_lty[3], curve_pch[3])
  add_styled_curve(res$r_grid, res$risk_bnn_do, curve_cols[4], curve_lty[4], curve_pch[4])
  legend(
    "bottomright",
    legend = curve_labels,
    col = curve_cols,
    lty = curve_lty,
    pch = curve_pch,
    lwd = 2.5,
    pt.cex = 0.8,
    cex = 0.8,
    bty = "n"
  )
}