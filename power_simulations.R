###############################################################################
# Power Analysis for Secondary Analysis in No-Report Paradigms
# 
# Design: 
#   Stage 1 (Primary): Each trial is classified as percept A or B using a
#     physiological measure (e.g., OKN, pupil dilation) with accuracy `acc`.
#   Stage 2 (Secondary): A secondary neural variable (e.g., BOLD signal) is
#     compared between classified-A and classified-B trials. For each of n
#     participants, the mean difference D_hat is computed across m trials.
#     A one-sample t-test tests whether the group mean of D_hat differs from 0.
#
# Variance structure:
#   - sigma_within: trial-level (residual) SD of secondary measure
#   - sigma_between: between-subject SD of the true effect Delta
#   - Delta: population mean true effect size (mu_1 - mu_0)
#   - acc: classification accuracy (probability of correct label)
#
# Continuous viewing (Section 4) parameterized in physical units
# (session length, sample rate, mean dominance duration) 
# Two versions: 
#   (i)  power_simulate_autocorr  -- one decision per dominance period. This is
#        the reference that instantiates the assumption behind Eq. (3)/(A5),
#        i.e. m_class = number of periods. Use it to validate the formula
#   (ii) power_simulate_persample -- a realistic decoder that classifies each
#        sample from a distinct feature channel (feature signal + AR(1) feature
#        noise, thresholded). Here `acc` and the effective m_class emerge from
#        the feature's autocorrelation timescale, and Section 7b tests whether
#        the formula's "m_class = periods" is conservative or optimistic.
#
# NOISE TIMESCALE. The number of independent noise samples is set by the noise
# correlation time (seconds). The functions therefore accept a physical 
# `tau_noise_s` that is converted to the AR(1) lag-1 autocorrelation at the given 
# sampling rate via rho = exp(-1/(sampling_hz * tau_noise_s)). 

###############################################################################

rm(list = ls())
setwd("/Users/daniel/Documents/Arbeit/PHD/Research/review/review\ power\ derivation/final")

library(dplyr)   
library(tidyr)  

# ─────────────────────────────────────────────────────────────────────────────
# 1. Analytic power function
# ─────────────────────────────────────────────────────────────────────────────

#' Compute analytic power for the secondary analysis
#'
#' @param n Number of participants
#' @param m Number of (effective) trials per participant
#' @param acc Classification accuracy (0.5 = chance, 1 = perfect)
#' @param delta True effect size (mu_1 - mu_0) on the secondary measure
#' @param sigma_within Within-subject, trial-level SD of secondary measure
#' @param sigma_between Between-subject SD of true effect (default 0)
#' @param alpha Significance level (two-sided)
#' @return Analytic power estimate
power_analytic <- function(n, m, acc, delta, sigma_within, 
                           sigma_between = 0, alpha = 0.05) {
  
  eps <- 1 - acc
  
  # Expected observed difference per subject
  E_Dhat <- (2 * acc - 1) * delta
  
  # Within-subject variance of D_hat (across m trials, balanced classes)
  s2 <- sigma_within^2 + eps * (1 - eps) * (delta^2 + sigma_between^2)
  var_Dhat_within <- 4 * s2 / m
  
  # Total variance of D_hat across subjects includes between-subject variance
  # in the true effect. If sigma_between > 0, each subject has their own delta_i
  # drawn from N(delta, sigma_between^2), so the observed D_hat_i has additional
  # variance from (2*acc - 1)^2 * sigma_between^2
  var_Dhat_total <- var_Dhat_within + (2 * acc - 1)^2 * sigma_between^2
  
  SD_Dhat <- sqrt(var_Dhat_total)
  
  # Group-level SE and noncentrality parameter
  SE_group <- SD_Dhat / sqrt(n)
  lambda <- E_Dhat / SE_group
  
  # Power: P(|T| > t_crit) under noncentral t
  df <- n - 1
  t_crit <- qt(1 - alpha / 2, df)
  
  power <- pt(-t_crit, df, ncp = lambda) + 1 - pt(t_crit, df, ncp = lambda)
  return(power)
}


# ─────────────────────────────────────────────────────────────────────────────
# 1b. Study-planning helper for continuous viewing (two effective counts)
# ─────────────────────────────────────────────────────────────────────────────
#
# For continuous viewing the effective sample size is TWO counts, not one
# (cf. Eq. 3 / A5):
#   * m_class = dominance periods per subject. Governs the classification-
#               attenuation term (one decode decision per percept).
#               Estimate as  session_time / mean_dominance_duration.
#               From active-report pilots this is optimistic (report speeds
#               reversals -> too many periods); deflate it or treat as a bound.
#   * m_noise = effective independent noise samples. Governs the measurement-
#               variance term. Estimate from the residual autocorrelation of the
#               secondary measure:  m_noise = m / tau_int,
#               tau_int = 1 + 2*sum_k rho_k  (integrated autocorrelation time);
#               for AR(1):  m_noise = m*(1-rho)/(1+rho),  m = time * sampling_rate.
#               Estimate rho on RESIDUALS, not raw signal (the piecewise-constant
#               percept inflates raw autocorrelation).
#               Bounds:  m_class <= m_noise <= m  (slow/BOLD -> m_class;
#               fast eye/M-EEG -> ~m).

#' Two-component analytic power (continuous viewing), Eq. (3)/(A5).
#' Reduces to power_analytic() when m_class == m_noise == m.
power_analytic_2comp <- function(n, m_class, m_noise, acc, delta, sigma_within,
                                 sigma_between = 0, alpha = 0.05) {
  eps    <- 1 - acc
  E_Dhat <- (2 * acc - 1) * delta
  # noise term uses m_noise; classification-attenuation term uses m_class
  class_term      <- eps * (1 - eps) * (delta^2 + sigma_between^2)
  var_Dhat_within <- 4 * (sigma_within^2 / m_noise + class_term / m_class)
  var_Dhat_total  <- var_Dhat_within + (2 * acc - 1)^2 * sigma_between^2
  lambda <- E_Dhat / sqrt(var_Dhat_total / n)
  df     <- n - 1
  t_crit <- qt(1 - alpha / 2, df)
  pt(-t_crit, df, ncp = lambda) + 1 - pt(t_crit, df, ncp = lambda)
}


power_analytic_2comp(10, 10, 20, 0.7, 0.2, 0.1, 0.1)

#' Smallest number of subjects reaching `target_power` for given counts.
#' Returns NA if unreachable within n_max (e.g. a tau-limited design).
required_n_2comp <- function(m_class, m_noise, acc, delta, sigma_within,
                             sigma_between = 0, target_power = 0.8,
                             alpha = 0.05, n_max = 1000) {
  for (n in 2:n_max) {
    if (power_analytic_2comp(n, m_class, m_noise, acc, delta, sigma_within,
                             sigma_between, alpha) >= target_power)
      return(n)
  }
  NA
}

required_n_2comp(10, 20, 0.6, 0.2, 0.1, 0.1)

#' Convert a physical noise correlation time (seconds) to the AR(1) lag-1
#' autocorrelation at a given sampling rate:  rho = exp(-1 / (fs * tau_s)).
#'
#' This is the sample-rate-invariant way to specify temporal dependence. The
#' number of independent noise samples in a fixed-length recording is fixed by
#' the correlation time, not by the polling density: doubling the sample rate
#' interleaves more highly-correlated samples, it does not create independent
#' information. Holding a fixed `rho` while varying `sampling_hz` (the old
#' behaviour) therefore manufactures independent samples out of nothing. With
#' rho derived here, m_noise = m*(1-rho)/(1+rho) saturates as the rate grows,
#' approaching m_noise -> T_s / (2 * tau_s) in the well-sampled limit.
rho_from_tau <- function(tau_noise_s, sampling_hz) {
  exp(-1 / (sampling_hz * tau_noise_s))
}

#' Plan a continuous-viewing no-report study.
#'
#' Physical inputs -> effective counts -> required sample size, with a
#' conservative/optimistic bracket on m_noise.
#'
#' Effect is standardized: d = Delta/sigma (per trial) and tau = between-subject
#' SD of that standardized effect (so internally sigma_within = 1, delta = d,
#' sigma_between = tau). Any argument may be a vector; all combinations are
#' evaluated (expand.grid) for sensitivity analysis.
#'
#' NOISE PARAMETERIZATION. m_noise is the number of *independent* noise samples,
#' which is set by the noise correlation time (seconds). Supply `tau_noise_s`; the 
#' lag-1 autocorrelation is derived per row as rho = exp(-1 / (sampling_hz * tau_noise_s)) 
#' (see rho_from_tau). 
#'
#' @param session_min  viewing time per subject (minutes)
#' @param mean_dom_s   mean dominance (percept) duration (seconds)
#' @param sampling_hz  effective sampling rate of the secondary measure (Hz)
#' @param tau_noise_s  residual-noise correlation time (seconds); preferred
#' @param rho          fixed lag-1 autocorrelation (legacy; overrides
#'                     tau_noise_s and makes m_noise sample-rate dependent)
#' @param acc          expected primary classification accuracy
#' @param d            per-trial standardized effect (Delta/sigma); default 0.5
#' @param tau          between-subject SD of the standardized effect; default 0.2
#' @param target_power desired power; default 0.8
#' @return data.frame: derived rho, m, m_class, m_noise, and required n at the
#'   point estimate of m_noise plus conservative (m_noise = m_class) and
#'   optimistic (m_noise = m) brackets.
plan_continuous <- function(session_min, mean_dom_s, sampling_hz,
                            tau_noise_s = NULL, rho = NULL,
                            acc, d = 0.5, tau = 0.2,
                            target_power = 0.8, alpha = 0.05, n_max = 1000) {
  
  use_rho <- !is.null(rho)
  if (use_rho && !is.null(tau_noise_s))
    warning("Both `rho` and `tau_noise_s` supplied; using fixed `rho` and ",
            "ignoring `tau_noise_s`. A fixed rho makes m_noise scale with ",
            "sampling_hz, which is a modelling artifact.")
  if (!use_rho && is.null(tau_noise_s))
    stop("Supply `tau_noise_s` (noise correlation time in seconds) or `rho`.")
  
  base <- list(session_min = session_min, mean_dom_s = mean_dom_s,
               sampling_hz = sampling_hz, acc = acc, d = d, tau = tau)
  base[[if (use_rho) "rho" else "tau_noise_s"]] <-
    if (use_rho) rho else tau_noise_s
  g <- do.call(expand.grid, c(base, list(KEEP.OUT.ATTRS = FALSE)))
  
  T_s       <- g$session_min * 60
  g$m       <- round(T_s * g$sampling_hz)               # total samples
  g$m_class <- pmax(2, round(T_s / g$mean_dom_s))       # dominance periods
  
  # Derive rho from the physical correlation time so m_noise does not depend on
  # sampling_hz above the resolving rate (preferred); or use the fixed rho.
  if (!use_rho)
    g$rho <- rho_from_tau(g$tau_noise_s, g$sampling_hz)
  
  g$m_noise <- pmax(g$m_class, round(g$m * (1 - g$rho) / (1 + g$rho)))
  
  g$n_point   <- mapply(required_n_2comp, g$m_class, g$m_noise, g$acc, g$d,
                        1, g$tau, target_power, alpha, n_max)
  g$n_conserv <- mapply(required_n_2comp, g$m_class, g$m_class, g$acc, g$d,
                        1, g$tau, target_power, alpha, n_max)  # m_noise = m_class
  g$n_optim   <- mapply(required_n_2comp, g$m_class, g$m,       g$acc, g$d,
                        1, g$tau, target_power, alpha, n_max)  # m_noise = m
  g
}

# Examples (now in physical noise-timescale units):
#   tau_noise_s = 9.5   reproduces the old rho = 0.9 at 1 Hz;
#   tau_noise_s = 0.0025 reproduces the old rho = 0.2 at 250 Hz.
plan_continuous(session_min = 20, mean_dom_s = c(2, 4, 8),
                sampling_hz = 1,   tau_noise_s = 9.5,    acc = 0.75)  # slow / BOLD-like
plan_continuous(session_min = 20, mean_dom_s = 4,
                sampling_hz = 250, tau_noise_s = 0.0025, acc = 0.75)  # fast / M-EEG-like

# Verify the fix: with tau_noise_s fixed, m_noise SATURATES with sampling rate
# (approaching T_s/(2*tau_noise_s)) instead of growing without bound.
plan_continuous(session_min = 20, mean_dom_s = 4,
                sampling_hz = c(0.5, 1, 2, 5, 20), tau_noise_s = 1, acc = 0.75)


# ─────────────────────────────────────────────────────────────────────────────
# 2. Simulation function
# ─────────────────────────────────────────────────────────────────────────────

#' Simulate power via Monte Carlo
#'
#' @param n Number of participants
#' @param m Number of trials per participant (total; split ~equally)
#' @param acc Classification accuracy
#' @param delta True effect size
#' @param sigma_within Trial-level SD
#' @param sigma_between Between-subject SD of effect (default 0)
#' @param alpha Significance level
#' @param n_sim Number of Monte Carlo replications
#' @param balanced If TRUE, true classes are exactly balanced
#' @param return_pvalues If TRUE, return the full vector of p-values instead of
#'   the scalar power (useful for null calibration checks against Uniform(0,1)).
#' @return Simulated power (proportion of significant tests), or the p-value
#'   vector if return_pvalues = TRUE.
power_simulate <- function(n, m, acc, delta, sigma_within,
                           sigma_between = 0, alpha = 0.05, 
                           n_sim = 5000, balanced = TRUE,
                           return_pvalues = FALSE) {
  
  p_values <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    Dhat <- numeric(n)
    
    for (i in 1:n) {
      # Subject-specific true effect
      delta_i <- delta + rnorm(1, 0, sigma_between)
      
      # True labels: balanced classes
      if (balanced) {
        m0 <- floor(m / 2)
        m1 <- m - m0
      } else {
        # Random class proportions per subject (mildly unbalanced)
        m1 <- rbinom(1, m, 0.5)
        m0 <- m - m1
        if (m0 == 0) { m0 <- 1; m1 <- m - 1 }
        if (m1 == 0) { m1 <- 1; m0 <- m - 1 }
      }
      
      # Generate secondary measure for each trial
      x0 <- rnorm(m0, mean = 0, sd = sigma_within)            # true class 0
      x1 <- rnorm(m1, mean = delta_i, sd = sigma_within)      # true class 1
      
      # True labels
      y_true <- c(rep(0, m0), rep(1, m1))
      x_all  <- c(x0, x1)
      
      # Misclassify with probability (1 - acc), independently
      flip <- rbinom(m, 1, 1 - acc)
      y_hat <- ifelse(flip == 1, 1 - y_true, y_true)
      
      # Compute mean difference based on classified labels
      idx1 <- which(y_hat == 1)
      idx0 <- which(y_hat == 0)
      
      if (length(idx1) < 2 || length(idx0) < 2) {
        # Degenerate case: skip (treat as non-significant)
        Dhat[i] <- 0
      } else {
        Dhat[i] <- mean(x_all[idx1]) - mean(x_all[idx0])
      }
    }
    
    # One-sample t-test: is the group mean of Dhat different from 0?
    tt <- t.test(Dhat, mu = 0)
    p_values[sim] <- tt$p.value
  }
  
  if (return_pvalues) return(p_values)
  mean(p_values < alpha)
}


# ─────────────────────────────────────────────────────────────────────────────
# 3. Simulation with correlated classification errors (assumption violation)
# ─────────────────────────────────────────────────────────────────────────────

#' Simulate power when classification errors correlate with secondary measure
#' 
#' This violates the assumption that misclassification is independent of x.
#' Trials closer to the class boundary (x near grand mean) are more likely 
#' to be misclassified.
#'
#' @inheritParams power_simulate
#' @param return_pvalues If TRUE, return the full vector of p-values.
power_simulate_correlated <- function(n, m, acc, delta, sigma_within,
                                      sigma_between = 0, alpha = 0.05,
                                      n_sim = 5000, return_pvalues = FALSE) {
  p_values <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    Dhat <- numeric(n)
    
    for (i in 1:n) {
      delta_i <- delta + rnorm(1, 0, sigma_between)
      m0 <- floor(m / 2)
      m1 <- m - m0
      
      x0 <- rnorm(m0, mean = 0, sd = sigma_within)
      x1 <- rnorm(m1, mean = delta_i, sd = sigma_within)
      
      y_true <- c(rep(0, m0), rep(1, m1))
      x_all  <- c(x0, x1)
      
      # Misclassification probability depends on distance from boundary
      # Boundary is at delta_i / 2 (midpoint of the two class means)
      boundary <- delta_i / 2
      dist_from_boundary <- abs(x_all - boundary)
      
      # Sigmoid: trials near boundary have higher error rate
      # Scale so that average error ~ (1 - acc)
      # Use logistic function: p(error) = 1 / (1 + exp(k * dist))
      # Calibrate k so that mean error rate ≈ (1 - acc)
      target_eps <- 1 - acc
      
      # Binary search for k
      k_low <- 0; k_high <- 20
      for (iter in 1:30) {
        k_mid <- (k_low + k_high) / 2
        p_err <- 1 / (1 + exp(k_mid * dist_from_boundary))
        if (mean(p_err) > target_eps) k_low <- k_mid else k_high <- k_mid
      }
      p_err <- 1 / (1 + exp(k_mid * dist_from_boundary))
      
      flip <- rbinom(m, 1, p_err)
      y_hat <- ifelse(flip == 1, 1 - y_true, y_true)
      
      idx1 <- which(y_hat == 1)
      idx0 <- which(y_hat == 0)
      
      if (length(idx1) < 2 || length(idx0) < 2) {
        Dhat[i] <- 0
      } else {
        Dhat[i] <- mean(x_all[idx1]) - mean(x_all[idx0])
      }
    }
    
    tt <- t.test(Dhat, mu = 0)
    p_values[sim] <- tt$p.value
  }
  
  if (return_pvalues) return(p_values)
  mean(p_values < alpha)
}


#' Simulate power when classification errors concentrate FAR from the boundary
#' (the inversion of power_simulate_correlated).
#'
#' Where power_simulate_correlated misclassifies ambiguous trials (x near the
#' class midpoint). This is a benign error pattern, because flipping a near-central
#' value barely moves the group means. 
#'
#' Error probability is made proportional to distance from the boundary and
#' rescaled so the mean error rate matches (1 - acc), clipped to [0, 1].
#'
#' @inheritParams power_simulate
#' @param return_pvalues If TRUE, return the full vector of p-values.
power_simulate_correlated_far <- function(n, m, acc, delta, sigma_within,
                                          sigma_between = 0, alpha = 0.05,
                                          n_sim = 5000, return_pvalues = FALSE) {
  p_values   <- numeric(n_sim)
  target_eps <- 1 - acc
  
  for (sim in 1:n_sim) {
    Dhat <- numeric(n)
    
    for (i in 1:n) {
      delta_i <- delta + rnorm(1, 0, sigma_between)
      m0 <- floor(m / 2)
      m1 <- m - m0
      
      x0 <- rnorm(m0, mean = 0, sd = sigma_within)
      x1 <- rnorm(m1, mean = delta_i, sd = sigma_within)
      
      y_true <- c(rep(0, m0), rep(1, m1))
      x_all  <- c(x0, x1)
      
      # Error probability INCREASES with distance from the boundary.
      boundary           <- delta_i / 2
      dist_from_boundary <- abs(x_all - boundary)
      w <- dist_from_boundary / mean(dist_from_boundary)   # mean 1, larger when far
      p_err <- pmin(target_eps * w, 0.5) #  Cap per-trial error at 0.5
      # one rescaling pass so the realised mean error stays ~ target after clip
      mp <- mean(p_err)
      if (mp > 0) p_err <- pmin(p_err * (target_eps / mp), 0.5)
      
      flip  <- rbinom(m, 1, p_err)
      y_hat <- ifelse(flip == 1, 1 - y_true, y_true)
      
      idx1 <- which(y_hat == 1)
      idx0 <- which(y_hat == 0)
      
      if (length(idx1) < 2 || length(idx0) < 2) {
        Dhat[i] <- 0
      } else {
        Dhat[i] <- mean(x_all[idx1]) - mean(x_all[idx0])
      }
    }
    
    tt <- t.test(Dhat, mu = 0)
    p_values[sim] <- tt$p.value
  }
  
  if (return_pvalues) return(p_values)
  mean(p_values < alpha)
}


# ─────────────────────────────────────────────────────────────────────────────
# 4. Continuous-viewing simulations (physical units) Both functions take session
# length, sample rate and mean dominance duration, and share the same
# percept-sequence generator, so they are directly comparable. They differ only
# in how the decoder labels the data.
#
# Temporal dependence may be specified either as a fixed lag-1 rho / feat_rho
# (default) OR, sample-rate-invariantly, as a correlation time in seconds via
# tau_noise_s / feat_tau_s (converted with rho_from_tau). Within a single
# fixed-rate scenario the two are equivalent; the tau_* form matters only if you
# sweep sampling_hz (see plan_continuous for the rationale).
# ─────────────────────────────────────────────────────────────────────────────

#' Fast stationary AR(1) series with the correct marginal SD from sample 1.
#' Uses stats::filter (C-level) instead of an R for-loop.
#' @param m length; @param rho lag-1 correlation; @param marg_sd marginal SD.
ar1_series <- function(m, rho, marg_sd = 1) {
  if (m <= 0) return(numeric(0))
  innov      <- rnorm(m, 0, marg_sd * sqrt(1 - rho^2))
  innov[1]   <- rnorm(1, 0, marg_sd)                 # stationary start (Var = marg_sd^2)
  as.numeric(stats::filter(innov, filter = rho, method = "recursive"))
}

#' AR(1) series with a non-Gaussian, standardized marginal (NORTA / Gaussian
#' copula). A Gaussian AR(1) supplies the temporal dependence; its marginal is
#' pushed through the probability integral transform to the target shape, then
#' standardized to mean 0 and SD `marg_sd`. The realized linear autocorrelation is a
#' little below rho (copula attenuation), but strong dependence is retained.
#'
#' @param dist "gaussian", "t" (heavy-tailed, symmetric) or "lognorm" (skewed).
#' @param df   degrees of freedom for the t marginal (must be > 2).
#' @param sdlog shape of the log-normal marginal (larger = more skew; 
#' 3 roughly reported actual shape in the literature).
ar1_series_marg <- function(m, rho, marg_sd = 1,
                            dist = c("gaussian", "t", "lognorm"),
                            df = 3, sdlog = 1) {
  dist <- match.arg(dist)
  z <- ar1_series(m, rho, 1)                 # standard Gaussian AR(1), marginal N(0,1)
  if (dist == "gaussian") return(z * marg_sd)
  u <- pmin(pmax(pnorm(z), 1e-6), 1 - 1e-6)  # uniform via copula (clamped for safety)
  x <- switch(dist,
              t = {                                     # Student-t, standardized to unit variance
                stopifnot(df > 2)
                qt(u, df) / sqrt(df / (df - 2))
              },
              lognorm = {                               # log-normal, centred & scaled to unit variance
                raw     <- qlnorm(u, meanlog = 0, sdlog = sdlog)
                mu_raw  <- exp(sdlog^2 / 2)
                var_raw <- (exp(sdlog^2) - 1) * exp(sdlog^2)
                (raw - mu_raw) / sqrt(var_raw)
              }
  )
  x * marg_sd
}


#' Draw an alternating dominance-period sequence covering >= m samples.
#' Period lengths (in samples) are Gamma(shape, scale) with the requested mean.
#' Returns a list with the per-sample true label y, the sign code d = 2y-1, and
#' the per-sample period id.
draw_percept_sequence <- function(m, mean_period_samp, gamma_shape = 4) {
  scale     <- mean_period_samp / gamma_shape
  durations <- integer(0); total <- 0L
  while (total < m) {
    dur       <- max(1L, round(rgamma(1, shape = gamma_shape, scale = scale)))
    durations <- c(durations, dur); total <- total + dur
  }
  start <- sample(0:1, 1)
  labs  <- (start + seq_along(durations) - 1) %% 2       # strict alternation
  y     <- rep(labs, durations)[1:m]
  pid   <- rep(seq_along(durations), durations)[1:m]
  list(y = y, d = 2 * y - 1, period_id = pid, n_periods = length(durations))
}

hist(diff(which(abs(diff(draw_percept_sequence(1000, 5)$y)) == 1)))

#' (i) Continuous viewing, one decision per dominance period.
#'
#' Reference simulation that instantiates the Eq. (3)/(A5) assumption
#' m_class = number of periods: the decoder emits a single label per percept,
#' correct with probability `acc`, held constant across the period. The
#' secondary measure carries AR(1) noise (lag-1 `rho`) so m_noise behaves as
#' m*(1-rho)/(1+rho). Matches power_analytic_2comp with those two counts.
#'
#' @param n subjects
#' @param session_min viewing time per subject (minutes)
#' @param sampling_hz sample rate of the secondary measure (Hz)
#' @param mean_dom_s mean dominance duration (seconds)
#' @param acc per-period classification accuracy (INPUT here)
#' @param delta per-trial effect; @param sigma_within marginal noise SD
#' @param rho lag-1 autocorrelation of the AR(1) secondary-measure noise
#' @param sigma_between between-subject SD of the effect
#' @param gamma_shape shape of the period-length distribution
#' @param return_pvalues if TRUE return the p-value vector (NAs retained)
#' @param tau_noise_s optional noise correlation time (s); if given, overrides
#'   `rho` via rho = exp(-1/(sampling_hz * tau_noise_s)).
power_simulate_autocorr <- function(n, session_min, sampling_hz, mean_dom_s,
                                    acc, delta, sigma_within,
                                    rho = 0.3, sigma_between = 0,
                                    gamma_shape = 2, alpha = 0.05,
                                    n_sim = 2000, return_pvalues = FALSE,
                                    noise_dist = "gaussian", noise_df = 3,
                                    tau_noise_s = NULL) {
  session_s        <- session_min * 60
  m                <- round(session_s * sampling_hz)
  mean_period_samp <- mean_dom_s * sampling_hz
  if (!is.null(tau_noise_s)) rho <- rho_from_tau(tau_noise_s, sampling_hz)
  p_values         <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    Dhat <- numeric(n)
    for (i in 1:n) {
      delta_i <- delta + rnorm(1, 0, sigma_between)
      seqn    <- draw_percept_sequence(m, mean_period_samp, gamma_shape)
      
      # secondary measure = percept effect + AR(1) noise (marginal per noise_dist)
      x <- seqn$y * delta_i +
        ar1_series_marg(m, rho, sigma_within, noise_dist, noise_df)
      
      # ONE decision per period, correct w.p. acc, persisting within the period
      period_flip <- rbinom(seqn$n_periods, 1, 1 - acc)
      flip        <- period_flip[seqn$period_id]
      y_hat       <- ifelse(flip == 1, 1 - seqn$y, seqn$y)
      
      idx1 <- which(y_hat == 1); idx0 <- which(y_hat == 0)
      Dhat[i] <- if (length(idx1) < 2 || length(idx0) < 2) NA_real_
      else mean(x[idx1]) - mean(x[idx0])
    }
    Dhat <- Dhat[is.finite(Dhat)]
    if (length(Dhat) < 2) { p_values[sim] <- NA_real_; next }
    p_values[sim] <- t.test(Dhat, mu = 0)$p.value
  }
  
  if (return_pvalues) return(p_values)
  mean(p_values < alpha, na.rm = TRUE)
}

#' (ii) Continuous viewing, PER-SAMPLE decoder on a distinct feature channel.
#'
#' The decoder does not know the period boundaries. It classifies each sample by
#' thresholding a separate feature channel:
#'      f_t = d_t * (b/2) + feature_noise_t ,   feature_noise ~ AR(1, feat_rho),
#' with marginal SD 1, so the MARGINAL per-sample accuracy is pnorm(b/2). To hit
#' a target marginal accuracy `acc`, set b = 2*qnorm(acc). The feature channel is
#' independent of the secondary measure (avoids the correlated-error confound),
#' and its autocorrelation `feat_rho` controls how fast decoding errors
#' decorrelate:
#'   * feat_rho -> 0            : per-sample errors independent -> m_class ~ m
#'                                (formula's "m_class = periods" is CONSERVATIVE)
#'   * feat_rho ~ period scale  : errors decorrelate on the percept timescale
#'                                -> m_class ~ number of periods (formula OK)
#' `rho` still governs the secondary-measure noise (-> m_noise), separately.
#'
#' @inheritParams power_simulate_autocorr
#' @param acc target marginal per-sample accuracy (sets the feature separation b)
#' @param feat_rho lag-1 autocorrelation of the AR(1) feature noise
#' @param tau_noise_s optional secondary-noise correlation time (s); overrides `rho`.
#' @param feat_tau_s  optional feature-noise correlation time (s); overrides `feat_rho`.
power_simulate_persample <- function(n, session_min, sampling_hz, mean_dom_s,
                                     acc, delta, sigma_within,
                                     rho = 0.3, feat_rho = 0.3,
                                     sigma_between = 0, gamma_shape = 2,
                                     alpha = 0.05, n_sim = 2000,
                                     return_pvalues = FALSE,
                                     noise_dist = "gaussian", noise_df = 3,
                                     tau_noise_s = NULL, feat_tau_s = NULL) {
  session_s        <- session_min * 60
  m                <- round(session_s * sampling_hz)
  mean_period_samp <- mean_dom_s * sampling_hz
  b                <- 2 * qnorm(acc)                 # feature separation for target acc
  if (!is.null(tau_noise_s)) rho      <- rho_from_tau(tau_noise_s, sampling_hz)
  if (!is.null(feat_tau_s))  feat_rho <- rho_from_tau(feat_tau_s,  sampling_hz)
  p_values         <- numeric(n_sim)
  
  for (sim in 1:n_sim) {
    Dhat <- numeric(n)
    for (i in 1:n) {
      delta_i <- delta + rnorm(1, 0, sigma_between)
      seqn    <- draw_percept_sequence(m, mean_period_samp, gamma_shape)
      
      # secondary measure (independent AR(1) noise; marginal per noise_dist)
      x <- seqn$y * delta_i +
        ar1_series_marg(m, rho, sigma_within, noise_dist, noise_df)
      
      # per-sample classification on a DISTINCT (Gaussian) feature channel
      f     <- seqn$d * (b / 2) + ar1_series(m, feat_rho, 1)
      y_hat <- as.integer(f > 0)
      
      idx1 <- which(y_hat == 1); idx0 <- which(y_hat == 0)
      Dhat[i] <- if (length(idx1) < 2 || length(idx0) < 2) NA_real_
      else mean(x[idx1]) - mean(x[idx0])
    }
    Dhat <- Dhat[is.finite(Dhat)]
    if (length(Dhat) < 2) { p_values[sim] <- NA_real_; next }
    p_values[sim] <- t.test(Dhat, mu = 0)$p.value
  }
  
  if (return_pvalues) return(p_values)
  mean(p_values < alpha, na.rm = TRUE)
}


#' Runs `n_seq` feature-channel sequences and measures, on average:
#'   - emergent_acc   : realised marginal per-sample accuracy (should ~ acc)
#'   - n_periods      : dominance periods per session (the formula's m_class)
#'   - tau_correctness: integrated autocorrelation time of the correctness sign
#'                      c_t = dhat_t * d_t  (initial-positive-sequence estimator)
#'   - m_class_eff    : m / tau_correctness --> effective classification count
#'   - m_noise_eff    : m*(1-rho)/(1+rho) --> effective noise count
#' Compare m_class_eff against n_periods to see whether "m_class = periods"
#' is conservative (m_class_eff > n_periods) or about right.
#'
#' @param tau_noise_s optional secondary-noise correlation time (s); overrides `rho`.
#' @param feat_tau_s  optional feature-noise correlation time (s); overrides `feat_rho`.
measure_persample <- function(session_min, sampling_hz, mean_dom_s, acc,
                              rho, feat_rho, gamma_shape = 2, n_seq = 100,
                              tau_noise_s = NULL, feat_tau_s = NULL) {
  session_s        <- session_min * 60
  m                <- round(session_s * sampling_hz)
  mean_period_samp <- mean_dom_s * sampling_hz
  b                <- 2 * qnorm(acc)
  if (!is.null(tau_noise_s)) rho      <- rho_from_tau(tau_noise_s, sampling_hz)
  if (!is.null(feat_tau_s))  feat_rho <- rho_from_tau(feat_tau_s,  sampling_hz)
  lag_max          <- min(m - 1, 5 * ceiling(mean_period_samp) + 50)
  
  emerg <- numeric(n_seq); tau <- numeric(n_seq); nper <- numeric(n_seq)
  for (s in 1:n_seq) {
    seqn <- draw_percept_sequence(m, mean_period_samp, gamma_shape)
    f    <- seqn$d * (b / 2) + ar1_series(m, feat_rho, 1)
    yhat <- as.integer(f > 0)
    cc   <- (2 * yhat - 1) * seqn$d          # correctness sign in {-1,+1}
    emerg[s] <- mean(yhat == seqn$y)
    nper[s]  <- seqn$n_periods
    if (sd(cc) == 0) { tau[s] <- 1; next }   # all-correct sequence (rare)
    ac  <- acf(cc, lag.max = lag_max, plot = FALSE)$acf[-1]
    neg <- which(ac <= 0)[1]                  # initial positive sequence
    kk  <- if (is.na(neg)) length(ac) else max(0, neg - 1)
    tau[s] <- 1 + 2 * sum(ac[seq_len(kk)])
  }
  data.frame(
    m               = m,
    n_periods       = mean(nper),
    emergent_acc    = mean(emerg),
    tau_correctness = mean(tau),
    m_class_eff     = m / mean(tau),
    m_noise_eff     = m * (1 - rho) / (1 + rho)
  )
}


# ─────────────────────────────────────────────────────────────────────────────
# 5. Run main validation: analytic vs simulation
# ─────────────────────────────────────────────────────────────────────────────

cat("=== Running main validation ===\n")

# Parameter grid
params <- expand.grid(
  n   = c(10, 20, 30),
  m   = c(50, 200),
  acc = seq(0.55, 0.95, by = 0.05),
  delta = 0.5,
  sigma_within = 1,
  sigma_between = 0
)

results <- params
results$power_analytic <- NA
results$power_sim <- NA

set.seed(42)

for (i in 1:nrow(results)) {
  r <- results[i, ]
  
  results$power_analytic[i] <- power_analytic(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, sigma_between = r$sigma_between
  )
  
  results$power_sim[i] <- power_simulate(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, sigma_between = r$sigma_between,
    n_sim = 5000
  )
  
  cat(sprintf("  n=%d, m=%d, acc=%.2f | analytic=%.3f, sim=%.3f\n",
              r$n, r$m, r$acc, results$power_analytic[i], results$power_sim[i]))
}


# ─────────────────────────────────────────────────────────────────────────────
# 6. Between-subject variability scenario
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Between-subject variability ===\n")

params_bw <- expand.grid(
  n = 20,
  m = c(50, 200),
  acc = seq(0.55, 0.95, by = 0.1),
  delta = 0.5,
  sigma_within = 1,
  sigma_between = c(0, 0.25, 0.5)
)

results_bw <- params_bw
results_bw$power_analytic <- NA
results_bw$power_sim <- NA

for (i in 1:nrow(results_bw)) {
  r <- results_bw[i, ]
  results_bw$power_analytic[i] <- power_analytic(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, sigma_between = r$sigma_between
  )
  results_bw$power_sim[i] <- power_simulate(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, sigma_between = r$sigma_between,
    n_sim = 5000
  )
  cat(sprintf("  m=%d, acc=%.2f, sigma_bw=%.2f | analytic=%.3f, sim=%.3f\n",
              r$m, r$acc, r$sigma_between, 
              results_bw$power_analytic[i], results_bw$power_sim[i]))
}


# ─────────────────────────────────────────────────────────────────────────────
# 7. Assumption violations (trial-based: independent, correlated, unbalanced)
#    Continuous viewing moved to its own physical-unit section (7b).
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Assumption violations ===\n")

params_viol <- expand.grid(
  n = 20,
  m = 100,
  acc = c(0.6, 0.7, 0.8, 0.9),
  delta = 0.5,
  sigma_within = 1
)

results_viol <- params_viol
results_viol$power_analytic <- NA
results_viol$power_independent <- NA
results_viol$power_correlated <- NA
results_viol$power_correlated_far <- NA
results_viol$power_unbalanced <- NA

for (i in 1:nrow(results_viol)) {
  r <- results_viol[i, ]
  
  results_viol$power_analytic[i] <- power_analytic(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within
  )
  
  results_viol$power_independent[i] <- power_simulate(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, n_sim = 5000
  )
  
  results_viol$power_correlated[i] <- power_simulate_correlated(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, n_sim = 5000
  )
  
  results_viol$power_correlated_far[i] <- power_simulate_correlated_far(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, n_sim = 5000
  )
  
  results_viol$power_unbalanced[i] <- power_simulate(
    n = r$n, m = r$m, acc = r$acc, delta = r$delta,
    sigma_within = r$sigma_within, n_sim = 5000, balanced = FALSE
  )
  
  cat(sprintf(paste0("  acc=%.2f | analytic=%.3f indep=%.3f corr(near)=%.3f ",
                     "corr(far)=%.3f unbal=%.3f\n"),
              r$acc, results_viol$power_analytic[i], results_viol$power_independent[i],
              results_viol$power_correlated[i], results_viol$power_correlated_far[i],
              results_viol$power_unbalanced[i]))
}


# ─────────────────────────────────────────────────────────────────────────────
# 7b. Continuous viewing: one-per-period reference vs per-sample decoder Does
# the formula's "m_class = number of periods" hold for a decoder that labels
# every sample without knowing the boundaries? We compare, per physical scenario
# and feature timescale:
#       * per-sample SIM power
#       * one-per-period SIM power (uses the same emergent accuracy)
#       * analytic power with m_class = effective count (from correctness ACF)
#       * analytic power with m_class = number of perids (the formula's rule)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Continuous viewing: per-sample vs one-per-period ===\n")

set.seed(7)

cv_delta   <- 0.5
cv_delta   <- 0.2
cv_sigma_w <- 1
cv_tau     <- 0.2      # between-subject SD
cv_acc     <- 0.75     # target marginal per-sample / per-period accuracy
cv_n       <- 15
cv_nsim    <- 800
cv_nseq    <- 120

# Two acquisition regimes; for each, a fast (independent) and a slow
# (period-scale) feature channel. feat_rho_slow is chosen so the feature
# correlation time ~ one dominance period:  feat_rho = exp(-1 / period_samples).
cv_scen <- data.frame(
  scenario    = c("BOLD-like (slow)", "M/EEG-like (fast)"),
  session_min = c(20,   4),
  sampling_hz = c(1,    50),
  mean_dom_s  = c(4,    3),
  rho_sec     = c(0.6,  0.2),   # secondary-measure noise autocorrelation
  stringsAsFactors = FALSE
)
cv_scen$period_samp <- cv_scen$mean_dom_s * cv_scen$sampling_hz
cv_scen$feat_rho_slow <- exp(-1 / cv_scen$period_samp)
cv_scen$feat_rho_fast <- 0

cv_rows <- list()
for (s in 1:nrow(cv_scen)) {
  for (ftype in c("fast", "slow")) {
    sc      <- cv_scen[s, ]
    feat_rho <- if (ftype == "fast") sc$feat_rho_fast else sc$feat_rho_slow
    
    diag <- measure_persample(
      session_min = sc$session_min, sampling_hz = sc$sampling_hz,
      mean_dom_s  = sc$mean_dom_s,  acc = cv_acc,
      rho = sc$rho_sec, feat_rho = feat_rho, n_seq = cv_nseq
    )
    
    p_persample <- power_simulate_persample(
      n = cv_n, session_min = sc$session_min, sampling_hz = sc$sampling_hz,
      mean_dom_s = sc$mean_dom_s, acc = cv_acc, delta = cv_delta,
      sigma_within = cv_sigma_w, rho = sc$rho_sec, feat_rho = feat_rho,
      sigma_between = cv_tau, n_sim = cv_nsim
    )
    p_perperiod <- power_simulate_autocorr(
      n = cv_n, session_min = sc$session_min, sampling_hz = sc$sampling_hz,
      mean_dom_s = sc$mean_dom_s, acc = diag$emergent_acc, delta = cv_delta,
      sigma_within = cv_sigma_w, rho = sc$rho_sec,
      sigma_between = cv_tau, n_sim = cv_nsim
    )
    p_an_eff <- power_analytic_2comp(
      n = cv_n, m_class = diag$m_class_eff, m_noise = diag$m_noise_eff,
      acc = diag$emergent_acc, delta = cv_delta,
      sigma_within = cv_sigma_w, sigma_between = cv_tau
    )
    p_an_per <- power_analytic_2comp(
      n = cv_n, m_class = diag$n_periods, m_noise = diag$m_noise_eff,
      acc = diag$emergent_acc, delta = cv_delta,
      sigma_within = cv_sigma_w, sigma_between = cv_tau
    )
    
    cv_rows[[length(cv_rows) + 1]] <- data.frame(
      scenario = sc$scenario, feature = ftype,
      emergent_acc = diag$emergent_acc, n_periods = diag$n_periods,
      m_class_eff = diag$m_class_eff, m_noise_eff = diag$m_noise_eff,
      power_persample = p_persample, power_perperiod = p_perperiod,
      power_analytic_effcount = p_an_eff, power_analytic_periods = p_an_per
    )
    cat(sprintf(paste0("  %-18s feat=%-4s | acc~%.2f  periods=%.0f  ",
                       "m_class_eff=%.0f  m_noise=%.0f | ps=%.3f pp=%.3f ",
                       "an_eff=%.3f an_per=%.3f\n"),
                sc$scenario, ftype, diag$emergent_acc, diag$n_periods,
                diag$m_class_eff, diag$m_noise_eff,
                p_persample, p_perperiod, p_an_eff, p_an_per))
  }
}
cv_results <- do.call(rbind, cv_rows)


# ─────────────────────────────────────────────────────────────────────────────
# 8. Type I error rate check (delta = 0)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Type I error rate (delta = 0) ===\n")

type1_params <- expand.grid(
  n = c(10, 20, 30),
  m = 100,
  acc = c(0.6, 0.8, 0.95)
)

for (i in 1:nrow(type1_params)) {
  r <- type1_params[i, ]
  rate <- power_simulate(
    n = r$n, m = r$m, acc = r$acc, delta = 0,
    sigma_within = 1, n_sim = 10000
  )
  cat(sprintf("  n=%d, acc=%.2f | Type I error = %.4f\n", r$n, r$acc, rate))
}


# ─────────────────────────────────────────────────────────────────────────────
# 8b. Null p-value distribution vs. Uniform(0,1) p-value vectors from each Monte
# Carlo variant compared them to the uniform reference via KS test + QQ plot.
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Null p-value distribution vs Uniform(0,1) (delta = 0) ===\n")

set.seed(123)
n_sim_null <- 5000
null_n     <- 20
null_m     <- 100
null_acc   <- 0.75

# physical scenario for the continuous-viewing null entries 
null_session_min <- 5
null_sampling_hz <- 1
null_mean_dom_s  <- 15

pval_null <- list(
  "Independent errors" = power_simulate(
    n = null_n, m = null_m, acc = null_acc, delta = 0, sigma_within = 1,
    n_sim = n_sim_null, return_pvalues = TRUE
  ),
  "Correlated errors" = power_simulate_correlated(
    n = null_n, m = null_m, acc = null_acc, delta = 0, sigma_within = 1,
    n_sim = n_sim_null, return_pvalues = TRUE
  ),
  "Continuous: one-per-period" = power_simulate_autocorr(
    n = null_n, session_min = null_session_min, sampling_hz = null_sampling_hz,
    mean_dom_s = null_mean_dom_s, acc = null_acc, delta = 0, sigma_within = 1,
    rho = 0.6, n_sim = n_sim_null, return_pvalues = TRUE
  ),
  "Continuous: per-sample" = power_simulate_persample(
    n = null_n, session_min = null_session_min, sampling_hz = null_sampling_hz,
    mean_dom_s = null_mean_dom_s, acc = null_acc, delta = 0, sigma_within = 1,
    rho = 0.6, feat_rho = 0.6, n_sim = n_sim_null, return_pvalues = TRUE
  )
)

# One-sample KS test against the uniform CDF, plus the empirical rejection rate.
for (nm in names(pval_null)) {
  p  <- pval_null[[nm]]
  p  <- p[is.finite(p)]
  ks <- suppressWarnings(ks.test(p, "punif"))
  cat(sprintf("  %-28s | KS D = %.4f, p = %.3f | P(p<0.05) = %.4f (N = %d)\n",
              nm, ks$statistic, ks$p.value, mean(p < 0.05), length(p)))
}

# Long data frame: sorted p-values vs their expected uniform quantiles,
# plus the raw values for a histogram.
qq_null <- bind_rows(lapply(names(pval_null), function(nm) {
  p <- sort(pval_null[[nm]][is.finite(pval_null[[nm]])])
  data.frame(condition   = nm,
             theoretical = ppoints(length(p)),   # expected uniform quantiles
             empirical   = p)
}))

hist_null <- bind_rows(lapply(names(pval_null), function(nm) {
  p <- pval_null[[nm]][is.finite(pval_null[[nm]])]
  data.frame(condition = nm, pvalue = p)
}))


# ─────────────────────────────────────────────────────────────────────────────
# 8c. Stress test: non-Gaussian trial noise at low effective count The CLT
# protects D_hat's normality only when the effective sample count is large.
# Short session, a slow feature and slow secondary noise drive the effective
# count down. Furtherk, the secondary-measure marginal is made heavy-tailed (t3)
# or skewed (log-normal). Under the null (delta = 0) a calibrated group t-test
# still yields Uniform(0,1) p-values; deviation here shows exactly where the
# trial-level Gaussian assumption starts to matter. The high-count row is the
# control: with many effective samples the CLT should rescue even t3 noise.
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Null p-values: non-Gaussian noise, low vs high effective count ===\n")

set.seed(2024)
n_sim_stress <- 5000
stress_n     <- 20
stress_acc   <- 0.75

# LOW effective count: short session, slow feature, slow noise, few periods.
low <- list(session_min = 2, sampling_hz = 1, mean_dom_s = 20,
            rho = 0.8, feat_rho = 0.9)
# HIGH effective count control: long session, fast feature, fast noise.
high <- list(session_min = 20, sampling_hz = 1, mean_dom_s = 4,
             rho = 0.2, feat_rho = 0.0)

# Report the effective counts of the low-count scenario so the regime is explicit.
low_diag <- measure_persample(
  session_min = low$session_min, sampling_hz = low$sampling_hz,
  mean_dom_s = low$mean_dom_s, acc = stress_acc,
  rho = low$rho, feat_rho = low$feat_rho, n_seq = 200
)
cat(sprintf(paste0("  LOW-count scenario: m = %.0f, periods = %.1f, ",
                   "m_class_eff = %.1f, m_noise_eff = %.1f\n"),
            low_diag$m, low_diag$n_periods,
            low_diag$m_class_eff, low_diag$m_noise_eff))

run_stress <- function(scen, dist, df = 3) {
  power_simulate_persample(
    n = stress_n, session_min = scen$session_min, sampling_hz = scen$sampling_hz,
    mean_dom_s = scen$mean_dom_s, acc = stress_acc, delta = 0, sigma_within = 1,
    rho = scen$rho, feat_rho = scen$feat_rho, n_sim = n_sim_stress,
    return_pvalues = TRUE, noise_dist = dist, noise_df = df
  )
}

pval_stress <- list(
  "Gaussian, low count"        = run_stress(low,  "gaussian"),
  "Heavy-tailed t3, low count" = run_stress(low,  "t", df = 3),
  "Skewed lognorm, low count"  = run_stress(low,  "lognorm"),
  "Heavy-tailed t3, HIGH count"= run_stress(high, "t", df = 3)
)

for (nm in names(pval_stress)) {
  p  <- pval_stress[[nm]]; p <- p[is.finite(p)]
  ks <- suppressWarnings(ks.test(p, "punif"))
  cat(sprintf("  %-28s | KS D = %.4f, p = %.3f | P(p<0.05) = %.4f (N = %d)\n",
              nm, ks$statistic, ks$p.value, mean(p < 0.05), length(p)))
}

qq_stress <- bind_rows(lapply(names(pval_stress), function(nm) {
  p <- sort(pval_stress[[nm]][is.finite(pval_stress[[nm]])])
  data.frame(condition = nm, theoretical = ppoints(length(p)), empirical = p)
}))
hist_stress <- bind_rows(lapply(names(pval_stress), function(nm) {
  p <- pval_stress[[nm]][is.finite(pval_stress[[nm]])]
  data.frame(condition = nm, pvalue = p)
}))

# Fix a sensible facet/legend order (control last).
stress_levels <- c("Gaussian, low count", "Heavy-tailed t3, low count",
                   "Skewed lognorm, low count", "Heavy-tailed t3, HIGH count")
qq_stress$condition   <- factor(qq_stress$condition,   levels = stress_levels)
hist_stress$condition <- factor(hist_stress$condition, levels = stress_levels)


# ─────────────────────────────────────────────────────────────────────────────
# 9. Plots  

power_given_acc <- function(acc, n, m_eff, Delta, sigma_within = 1, tau = 0.2) {
  power_analytic(n = n, m = m_eff, acc = acc, delta = Delta,
                 sigma_within = sigma_within, sigma_between = tau)
}

trials <- c(50, 200)
cols   <- c("red", "blue")
ns     <- c(10, 20, 30)


# --- Plot 1: Analytic (lines) vs simulated (points), main validation ---------
pdf("plots/plot_validation.pdf", width = 5, height = 5)
v_acc <- seq(0.55, 0.95, 0.001)
plot(NA, axes = FALSE, xlim = c(0.55, 0.95), ylim = c(0, 1),
     xlab = "Primary analysis accuracy", ylab = "Secondary analysis power")
axis(1); axis(2)
abline(h = 0.8, lty = 2, lwd = 1)
for (ti in seq_along(trials)) {
  for (ni in seq_along(ns)) {
    yv <- sapply(v_acc, function(a)
      power_given_acc(a, ns[ni], trials[ti], Delta = 0.5, tau = 0))
    lines(v_acc, yv, col = cols[ti], lty = ni, lwd = 2)
    sub <- results[results$m == trials[ti] & results$n == ns[ni], ]
    points(sub$acc, sub$power_sim, col = cols[ti], pch = ni)
  }
}
legend("bottomright",
       c("50 trials", "200 trials", paste("n =", ns)),
       col = c("red", "blue", rep("black", 3)),
       pch = c(15, 15, seq_along(ns)),
       lty = c(NA, NA, seq_along(ns)),
       lwd = c(NA, NA, rep(2, length(ns))), bg = "white")
title(main = "Analytic (lines) vs. simulation (points)", font.main = 1)
dev.off()


# --- Plot 2: Between-subject variability -------------------------------------
sbs <- sort(unique(results_bw$sigma_between))
pdf("plots/plot_between_subject.pdf", width = 5, height = 5)
v_acc <- seq(0.55, 0.95, 0.001)
plot(NA, axes = FALSE, xlim = c(0.55, 0.95), ylim = c(0, 1),
     xlab = "Primary analysis accuracy", ylab = "Secondary analysis power")
axis(1); axis(2)
abline(h = 0.8, lty = 2, lwd = 1)
for (ti in seq_along(trials)) {
  for (si in seq_along(sbs)) {
    yv <- sapply(v_acc, function(a)
      power_given_acc(a, 20, trials[ti], Delta = 0.5, tau = sbs[si]))
    lines(v_acc, yv, col = cols[ti], lty = si, lwd = 2)
    sub <- results_bw[results_bw$m == trials[ti] &
                        results_bw$sigma_between == sbs[si], ]
    points(sub$acc, sub$power_sim, col = cols[ti], pch = si)
  }
}
legend("bottomright",
       c("50 trials", "200 trials", paste("tau =", sbs)),
       col = c("red", "blue", rep("black", length(sbs))),
       pch = c(15, 15, seq_along(sbs)),
       lty = c(NA, NA, seq_along(sbs)),
       lwd = c(NA, NA, rep(2, length(sbs))), bg = "white")
title(main = "Between-subject variability", font.main = 1)
dev.off()


# --- Plot 3: Assumption violations -------------------------------------------
# Each condition a separate colour; lines connect the discrete accuracy points.
# The two error-correlation conditions are the mirror pair: "near" (benign,
# errors on ambiguous trials) tends to sit at/above analytic, "far" (damaging,
# errors on confident trials) falls below it.
viol_cols  <- c("black", "red", "blue", "purple", "darkgreen")
viol_names <- c("Analytic", "Independent errors",
                "Correlated (near boundary)", "Correlated (far boundary)",
                "Unbalanced classes")
viol_cols_map <- c(power_analytic = viol_cols[1], power_independent = viol_cols[2],
                   power_correlated = viol_cols[3], power_correlated_far = viol_cols[4],
                   power_unbalanced = viol_cols[5])
viol_lty <- c(1, 2, 2, 3, 2)
pdf("plots/plot_violations.pdf", width = 8, height = 6)
plot(NA, axes = FALSE, xlim = c(0.6, 0.9), ylim = c(0, 1),
     xlab = "Classification accuracy", ylab = "Power")
axis(1); axis(2); abline(h = 0.8, lty = 2)
vcols <- names(viol_cols_map)
for (k in seq_along(vcols)) {
  o <- order(results_viol$acc)
  lines(results_viol$acc[o], results_viol[[vcols[k]]][o],
        col = viol_cols_map[k], lwd = 2, lty = viol_lty[k])
  points(results_viol$acc, results_viol[[vcols[k]]],
         col = viol_cols_map[k], pch = 16)
}
legend("bottomright", viol_names, col = viol_cols, lwd = 2,
       lty = viol_lty, pch = 16, bg = "white")
title(main = "Robustness to assumption violations", font.main = 1)
dev.off()


# --- Plot 4: Required accuracy for 0.8 power  -----
Delta_fig4 <- 0.5
tau_fig4   <- 0.2
pdf("plots/power.pdf", width = 5, height = 5)
plot(NA, axes = FALSE, xlim = c(0.5, 0.9), ylim = c(0, 1),
     ylab = "Secondary analysis power", xlab = "Primary analysis accuracy")
v_acc <- seq(0.5, 0.99, 0.001)
acc_vals_08 <- matrix(NA, nrow = length(trials), ncol = length(ns))
y_pos       <- matrix(NA, nrow = length(trials), ncol = length(ns))
dx <- 0.02; dy <- 0.03

for (trial_idx in seq_along(trials)) {
  v_p <- list()
  for (n in ns) {
    vals <- sapply(v_acc, function(a)
      power_given_acc(acc = a, n = n, m_eff = trials[trial_idx],
                      Delta = Delta_fig4, tau = tau_fig4))
    v_p[[toString(n)]] <- vals
  }
  for (n in seq_along(ns)) {
    lines(v_acc, v_p[[toString(ns[n])]], lwd = 2, lty = n, col = cols[trial_idx])
  }
  axis(1); axis(2)
  abline(h = 0.8, lty = 2, lwd = 1)
  acc_vals_08[trial_idx, ] <- sapply(ns, function(x)
    v_acc[which.min(abs(v_p[[toString(x)]] - 0.8))])
  abline(v = acc_vals_08[trial_idx, ], lwd = 1, lty = seq_along(ns),
         col = cols[trial_idx])
  y_pos[trial_idx, ] <- seq(2, 1, length.out = length(ns))^5 / 50
}
for (idx_2 in seq_along(trials)) {
  for (idx_3 in seq_along(acc_vals_08[idx_2, ])) {
    xc <- acc_vals_08[idx_2, idx_3]; yc <- y_pos[idx_2, idx_3]
    polygon(c(xc - dx, xc - dx, xc + dx, xc + dx),
            c(yc - dy, yc + dy, yc + dy, yc - dy),
            col = "white", border = "white")
  }
  text(acc_vals_08[idx_2, ], y = y_pos[idx_2, ],
       round(acc_vals_08[idx_2, ], 2), cex = 1, col = cols[idx_2])
}
legend("bottomright", c("50 trials", "200 trials", paste("n =", ns)),
       col = c("red", "blue", rep("black", 3)),
       pch = c(15, 15, NA, NA, NA),
       lty = c(NA, NA, seq_along(ns)), lwd = 1, bg = "white")
dev.off()


# --- Plot 5: Null p-value QQ vs Uniform(0,1) ---------------------------------
# On the diagonal => uniform => calibrated. One coloured curve per variant.
qq_cols <- c("black", "red", "blue", "darkgreen")
pdf("plots/plot_pvalue_qq.pdf", width = 6, height = 5)
plot(NA, axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
     xlab = "Theoretical uniform quantile", ylab = "Empirical p-value quantile")
axis(1); axis(2)
abline(0, 1, lty = 2)
nm_null <- names(pval_null)
for (k in seq_along(nm_null)) {
  p <- sort(pval_null[[k]][is.finite(pval_null[[k]])])
  lines(ppoints(length(p)), p, col = qq_cols[k], lwd = 2)
}
legend("bottomright", nm_null, col = qq_cols[seq_along(nm_null)],
       lwd = 2, bg = "white")
title(main = "Null p-values vs. Uniform(0,1)", font.main = 1)
dev.off()


# --- Plot 5b: Null p-value histograms (flat under the null) ------------------
pdf("plots/plot_pvalue_hist.pdf", width = 10, height = 3.2)
op <- par(mfrow = c(1, length(pval_null)), mar = c(4, 4, 3, 1))
for (k in seq_along(pval_null)) {
  p <- pval_null[[k]][is.finite(pval_null[[k]])]
  hist(p, breaks = seq(0, 1, 0.05), col = "grey70", border = "white",
       main = names(pval_null)[k], font.main = 1,
       xlab = "p-value", ylab = "Count")
  abline(h = length(p) / 20, lty = 2)
}
par(op)
dev.off()


# --- Plot 6: Continuous viewing -- per-sample vs one-per-period vs formula ----
# Grouped bars: 4 methods within each scenario x feature design.
meth_cols <- c("firebrick", "orange", "steelblue", "gray40")
bar_mat <- t(as.matrix(cv_results[, c("power_persample", "power_perperiod",
                                      "power_analytic_effcount",
                                      "power_analytic_periods")]))
rownames(bar_mat) <- c("Sim: per-sample", "Sim: one-per-period",
                       "Analytic: m_class=eff", "Analytic: m_class=#periods")
colnames(bar_mat) <- paste0(cv_results$scenario, "\n", cv_results$feature)
pdf("plots/plot_continuous_mclass.pdf", width = 9, height = 6)
op <- par(mar = c(7, 4, 3, 1))
bp <- barplot(bar_mat, beside = TRUE, col = meth_cols, ylim = c(0, 1),
              ylab = "Secondary-analysis power", las = 2, cex.names = 0.7,
              main = "Continuous viewing: does m_class = #periods hold?",
              font.main = 1)
abline(h = 0.8, lty = 2)
legend("bottomright", rownames(bar_mat), fill = meth_cols, bg = "white", cex = 0.8)
par(op)
dev.off()


# --- Plot 7: Stress-test QQ -- non-Gaussian noise, low vs high count ---------
stress_cols <- c("black", "red", "blue", "darkgreen")
pdf("plots/plot_stress_qq.pdf", width = 6, height = 5)
plot(NA, axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
     xlab = "Theoretical uniform quantile", ylab = "Empirical p-value quantile")
axis(1); axis(2); abline(0, 1, lty = 2)
nm_stress <- names(pval_stress)
for (k in seq_along(nm_stress)) {
  p <- sort(pval_stress[[k]][is.finite(pval_stress[[k]])])
  lines(ppoints(length(p)), p, col = stress_cols[k],
        lwd = 2, lty = if (k == length(nm_stress)) 2 else 1)
}
legend("bottomright", nm_stress, col = stress_cols[seq_along(nm_stress)],
       lwd = 2, lty = c(rep(1, length(nm_stress) - 1), 2), bg = "white")
title(main = "Null p-values under non-Gaussian noise", font.main = 1)
dev.off()


# --- Plot 7b: Stress-test histograms ----------------------------------------
pdf("plots/plot_stress_hist.pdf", width = 12, height = 3.2)
op <- par(mfrow = c(1, length(pval_stress)), mar = c(4, 4, 3, 1))
for (k in seq_along(pval_stress)) {
  p <- pval_stress[[k]][is.finite(pval_stress[[k]])]
  hist(p, breaks = seq(0, 1, 0.05), col = "grey70", border = "white",
       main = names(pval_stress)[k], font.main = 1,
       xlab = "p-value", ylab = "Count")
  abline(h = length(p) / 20, lty = 2)
}
par(op)
dev.off()


# ─────────────────────────────────────────────────────────────────────────────
# 10. Summary tables
# ─────────────────────────────────────────────────────────────────────────────

cat("\n\n=== Summary: Analytic vs Simulation (selected) ===\n")
summary_main <- results %>%
  filter(acc %in% c(0.6, 0.7, 0.8, 0.9)) %>%
  mutate(abs_diff = abs(power_analytic - power_sim)) 
print(summary_main, digits = 3)

cat("\n=== Summary: Assumption violations ===\n")
print(results_viol, digits = 3)
cat("\n=== Summary: Continuous viewing (per-sample vs formula) ===\n")
print(cv_results, digits = 3)

cat("Done.\n")