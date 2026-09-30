# Dunnett's minimum significant difference ------------------------------------
#
# Shared by compare_nsec_lecx.qmd and case_studies/ce_nsec_analysis.R. The
# percent minimum significant difference (PMSD) of a one-sided Dunnett's test
# is the effect size x of the first LECx method of de Bruyn and Elphick (2013),
# and sets delta for the NSEC-MSD.

# One-sided critical value of Dunnett's many-to-one test. The comparisons share
# the control mean, so their correlations factor as rho_ij = l_i * l_j, with
# l_i = sqrt(n_i / (n_i + n_control)). Conditional on one standard normal
# variable the comparisons are then independent, and the probability that all
# fall below d reduces to a two-dimensional integral, over that variable and
# the scaled chi distribution of the pooled SD. It is evaluated by Simpson's
# rule on a fixed 401 x 401 grid, which is deterministic. mvtnorm::qmvt() was
# rejected because its randomised integration returned values differing in the
# fourth decimal between seeds. The grid value agrees to six decimals with
# adaptive integrate() and with an 801-point grid, lies within 0.0005 of the
# mean of ten qmvt() runs, equals Student's t for a single treatment, and gives
# 2.356 for four treatments with 15 df (tabulated as 2.36).
dunnett_crit <- function(n_control, n_treat, df, alpha, n_nodes = 401) {
  l <- sqrt(n_treat / (n_treat + n_control))
  simpson <- function(a, b) {
    w <- rep(c(2, 4), length.out = n_nodes)
    w[c(1, n_nodes)] <- 1
    list(x = seq(a, b, length.out = n_nodes), w = w * (b - a) / (n_nodes - 1) / 3)
  }
  u <- simpson(-8, 8)
  s <- simpson(sqrt(qchisq(1e-10, df) / df), sqrt(qchisq(1 - 1e-10, df) / df))
  dens_s <- exp(log(2) + (df / 2) * log(df / 2) - lgamma(df / 2) + (df - 1) * log(s$x) - df * s$x^2 / 2)
  weights <- outer(u$w * dnorm(u$x), s$w * dens_s)
  p_all_below <- function(d) {
    p <- 1
    for (li in l) p <- p * pnorm(outer(-li * u$x, d * s$x, "+") / sqrt(1 - li^2))
    sum(weights * p)
  }
  uniroot(\(d) p_all_below(d) - (1 - alpha), c(0.9 * qt(1 - alpha, df), 10), tol = 1e-10)$root
}

# MSD for a decrease from the control, from a one-way ANOVA on the raw
# response with every tested concentration as a group. The PMSD is expressed
# relative to the observed control mean, as in USEPA (2000). The treatment
# replication must be equal for a single MSD to apply to every comparison;
# unequal replication stops the analysis rather than choosing a convention.
dunnett_msd <- function(conc, response, alpha, control = min(conc)) {
  groups <- split(response, conc)
  n <- lengths(groups)
  is_control <- names(n) == as.character(control)
  n_control <- n[is_control][[1]]
  n_treat <- n[!is_control]
  if (length(unique(n_treat)) != 1) {
    stop("Unequal treatment replication; choose an MSD convention.")
  }
  df <- sum(n) - length(n)
  s_pooled <- sqrt(sum(vapply(groups, \(g) sum((g - mean(g))^2), numeric(1))) / df)
  d_crit <- dunnett_crit(n_control, n_treat, df, alpha)
  msd <- d_crit * s_pooled * sqrt(1 / n_control + 1 / n_treat[[1]])
  control_mean <- mean(groups[is_control][[1]])
  list(
    k = length(n_treat), n_control = n_control, n_treatment = n_treat[[1]],
    df = df, s_pooled = s_pooled, d_crit = d_crit, msd = msd,
    control_mean = control_mean, pmsd = 100 * msd / control_mean
  )
}
