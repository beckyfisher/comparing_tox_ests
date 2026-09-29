# NSEC and NSEC-cv for the cerium (chloride) case study --------------------
#
# Applies the NSEC and the control-variance anchored NSEC (NSEC-cv), as defined
# in compare_nsec_lecx.qmd, to the ten single-species cerium chloride toxicity
# tests supplied by J. Baker (workbook received 2026-09-09). The LECx, NOEC,
# LOEC, IC10 and IC25 for the same tests are being computed separately in CETIS,
# and the results table written here is laid out for direct comparison with
# that output.
#
# Run from anywhere inside the project:
#   Rscript case_studies/ce_nsec_analysis.R
#
# Outputs (case_studies/output/):
#   ce_data_long.csv     the imported data, one row per replicate
#   ce_nsec_results.csv  one row per test: fitted parameters, NSEC, NSEC-cv,
#                        their equivalent effect sizes, and EC10/EC25
#   ce_nsec_fits.png     data, both fitted curves, and both thresholds per test
#   ce_session_info.txt  R and package versions used for the run
#
# No step is stochastic, so no seed is set: repeated runs give identical output.

required <- c(
  "dplyr", "ggplot2", "here", "patchwork", "purrr", "readr", "readxl",
  "scales", "tibble", "tidyr"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Install these packages before running: ", paste(missing, collapse = ", "))
}

here::i_am("case_studies/ce_nsec_analysis.R")

data_file <- here::here(
  "case_studies", "Ce dataset - for statistical comparison.xlsx"
)
out_dir <- here::here("case_studies", "output")
dir.create(out_dir, showWarnings = FALSE)

# One-sided significance level for both thresholds, as in compare_nsec_lecx.qmd.
alpha <- 0.05


# Import -------------------------------------------------------------------

# The workbook holds ten tables on one sheet, laid out side by side and stacked
# at irregular offsets, with header wording that varies between tables
# ("Test Conc (ug/L)", "Conc-mg/L", "C-ug/L"). Each table is located from its
# concentration header cell rather than from hard-coded cell ranges, so that the
# import fails loudly if the layout changes and can be reused on the further
# workbooks expected in the same format. The anchors are:
#   - the concentration header: the only cells containing "ug/L" or "mg/L";
#   - replicate columns: contiguous cells to its right that begin "Rep";
#   - data rows: consecutive rows below it with a numeric concentration;
#   - the title ("14-d Hyalella azteca growth"): the nearest cell above it, in
#     the same column, that begins with a duration;
#   - the response label ("Dry weight (mg)"): the cell diagonally above the
#     first replicate column.

read_tox_table <- function(cells, header_row, header_col) {
  rep_cols <- integer(0)
  j <- header_col + 1
  while (j <= ncol(cells) && grepl("^Rep", cells[header_row, j])) {
    rep_cols <- c(rep_cols, j)
    j <- j + 1
  }

  data_rows <- integer(0)
  i <- header_row + 1
  while (i <= nrow(cells) && !is.na(suppressWarnings(as.numeric(cells[i, header_col])))) {
    data_rows <- c(data_rows, i)
    i <- i + 1
  }

  above <- rev(seq_len(header_row - 1))
  title_pattern <- "^(\\d+)\\s*-\\s*([dh])\\s+(\\S+\\s+\\S+)\\s+(.+)$"
  title_row <- above[grepl(title_pattern, cells[above, header_col])][1]
  title <- cells[title_row, header_col]
  parts <- regmatches(title, regexec(title_pattern, title))[[1]]

  if (length(rep_cols) == 0 || length(data_rows) == 0 || is.na(title)) {
    stop("Table at row ", header_row, ", column ", header_col, " could not be parsed.")
  }

  values <- cells[data_rows, rep_cols, drop = FALSE]
  tibble::tibble(
    test = gsub(" ", "_", tolower(parts[4])),
    species = parts[4],
    duration = paste0(parts[2], "-", parts[3]),
    endpoint = parts[5],
    response_label = cells[header_row - 1, rep_cols[1]],
    conc_unit = regmatches(cells[header_row, header_col], regexpr("(ug|mg)/L", cells[header_row, header_col])),
    conc = rep(as.numeric(cells[data_rows, header_col]), times = length(rep_cols)),
    rep = rep(seq_along(rep_cols), each = length(data_rows)),
    response = as.numeric(values),
    sheet_row = rep(data_rows, times = length(rep_cols)),
    sheet_col = rep(rep_cols, each = length(data_rows))
  ) |>
    # Unequal replication (e.g. seven control and four treatment replicates in
    # the algal tests) leaves empty cells, which are not observations.
    dplyr::filter(!is.na(response))
}

read_tox_workbook <- function(path, sheet = 1) {
  cells <- readxl::read_excel(
    path,
    sheet = sheet, col_names = FALSE, col_types = "text", .name_repair = "minimal"
  ) |>
    as.matrix()
  headers <- which(!is.na(cells) & grepl("\\b(ug|mg)/L\\b", cells), arr.ind = TRUE)
  headers <- headers[order(headers[, "row"], headers[, "col"]), , drop = FALSE]
  purrr::map2(headers[, "row"], headers[, "col"], \(r, c) read_tox_table(cells, r, c)) |>
    dplyr::bind_rows() |>
    dplyr::mutate(test = factor(test, levels = unique(test)))
}

ce <- read_tox_workbook(data_file)

design <- ce |>
  dplyr::group_by(test, duration, endpoint, response_label, conc_unit) |>
  dplyr::summarise(
    n_conc = dplyr::n_distinct(conc),
    n_control = sum(conc == 0),
    reps_per_treatment = paste(unique(table(conc[conc > 0])), collapse = "/"),
    conc_range = paste0(min(conc[conc > 0]), "-", max(conc)),
    control_mean = mean(response[conc == 0]),
    .groups = "drop"
  )

# Checks against the workbook as received. A failure here means the layout or
# content differs from the 2026-09-09 file and the import must be re-examined.
stopifnot(
  nrow(design) == 10,
  !anyDuplicated(design$test),
  all(design$n_control >= 3),
  all(design$n_conc >= 5),
  all(ce$response >= 0)
)

print(as.data.frame(design), digits = 4)
readr::write_csv(ce, file.path(out_dir, "ce_data_long.csv"))


# Model and thresholds -------------------------------------------------------

# The concentration-response model is the three-parameter exponential decay
# used for the NSEC and NSEC-cv in compare_nsec_lecx.qmd,
#   y = b0 * exp(-b1 * x^b2),
# fitted by ordinary least squares on the raw response scale for every test.
# This keeps the definitions identical to those already circulated. It does
# not model the variance falling with the mean, which applies to the neonate
# counts (Ceriodaphnia) and to tests where the response approaches zero.
#
# The model is fitted in the equivalent form y = b0 * exp(-(x / e)^b2), where
# b1 = e^-b2 and e is the concentration at which the response falls to
# exp(-1), or 37%, of b0. With x in ug/L up to 1760 and b2 near 3, b1 is of
# order 1e-10 and nls() is poorly conditioned in (b1, b2); e lies within the
# tested range and is estimated on the log scale. The fitted curve, SE(b0) and
# every threshold are identical in either form, and b1 is reported so that the
# parameters can be compared with the qmd.
decay <- function(x, b0, e, b2) b0 * exp(-(x / e)^b2)

# Concentration at which the curve reaches response y.
solve_conc <- function(y, b0, e, b2) e * (-log(y / b0))^(1 / b2)

# Upper bound on b2. A response that is flat and then falls between two
# adjacent concentrations has no finite least-squares b2: the residual sum of
# squares keeps falling as the curve approaches a vertical step. The bound
# keeps such a fit finite. At b2 = 50 the curve falls from 95% to 5% of b0
# over a 1.09-fold change in concentration, which is a step at the spacing of
# any dilution series here.
b2_upper <- 50

# Starting values come from the profile of the residual sum of squares over
# b2: for each b2 on a grid, the RSS is minimised over e, with b0 either fixed
# or at its exact least-squares value (b0 enters the model linearly). A single
# default start does not converge across response shapes that range from a
# gradual decline (Oncorhynchus) to a near step (Pimephales, Neocloeon), and a
# coarse two-dimensional grid was found to select a near-step start for
# Pimephales, where the gradient is flat, although its least-squares b2 is 5.6.
#
# The profile also decides whether b2 is identified. If its minimum lies at
# b2_upper, the RSS falls monotonically towards a step and nls() is not run
# with b2 free, since it cannot converge on a flat surface.
profile_b2 <- function(x, y, b0_fixed = NULL) {
  pos <- x[x > 0]
  log_e_grid <- seq(log(min(pos) / 10), log(max(pos) * 10), length.out = 200)
  b2_grid <- 10^seq(-1, log10(b2_upper), length.out = 60)
  rss_at <- function(log_e, b2) {
    f <- exp(-(x / exp(log_e))^b2)
    b0 <- if (is.null(b0_fixed)) sum(y * f) / sum(f^2) else b0_fixed
    c(rss = sum((y - b0 * f)^2), b0 = b0)
  }
  prof <- purrr::map(b2_grid, \(b2) {
    rss <- vapply(log_e_grid, \(le) rss_at(le, b2)[["rss"]], numeric(1))
    k <- which.min(rss)
    bracket <- log_e_grid[c(max(1, k - 1), min(length(log_e_grid), k + 1))]
    opt <- optimize(\(le) rss_at(le, b2)[["rss"]], bracket, tol = 1e-10)
    list(b0 = rss_at(opt$minimum, b2)[["b0"]], log_e = opt$minimum, b2 = b2, rss = opt$objective)
  })
  rss <- vapply(prof, \(p) p$rss, numeric(1))
  best <- prof[[which.min(rss)]]
  # A relative tolerance, because along a monotone profile the RSS at large b2
  # is equal to machine precision and which.min() may select any of them.
  best$at_bound <- rss[length(rss)] <= min(rss) * (1 + 1e-6)
  best
}

# Fits the model with b0 estimated (b0_fixed = NULL) or fixed. The port
# algorithm keeps b2 within (0.05, b2_upper); b2 > 0 is required for a monotone
# decline and for the curve to equal b0 at the control.
#
# Where b2 is not identified, the model is fitted with b2 fixed at b2_upper and
# flagged. The fitted curve is then a step, and any threshold defined against
# the control falls just below the concentration at which the step occurs,
# whatever effect size it is anchored on. Its exact value depends on b2_upper,
# and it should be read as lying between the two tested concentrations that
# bracket the step.
fit_decay <- function(d, b0_fixed = NULL) {
  start <- profile_b2(d$conc, d$response, b0_fixed)
  pars <- c(if (is.null(b0_fixed)) "b0", "log_e", if (!start$at_bound) "b2")
  if (!is.null(b0_fixed)) d$b0 <- b0_fixed
  if (start$at_bound) d$b2 <- b2_upper

  fit <- nls(
    response ~ decay(conc, b0, exp(log_e), b2),
    data = d, start = start[pars], algorithm = "port",
    lower = c(b0 = 0, log_e = -Inf, b2 = 0.05)[pars],
    upper = c(b0 = Inf, log_e = Inf, b2 = b2_upper)[pars],
    control = list(maxiter = 1000, eval.max = 2000)
  )
  cf <- coef(fit)
  b2 <- if (start$at_bound) b2_upper else cf[["b2"]]
  list(
    fit = fit,
    b0 = if (is.null(b0_fixed)) cf[["b0"]] else b0_fixed,
    e = exp(cf[["log_e"]]),
    b2 = b2,
    b2_fixed = start$at_bound || b2 >= b2_upper * 0.999
  )
}

analyse_test <- function(d) {
  ctrl <- d$response[d$conc == 0]

  # NSEC, b0 estimated: the threshold is the lower one-sided (1 - alpha)
  # confidence bound on b0, using SE(b0) from the fit and the residual df.
  m <- fit_decay(d)
  se_b0 <- sqrt(vcov(m$fit)["b0", "b0"])
  t_crit <- qt(1 - alpha, df.residual(m$fit))
  lower <- m$b0 - t_crit * se_b0
  nsec <- solve_conc(lower, m$b0, m$e, m$b2)

  # NSEC-cv: b0 fixed at the observed control mean, and the threshold is the
  # lower one-sided (1 - alpha) confidence bound on that mean, from the control
  # replicates alone (df = n_control - 1).
  #
  # This is a one-sample bound on the control mean. The LECx effect size (the
  # PMSD) is a two-sample minimum significant difference instead: it includes a
  # factor of sqrt(1/n_control + 1/n_treatment) rather than sqrt(1/n_control),
  # uses the variance pooled over all groups, and in CETIS uses the Dunnett
  # critical value rather than t. The components are written to the results so
  # the two effect sizes can be reconciled term by term against the CETIS
  # output.
  ctrl_mean <- mean(ctrl)
  ctrl_sd <- sd(ctrl)
  n_ctrl <- length(ctrl)
  t_crit_cv <- qt(1 - alpha, n_ctrl - 1)
  lower_cv <- ctrl_mean - t_crit_cv * ctrl_sd / sqrt(n_ctrl)
  m_cv <- fit_decay(d, b0_fixed = ctrl_mean)
  nsec_cv <- solve_conc(lower_cv, ctrl_mean, m_cv$e, m_cv$b2)

  # EC10 and EC25 from the b0-estimated curve. These are not NSEC quantities;
  # they are reported so that this curve can be checked against the CETIS
  # IC10 and IC25 before any difference between NSEC-cv and LECx is attributed
  # to the threshold definition rather than to the fitted curve.
  ecx <- \(p) solve_conc((1 - p) * m$b0, m$b0, m$e, m$b2)

  # A threshold outside the tested range is set by the shape of the model
  # where there are no data. Below the lowest treatment this depends mainly on
  # b2: with b2 < 1 the curve has infinite slope at zero concentration, and the
  # threshold can fall far below any tested concentration.
  x_min <- min(d$conc[d$conc > 0])
  x_max <- max(d$conc)
  notes <- c(
    if (m$b2_fixed) "b2 not identified with b0 estimated; fitted at bound",
    if (m_cv$b2_fixed) "b2 not identified with b0 fixed; fitted at bound",
    if (nsec < x_min) "NSEC below lowest tested concentration",
    if (nsec_cv < x_min) "NSEC-cv below lowest tested concentration",
    if (nsec > x_max) "NSEC above highest tested concentration",
    if (nsec_cv > x_max) "NSEC-cv above highest tested concentration"
  )

  tibble::tibble(
    conc_min = x_min,
    conc_max = x_max,
    n_obs = nrow(d),
    b0 = m$b0,
    b1 = m$e^-m$b2,
    b2 = m$b2,
    e = m$e,
    se_b0 = se_b0,
    resid_sd = summary(m$fit)$sigma,
    df_resid = df.residual(m$fit),
    t_crit = t_crit,
    nsec_threshold = lower,
    nsec = nsec,
    nsec_ecx_pct = 100 * (1 - lower / m$b0),
    control_mean = ctrl_mean,
    control_sd = ctrl_sd,
    n_control = n_ctrl,
    t_crit_cv = t_crit_cv,
    nsec_cv_threshold = lower_cv,
    b1_cv = m_cv$e^-m_cv$b2,
    b2_cv = m_cv$b2,
    e_cv = m_cv$e,
    resid_sd_cv = summary(m_cv$fit)$sigma,
    nsec_cv = nsec_cv,
    nsec_cv_ecx_pct = 100 * (1 - lower_cv / ctrl_mean),
    ec10 = ecx(0.10),
    ec25 = ecx(0.25),
    fit_note = paste(notes, collapse = "; ")
  )
}

meta <- ce |>
  dplyr::distinct(test, species, duration, endpoint, response_label, conc_unit)

results <- ce |>
  dplyr::group_by(test) |>
  dplyr::group_modify(\(d, key) analyse_test(d)) |>
  dplyr::ungroup() |>
  dplyr::left_join(meta, by = "test") |>
  dplyr::relocate(species, duration, endpoint, response_label, conc_unit, .after = test)

readr::write_csv(results, file.path(out_dir, "ce_nsec_results.csv"))

results |>
  dplyr::transmute(
    species, conc_unit,
    nsec = signif(nsec, 3), nsec_ecx_pct = round(nsec_ecx_pct, 1),
    nsec_cv = signif(nsec_cv, 3), nsec_cv_ecx_pct = round(nsec_cv_ecx_pct, 1),
    ec10 = signif(ec10, 3), ec25 = signif(ec25, 3),
    fit_note
  ) |>
  as.data.frame() |>
  print()


# Figure ---------------------------------------------------------------------

# Concentrations are shown on a log axis, which cannot place the control at
# zero. The control is drawn at half the smallest positive value in the panel
# (lowest treatment or threshold) and labelled "0"; the curves are evaluated at
# true concentrations throughout.
fit_colours <- c(nsec = "#2a78d6", nsec_cv = "#eb6834")
fit_labels <- c(
  nsec = "NSEC: b0 estimated",
  nsec_cv = "NSEC-cv: b0 fixed at control mean"
)

format_conc <- function(x) {
  vapply(x, \(v) format(signif(v, 3), scientific = FALSE, drop0trailing = TRUE, trim = TRUE), "")
}

plot_test <- function(d, r) {
  thresholds <- c(r$nsec, r$nsec_cv)
  thresholds <- thresholds[is.finite(thresholds)]
  ctrl_pos <- min(c(r$conc_min, thresholds)) / 2
  x_hi <- max(c(r$conc_max, thresholds))

  d <- dplyr::mutate(d, plot_conc = ifelse(conc == 0, ctrl_pos, conc))
  grid <- 10^seq(log10(ctrl_pos), log10(x_hi), length.out = 300)
  curves <- dplyr::bind_rows(
    tibble::tibble(fit = "nsec", conc = grid, response = decay(grid, r$b0, r$e, r$b2)),
    tibble::tibble(fit = "nsec_cv", conc = grid, response = decay(grid, r$control_mean, r$e_cv, r$b2_cv))
  )
  marks <- tibble::tibble(
    fit = c("nsec", "nsec_cv"),
    conc = c(r$nsec, r$nsec_cv),
    response = c(r$nsec_threshold, r$nsec_cv_threshold)
  ) |>
    dplyr::filter(is.finite(conc))

  log_breaks <- scales::breaks_log(n = 5)(c(r$conc_min, x_hi))
  log_breaks <- log_breaks[log_breaks > ctrl_pos * 1.5 & log_breaks <= x_hi * 1.05]
  unit_label <- sub("ug", "\u00b5g", r$conc_unit)

  ggplot2::ggplot(d, ggplot2::aes(plot_conc, response)) +
    ggplot2::geom_point(colour = "#52514e", alpha = 0.6, size = 1.4, shape = 16) +
    ggplot2::geom_line(
      data = curves,
      ggplot2::aes(conc, response, colour = fit, linetype = fit),
      linewidth = 0.7
    ) +
    ggplot2::geom_segment(
      data = marks,
      ggplot2::aes(x = conc, xend = conc, y = -Inf, yend = response, colour = fit),
      linetype = "dotted", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = marks,
      ggplot2::aes(conc, response, colour = fit),
      size = 2.2, shape = 21, fill = "white", stroke = 1
    ) +
    ggplot2::scale_x_log10(
      breaks = c(ctrl_pos, log_breaks),
      labels = c("0", format_conc(log_breaks))
    ) +
    ggplot2::scale_colour_manual(values = fit_colours, labels = fit_labels, name = NULL) +
    ggplot2::scale_linetype_manual(
      values = c(nsec = "solid", nsec_cv = "longdash"), labels = fit_labels, name = NULL
    ) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      title = bquote(italic(.(r$species)) ~ "(" * .(r$duration) ~ .(r$endpoint) * ")"),
      x = paste0("Ce (", unit_label, ")"),
      y = r$response_label
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      axis.line = ggplot2::element_line(colour = "grey40", linewidth = 0.3),
      plot.title = ggplot2::element_text(size = 9.5),
      legend.position = "bottom"
    )
}

panels <- purrr::map(levels(ce$test), \(t) {
  plot_test(dplyr::filter(ce, test == t), dplyr::filter(results, test == t))
})

# The collected legend takes its position from the panel themes. The patchwork
# `&` operator and plot_annotation(theme = ) are not used because patchwork
# 1.3.0 does not accept ggplot2 4.0 theme objects in either.
figure <- patchwork::wrap_plots(panels, ncol = 2) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    caption = paste(
      "Points are replicates; the control is plotted at the left-most tick, labelled 0.",
      "Open circles mark each threshold: the concentration at which the curve",
      "reaches the lower one-sided 95% bound on the control response.",
      sep = "\n"
    )
  )

ggplot2::ggsave(
  file.path(out_dir, "ce_nsec_fits.png"), figure,
  width = 8.5, height = 13, dpi = 200, bg = "white"
)

writeLines(
  c(paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "", capture.output(sessionInfo())),
  file.path(out_dir, "ce_session_info.txt")
)
