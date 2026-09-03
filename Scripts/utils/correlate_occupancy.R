#!/usr/bin/env Rscript

# =============================================================================
# correlate_occupancy.R
#
# Quantifies the occupancy / enrichment correlation between CBX2/7 targets and
# H3K27me2 / EZH2 over pre-specified regions or genome-wide bins.
#
# Unlike delta vs delta analysis (which examines logFC changes upon EZH2i),
# this script analyzes the baseline abundance / occupancy itself (logCPM / aveLogCPM).
#
# Key Features:
#   1. No low-read filtering (retains regions where EZH2 is absent/low to test whether
#      CBX2 and H3K27me2 colocalize independently of EZH2).
#   2. Pairwise Scatterplots with both OLS (method="lm") and GAM (method="gam") fits.
#   3. Multiple Linear Regression (lm): Target_occ ~ K27me2_occ + EZH2_occ with 1D Spatial HAC CIs,
#      standardized partial beta ratios, and semi-partial correlation squared (sr^2) variance decomposition.
#   4. Multiple Generalized Additive Model (GAM): Target_occ ~ s(K27me2_occ) + s(EZH2_occ) to test
#      non-linear saturation, background thresholds, and compare deviance explained to linear OLS.
#   5. Steiger's Z-tests for dependent overlapping correlations.
#   6. EZH2 Occupancy Tertile Stratification (Low EZH2/Absent vs Moderate vs High EZH2) to directly
#      test the independence hypothesis inside EZH2-depleted/background regions.
#
# CLI arguments:
#   --bin_size             INT    Primary bin size (default: 10000)
#   --bin_size_sensitivity INT    Sensitivity bin size (default: 2000)
#   --out_dir              PATH   Output directory (default: Analysis_Data/occupancy_vs_occupancy)
#   --lfc_dir              PATH   Directory containing existing LFC RDS files (default: Analysis_Data/delta_vs_delta/k27me2_dmso_peaks)
#   --occ_type             STR    Occupancy metric to extract ('avelogcpm' or 'dmso_logcpm', default: 'avelogcpm')
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(showtext)
  library(ppcor)
  library(cocor)
  library(sandwich)
  library(mgcv)
})

# ---------------------------------------------------------------------------
# Font setup (Helvetica -> Arial -> DejaVuSans fallback)
# ---------------------------------------------------------------------------
find_font_path <- function(pattern) {
  res <- tryCatch(
    system(sprintf("fc-list : file family | grep -i '%s' | grep '\\.ttf' | head -1 | awk -F: '{print $1}' | xargs", pattern), intern = TRUE),
    error = function(e) ""
  )
  if (length(res) > 0 && nzchar(res[1]) && file.exists(res[1])) {
    return(res[1])
  }
  return("")
}

helv_path  <- find_font_path("Helvetica")
arial_path <- find_font_path("Arial")
fallback_path <- find_font_path("DejaVuSans.ttf")
if (!nzchar(fallback_path) && file.exists("$HOME/miniconda3/envs/chipseq/fonts/DejaVuSans.ttf")) {
  fallback_path <- "$HOME/miniconda3/envs/chipseq/fonts/DejaVuSans.ttf"
}
if (!nzchar(fallback_path) && file.exists("/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf")) {
  fallback_path <- "/usr/share/fonts/dejavu-sans-fonts/DejaVuSans.ttf"
}

if (nzchar(helv_path)) {
  font_add("Helvetica", regular = helv_path)
  PLOT_FONT <- "Helvetica"
} else if (nzchar(arial_path)) {
  font_add("Arial", regular = arial_path)
  PLOT_FONT <- "Arial"
} else if (nzchar(fallback_path)) {
  font_add("Helvetica", regular = fallback_path)
  font_add("Arial", regular = fallback_path)
  PLOT_FONT <- "Helvetica"
} else {
  PLOT_FONT <- "sans"
}
showtext_auto()

# ---------------------------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------------------------
args_cli <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] + 1 > length(args)) return(default)
  args[idx[1] + 1]
}

bin_size             <- as.integer(get_flag("--bin_size",             args_cli, 10000))
bin_size_sensitivity <- as.integer(get_flag("--bin_size_sensitivity", args_cli, 2000))
out_dir              <- get_flag("--out_dir",                         args_cli, file.path("Analysis_Data", "occupancy_vs_occupancy"))
lfc_dir              <- get_flag("--lfc_dir",                         args_cli, file.path("Analysis_Data", "delta_vs_delta", "k27me2_dmso_peaks"))
occ_type             <- tolower(get_flag("--occ_type",                args_cli, "avelogcpm"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== correlate_occupancy.R ===\n")
cat(sprintf("  Primary bin size    : %d bp\n", bin_size))
cat(sprintf("  Sensitivity bin size: %d bp\n", bin_size_sensitivity))
cat(sprintf("  Occupancy metric    : %s\n", occ_type))
cat(sprintf("  LFC / source dir    : %s\n", lfc_dir))
cat(sprintf("  Output dir          : %s\n\n", out_dir))

# ---------------------------------------------------------------------------
# Plotting theme
# ---------------------------------------------------------------------------
my_theme <- theme_bw(base_size = 13, base_family = PLOT_FONT) +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.2),
    strip.background = element_blank(),
    strip.text       = element_text(size = 15, face = "bold"),
    plot.title       = element_text(face = "bold", size = 15),
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    axis.text        = element_text(color = "black", size = 13),
    axis.title       = element_text(color = "black", size = 14),
    axis.line        = element_blank()
  )

save_svg <- function(plot_obj, filename, width = 6.5, height = 5.5) {
  svg(filename, width = width, height = height)
  print(plot_obj)
  invisible(dev.off())
  cat("Saved:", filename, "\n")
}

# ---------------------------------------------------------------------------
# Helper: 1D Spatial HAC (Conley / Newey-West) Standard Errors along Chromosomes
# ---------------------------------------------------------------------------
calc_spatial_hac_vcov <- function(mod, seqnames, pos, max_lag_dist = 500000) {
  X <- model.matrix(mod)
  w <- weights(mod)
  if (is.null(w)) w <- rep(1, nrow(X))
  e <- residuals(mod)
  
  N <- nrow(X)
  k <- ncol(X)
  u <- X * w * e
  
  chrs <- unique(seqnames)
  M_hac <- matrix(0, nrow = k, ncol = k)
  colnames(M_hac) <- colnames(X)
  rownames(M_hac) <- colnames(X)
  
  for (chr in chrs) {
    idx <- which(seqnames == chr)
    if (length(idx) <= 1) next
    
    ord <- order(pos[idx])
    idx <- idx[ord]
    
    u_chr <- u[idx, , drop = FALSE]
    p_chr <- pos[idx]
    n_c   <- length(idx)
    
    M_hac <- M_hac + crossprod(u_chr)
    
    for (l in 1:(n_c - 1)) {
      dists <- p_chr[(1 + l):n_c] - p_chr[1:(n_c - l)]
      valid_l <- which(dists <= max_lag_dist)
      if (length(valid_l) == 0) break
      
      wt_l <- 1 - dists[valid_l] / max_lag_dist
      u1 <- u_chr[valid_l, , drop = FALSE] * sqrt(wt_l)
      u2 <- u_chr[valid_l + l, , drop = FALSE] * sqrt(wt_l)
      M_hac <- M_hac + crossprod(u1, u2) + crossprod(u2, u1)
    }
  }
  
  B <- summary(mod)$cov.unscaled
  V_hac <- B %*% M_hac %*% B
  se_hac <- sqrt(pmax(0, diag(V_hac)))
  names(se_hac) <- names(coef(mod))
  
  t_crit <- qt(0.975, df = max(1, mod$df.residual))
  ci_low  <- coef(mod) - t_crit * se_hac
  ci_high <- coef(mod) + t_crit * se_hac
  
  list(se = se_hac, ci_low = ci_low, ci_high = ci_high, df = mod$df.residual, vcov = V_hac)
}

calc_block2mb_se <- function(mod, df_use) {
  if (!"block_2mb" %in% names(df_use)) {
    df_use$block_2mb <- paste0(df_use$seqnames, "_", floor(df_use$start / 2e6))
  }
  cluster <- as.factor(df_use$block_2mb)
  G <- nlevels(cluster)
  V <- sandwich::vcovCL(mod, cluster = cluster, type = "HC1")
  se <- sqrt(pmax(0, diag(V)))
  names(se) <- names(coef(mod))
  t_crit <- qt(0.975, df = max(1, G - 1))
  list(se = se, ci_low = coef(mod) - t_crit * se, ci_high = coef(mod) + t_crit * se)
}

gam_cv_dev_expl <- function(formula_str, data, k = 5) {
  set.seed(42)
  n <- nrow(data)
  if (n < 15) return(NA_real_)
  folds <- sample(rep(1:k, length.out = n))
  dev_res  <- numeric(n)
  dev_null <- numeric(n)
  y_var <- all.vars(as.formula(formula_str))[1]
  for (i in 1:k) {
    train_dat <- data[folds != i, , drop = FALSE]
    test_dat  <- data[folds == i, , drop = FALSE]
    fit <- tryCatch(gam(as.formula(formula_str), data = train_dat), error = function(e) NULL)
    if (is.null(fit)) next
    preds <- tryCatch(predict(fit, newdata = test_dat), error = function(e) rep(NA_real_, nrow(test_dat)))
    dev_res[folds == i]  <- (test_dat[[y_var]] - preds)^2
    dev_null[folds == i] <- (test_dat[[y_var]] - mean(train_dat[[y_var]], na.rm = TRUE))^2
  }
  valid <- !is.na(dev_res) & !is.na(dev_null)
  if (sum(valid) == 0 || sum(dev_null[valid]) == 0) return(NA_real_)
  1 - sum(dev_res[valid]) / sum(dev_null[valid])
}

compute_spatial_hac_summary <- function(mod, df_use, weighting_name = "Unweighted OLS") {
  seq_vec <- df_use$seqnames
  pos_vec <- df_use$start
  hac_res <- calc_spatial_hac_vcov(mod, seqnames = seq_vec, pos = pos_vec, max_lag_dist = 500000)
  blk_res <- calc_block2mb_se(mod, df_use)
  
  b   <- coef(mod)
  se  <- hac_res$se
  t   <- b / pmax(se, 1e-12)
  p   <- 2 * pt(-abs(t), df = mod$df.residual)
  
  sd_y <- sd(model.response(model.frame(mod)), na.rm = TRUE)
  X_mat <- model.matrix(mod)[, -1, drop = FALSE]
  sd_x <- apply(X_mat, 2, sd, na.rm = TRUE)
  std_b <- b[-1] * (sd_x / sd_y)
  
  if (ncol(X_mat) >= 2) {
    vif_vals <- tryCatch(diag(solve(cor(X_mat))), error = function(e) rep(NA_real_, ncol(X_mat)))
    names(vif_vals) <- colnames(X_mat)
    kappa_val <- tryCatch(kappa(X_mat, exact = TRUE), error = function(e) NA_real_)
  } else {
    vif_vals <- setNames(1.0, colnames(X_mat))
    kappa_val <- 1.0
  }
  
  list(
    coef = b, se = se, p = p, std_b = std_b,
    ci_low = hac_res$ci_low, ci_high = hac_res$ci_high,
    se_block = blk_res$se, ci_block_low = blk_res$ci_low, ci_block_high = blk_res$ci_high,
    vif = vif_vals, kappa = kappa_val,
    r_squared = summary(mod)$r.squared,
    adj_r_squared = summary(mod)$adj.r.squared
  )
}


# ---------------------------------------------------------------------------
# Core analysis function for a given bin size
# ---------------------------------------------------------------------------
run_analysis <- function(bs, label = "") {
  cat(sprintf("\n============================================================\n"))
  cat(sprintf("### OCCUPANCY ANALYSIS — %d-bp bins %s ###\n", bs, if (nzchar(label)) paste0("[", label, "]") else ""))
  cat(sprintf("============================================================\n"))
  
  rds_k27  <- file.path(lfc_dir, sprintf("lfc_bins_k27me2_%dbp.rds", bs))
  rds_cbx2 <- file.path(lfc_dir, sprintf("lfc_bins_cbx2_%dbp.rds", bs))
  rds_cbx7 <- file.path(lfc_dir, sprintf("lfc_bins_gstcbx7_%dbp.rds", bs))
  rds_ezh2 <- file.path(lfc_dir, sprintf("lfc_bins_ezh2_%dbp.rds", bs))
  rds_gst2 <- file.path(lfc_dir, sprintf("lfc_bins_gstcbx2_%dbp.rds", bs))
  rds_k27me3 <- file.path(lfc_dir, sprintf("lfc_bins_k27me3_%dbp.rds", bs))
  
  if (!file.exists(rds_k27) || !file.exists(rds_cbx2) || !file.exists(rds_cbx7) || !file.exists(rds_ezh2) || !file.exists(rds_gst2)) {
    cat("  NOTE: Required RDS files not found in lfc_dir for bin_size:", bs, "\n")
    return(invisible(NULL))
  }
  has_k27me3 <- file.exists(rds_k27me3)
  
  df_k27  <- readRDS(rds_k27)
  df_cbx2 <- readRDS(rds_cbx2)
  df_cbx7 <- readRDS(rds_cbx7)
  df_ezh2 <- readRDS(rds_ezh2)
  df_gst2 <- readRDS(rds_gst2)
  if (has_k27me3) df_k27me3 <- readRDS(rds_k27me3)
  
  common_ids <- intersect(intersect(intersect(intersect(
    df_k27$bin_id, df_cbx2$bin_id), df_cbx7$bin_id), df_ezh2$bin_id), df_gst2$bin_id)
  if (has_k27me3) common_ids <- intersect(common_ids, df_k27me3$bin_id)
  
  if (length(common_ids) < 50) {
    stop(sprintf("Too few common regions (%d) across all targets at %d bp", length(common_ids), bs))
  }
  
  cat(sprintf("Common regions across all targets at %d bp: N = %s\n", bs, format(length(common_ids), big.mark = ",")))
  
  df_k27  <- df_k27[match(common_ids, df_k27$bin_id), ]
  df_cbx2 <- df_cbx2[match(common_ids, df_cbx2$bin_id), ]
  df_cbx7 <- df_cbx7[match(common_ids, df_cbx7$bin_id), ]
  df_ezh2 <- df_ezh2[match(common_ids, df_ezh2$bin_id), ]
  df_gst2 <- df_gst2[match(common_ids, df_gst2$bin_id), ]
  if (has_k27me3) df_k27me3 <- df_k27me3[match(common_ids, df_k27me3$bin_id), ]
  
  # Extract occupancy metrics without low-read filtering
  df_merged <- data.frame(
    seqnames    = df_k27$seqnames,
    start       = df_k27$start,
    end         = df_k27$end,
    bin_id      = df_k27$bin_id,
    occ_k27me2  = df_k27$aveLogCPM,
    occ_ezh2    = df_ezh2$aveLogCPM,
    occ_cbx2    = df_cbx2$aveLogCPM,
    occ_gstcbx2 = df_gst2$aveLogCPM,
    occ_gstcbx7 = df_cbx7$aveLogCPM,
    se_k27me2   = df_k27$lfcSE,
    se_ezh2     = df_ezh2$lfcSE,
    se_cbx2     = df_cbx2$lfcSE,
    se_gstcbx2  = df_gst2$lfcSE,
    se_gstcbx7  = df_cbx7$lfcSE,
    stringsAsFactors = FALSE
  )
  if (has_k27me3) {
    df_merged$occ_k27me3 <- df_k27me3$aveLogCPM
    df_merged$se_k27me3  <- df_k27me3$lfcSE
  }
  df_merged$block_2mb <- paste0(df_merged$seqnames, "_", floor(df_merged$start / 2e6))
  
  # Remove only infinite or NA values if any exist
  keep_cols <- c("occ_k27me2", "occ_ezh2", "occ_cbx2", "occ_gstcbx2", "occ_gstcbx7", if (has_k27me3) "occ_k27me3")
  keep_finite <- complete.cases(df_merged[, keep_cols])
  for (col in keep_cols) keep_finite <- keep_finite & is.finite(df_merged[[col]])
  df_merged <- df_merged[keep_finite, ]
  N_use <- nrow(df_merged)
  
  cat(sprintf("Final finite regions analyzed (no aveLogCPM threshold filtering applied): N = %s\n\n", format(N_use, big.mark = ",")))
  
  cat("--- Occupancy Summary (aveLogCPM) ---\n")
  for (t_name in keep_cols) {
    cat(sprintf("  %12s : min = %6.2f, max = %6.2f, median = %6.2f, mean = %6.2f\n",
                t_name, min(df_merged[[t_name]]), max(df_merged[[t_name]]),
                median(df_merged[[t_name]]), mean(df_merged[[t_name]])))
  }
  
  sfx <- if (nzchar(label)) paste0("_", tolower(gsub(" ", "_", label))) else ""
  
  # ---------------------------------------------------------------------------
  # 1. Scatterplots with OLS + GAM overlays
  # ---------------------------------------------------------------------------
  make_occ_plot <- function(x_col, y_col, x_label, y_label, title_str) {
    rho   <- cor(df_merged[[x_col]], df_merged[[y_col]], method = "spearman")
    r_val <- cor(df_merged[[x_col]], df_merged[[y_col]], method = "pearson")
    
    # Fit GAM to get deviance explained
    gam_fit <- gam(as.formula(paste(y_col, "~ s(", x_col, ", bs = 'cs')")), data = df_merged)
    dev_exp <- summary(gam_fit)$dev.expl
    
    corr_label <- sprintf(
      "Spearman \u03c1 = %.3f\nGAM Dev.Expl = %.1f%%\nN = %s",
      rho, dev_exp * 100, format(N_use, big.mark = ",")
    )
    
    x_rng <- range(df_merged[[x_col]], finite = TRUE)
    y_rng <- range(df_merged[[y_col]], finite = TRUE)
    ann_x <- x_rng[1] + 0.05 * diff(x_rng)
    ann_y <- y_rng[2] - 0.05 * diff(y_rng)
    
    ggplot(df_merged, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_hex(bins = 60) +
      scale_fill_viridis_c(option = "plasma", name = "# regions") +
      # Linear OLS fit (dashed purple)
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#39107b", linewidth = 0.8, linetype = "dashed") +
      # Non-linear GAM fit (solid orange)
      geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "#d95f02", linewidth = 0.9) +
      annotate("text", x = ann_x, y = ann_y, label = corr_label, hjust = 0, vjust = 1, size = 4, family = PLOT_FONT, color = "black") +
      labs(x = x_label, y = y_label, title = title_str) +
      my_theme
  }
  
  save_svg(make_occ_plot("occ_k27me2", "occ_cbx2", "H3K27me2 Occupancy (aveLogCPM)", "WT CBX2 Occupancy (aveLogCPM)", "WT CBX2 vs H3K27me2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_cbx2_vs_k27me2_%dbp%s.svg", bs, sfx)))
  save_svg(make_occ_plot("occ_ezh2", "occ_cbx2", "EZH2 Occupancy (aveLogCPM)", "WT CBX2 Occupancy (aveLogCPM)", "WT CBX2 vs EZH2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_cbx2_vs_ezh2_%dbp%s.svg", bs, sfx)))
  
  save_svg(make_occ_plot("occ_k27me2", "occ_gstcbx2", "H3K27me2 Occupancy (aveLogCPM)", "GST-CBX2 Occupancy (aveLogCPM)", "GST-CBX2 vs H3K27me2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_gstcbx2_vs_k27me2_%dbp%s.svg", bs, sfx)))
  save_svg(make_occ_plot("occ_ezh2", "occ_gstcbx2", "EZH2 Occupancy (aveLogCPM)", "GST-CBX2 Occupancy (aveLogCPM)", "GST-CBX2 vs EZH2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_gstcbx2_vs_ezh2_%dbp%s.svg", bs, sfx)))
  
  save_svg(make_occ_plot("occ_k27me2", "occ_gstcbx7", "H3K27me2 Occupancy (aveLogCPM)", "GST-CBX7 Occupancy (aveLogCPM)", "GST-CBX7 vs H3K27me2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_gstcbx7_vs_k27me2_%dbp%s.svg", bs, sfx)))
  save_svg(make_occ_plot("occ_ezh2", "occ_gstcbx7", "EZH2 Occupancy (aveLogCPM)", "GST-CBX7 Occupancy (aveLogCPM)", "GST-CBX7 vs EZH2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_gstcbx7_vs_ezh2_%dbp%s.svg", bs, sfx)))
  
  save_svg(make_occ_plot("occ_k27me2", "occ_ezh2", "H3K27me2 Occupancy (aveLogCPM)", "EZH2 Occupancy (aveLogCPM)", "EZH2 vs H3K27me2 Occupancy [OLS & GAM]"),
           file.path(out_dir, sprintf("occupancy_k27me2_vs_ezh2_%dbp%s.svg", bs, sfx)))
  
  # ---------------------------------------------------------------------------
  # 2. Partial Correlations & Multiple Linear Regression vs GAM
  # ---------------------------------------------------------------------------
  cat("\n--- Marginal & Partial Occupancy Correlations (controlling for EZH2 Occupancy) ---\n")
  for (t_pair in list(c("occ_cbx2", "WT CBX2"), c("occ_gstcbx2", "GST-CBX2"), c("occ_gstcbx7", "GST-CBX7"))) {
    t_col  <- t_pair[1]
    t_name <- t_pair[2]
    
    # Partial correlation
    pcor_res <- pcor.test(df_merged[[t_col]], df_merged$occ_k27me2, df_merged$occ_ezh2, method = "spearman")
    cat(sprintf("  %10s ~ K27me2 | EZH2 : rho_partial = %+6.4f  p = %.2e\n", t_name, pcor_res$estimate, pcor_res$p.value))
  }
  
  cat("\n--- Multiple Linear Regression vs GAM Comparison ---\n")
  indep_metrics <- data.frame()
  decomp_long   <- data.frame()
  
  pred_cols  <- if (has_k27me3) c("occ_k27me2", "occ_ezh2", "occ_k27me3") else c("occ_k27me2", "occ_ezh2")
  formula_rhs <- paste(pred_cols, collapse = " + ")
  gam_rhs     <- paste(sprintf("s(%s, bs='cs')", pred_cols), collapse = " + ")
  
  for (t_pair in list(c("occ_cbx2", "WT CBX2", "se_cbx2"), c("occ_gstcbx2", "GST-CBX2", "se_gstcbx2"), c("occ_gstcbx7", "GST-CBX7", "se_gstcbx7"))) {
    t_col  <- t_pair[1]
    t_name <- t_pair[2]
    t_se   <- t_pair[3]
    
    # 1. Unweighted OLS
    mod_ols <- lm(as.formula(paste(t_col, "~", formula_rhs)), data = df_merged)
    res_ols <- compute_spatial_hac_summary(mod_ols, df_merged, "Unweighted OLS")
    
    # Semi-partial correlation squared (sr^2)
    r2_full <- res_ols$r_squared
    mod_no_k27 <- lm(as.formula(paste(t_col, "~", paste(setdiff(pred_cols, "occ_k27me2"), collapse = " + "))), data = df_merged)
    mod_no_ezh <- lm(as.formula(paste(t_col, "~", paste(setdiff(pred_cols, "occ_ezh2"), collapse = " + "))), data = df_merged)
    sr2_k27_unique <- max(0, r2_full - summary(mod_no_k27)$r.squared)
    sr2_ezh_unique <- max(0, r2_full - summary(mod_no_ezh)$r.squared)
    
    if (has_k27me3) {
      mod_no_k27m3 <- lm(as.formula(paste(t_col, "~", paste(setdiff(pred_cols, "occ_k27me3"), collapse = " + "))), data = df_merged)
      sr2_k27m3_unique <- max(0, r2_full - summary(mod_no_k27m3)$r.squared)
      sr2_shared <- max(0, r2_full - (sr2_k27_unique + sr2_ezh_unique + sr2_k27m3_unique))
    } else {
      sr2_k27m3_unique <- 0
      sr2_shared <- max(0, r2_full - (sr2_k27_unique + sr2_ezh_unique))
    }
    
    # GAM Comparison + 5-fold CV
    mod_gam <- gam(as.formula(paste(t_col, "~", gam_rhs)), data = df_merged)
    gam_summary <- summary(mod_gam)
    dev_exp_gam <- gam_summary$dev.expl
    cv_dev_exp  <- gam_cv_dev_expl(paste(t_col, "~", gam_rhs), df_merged, k = 5)
    edf_k27     <- gam_summary$s.table["s(occ_k27me2)", "edf"]
    edf_ezh     <- gam_summary$s.table["s(occ_ezh2)", "edf"]
    edf_k27m3   <- if (has_k27me3) gam_summary$s.table["s(occ_k27me3)", "edf"] else NA
    
    cat(sprintf("  [%s] OLS R² = %.4f | GAM Dev.Expl = %.4f (5-fold CV = %.4f) [edf: K27me2=%.1f, EZH2=%.1f%s]\n",
                t_name, r2_full, dev_exp_gam, cv_dev_exp, edf_k27, edf_ezh, if (has_k27me3) sprintf(", K27me3=%.1f", edf_k27m3) else ""))
    cat(sprintf("    OLS Betas : K27me2 β = %+6.4f (Std β* = %+6.4f, SpatialHAC p = %.2e, Block2Mb p = %.2e)\n",
                res_ols$coef["occ_k27me2"], res_ols$std_b["occ_k27me2"], res_ols$p["occ_k27me2"],
                2 * pt(-abs(res_ols$coef["occ_k27me2"] / pmax(1e-12, res_ols$se_block["occ_k27me2"])), df = mod_ols$df.residual)))
    if (has_k27me3) {
      cat(sprintf("                K27me3 β = %+6.4f (Std β* = %+6.4f, SpatialHAC p = %.2e, Block2Mb p = %.2e)\n",
                  res_ols$coef["occ_k27me3"], res_ols$std_b["occ_k27me3"], res_ols$p["occ_k27me3"],
                  2 * pt(-abs(res_ols$coef["occ_k27me3"] / pmax(1e-12, res_ols$se_block["occ_k27me3"])), df = mod_ols$df.residual)))
    }
    cat(sprintf("    Collinearity: K27me2 VIF = %.2f, EZH2 VIF = %.2f%s | Condition Number \u03ba = %.1f %s\n",
                res_ols$vif["occ_k27me2"], res_ols$vif["occ_ezh2"], if (has_k27me3) sprintf(", K27me3 VIF = %.2f", res_ols$vif["occ_k27me3"]) else "",
                res_ols$kappa, if (any(res_ols$vif > 5)) "[FLAG: VIF > 5 detected]" else "[No multicollinearity flag]"))
    cat(sprintf("    Variance  : Unique K27me2 sr² = %.4f (%.2f%%) | Unique EZH2 sr² = %.4f (%.2f%%)%s | Shared = %.4f\n",
                sr2_k27_unique, sr2_k27_unique * 100, sr2_ezh_unique, sr2_ezh_unique * 100,
                if (has_k27me3) sprintf(" | Unique K27me3 sr² = %.4f (%.2f%%)", sr2_k27m3_unique, sr2_k27m3_unique * 100) else "", sr2_shared))
    
    # Precision-Weighted (Response-only IVW)
    w_ivw <- 1 / (pmax(df_merged[[t_se]], 0.01)^2)
    mod_ivw <- lm(as.formula(paste(t_col, "~", formula_rhs)), data = df_merged, weights = w_ivw)
    res_ivw <- compute_spatial_hac_summary(mod_ivw, df_merged, "Response-only IVW")
    
    r2_ivw <- res_ivw$r_squared
    mod_no_k27_ivw <- lm(as.formula(paste(t_col, "~", paste(setdiff(pred_cols, "occ_k27me2"), collapse = " + "))), data = df_merged, weights = w_ivw)
    mod_no_ezh_ivw <- lm(as.formula(paste(t_col, "~", paste(setdiff(pred_cols, "occ_ezh2"), collapse = " + "))), data = df_merged, weights = w_ivw)
    sr2_k27_ivw <- max(0, r2_ivw - summary(mod_no_k27_ivw)$r.squared)
    sr2_ezh_ivw <- max(0, r2_ivw - summary(mod_no_ezh_ivw)$r.squared)
    if (has_k27me3) {
      mod_no_k27m3_ivw <- lm(as.formula(paste(t_col, "~", paste(setdiff(pred_cols, "occ_k27me3"), collapse = " + "))), data = df_merged, weights = w_ivw)
      sr2_k27m3_ivw <- max(0, r2_ivw - summary(mod_no_k27m3_ivw)$r.squared)
      sr2_shr_ivw <- max(0, r2_ivw - (sr2_k27_ivw + sr2_ezh_ivw + sr2_k27m3_ivw))
    } else {
      sr2_k27m3_ivw <- 0
      sr2_shr_ivw <- max(0, r2_ivw - (sr2_k27_ivw + sr2_ezh_ivw))
    }
    
    row_ols <- data.frame(
      Target = t_name, Weighting = "Unweighted OLS",
      Beta_K27me2 = res_ols$coef["occ_k27me2"], SE_SpatialHAC_K27 = res_ols$se["occ_k27me2"], SE_Block2Mb_K27 = res_ols$se_block["occ_k27me2"], p_SpatialHAC_K27 = res_ols$p["occ_k27me2"], Std_Beta_K27me2 = res_ols$std_b["occ_k27me2"],
      Beta_EZH2   = res_ols$coef["occ_ezh2"],   SE_SpatialHAC_EZH2 = res_ols$se["occ_ezh2"],   SE_Block2Mb_EZH2 = res_ols$se_block["occ_ezh2"],   p_SpatialHAC_EZH2 = res_ols$p["occ_ezh2"],   Std_Beta_EZH2   = res_ols$std_b["occ_ezh2"],
      VIF_K27me2  = res_ols$vif["occ_k27me2"],  VIF_EZH2           = res_ols$vif["occ_ezh2"],  Condition_Kappa = res_ols$kappa,
      OLS_R2 = r2_full, GAM_Dev_Expl = dev_exp_gam, GAM_CV_Dev_Expl = cv_dev_exp, GAM_edf_K27 = edf_k27, GAM_edf_EZH2 = edf_ezh,
      Unique_sr2_K27me2 = sr2_k27_unique, Unique_sr2_EZH2 = sr2_ezh_unique, Shared_sr2 = sr2_shared,
      stringsAsFactors = FALSE
    )
    if (has_k27me3) {
      row_ols$Beta_K27me3     <- res_ols$coef["occ_k27me3"]
      row_ols$SE_SpatialHAC_K27m3 <- res_ols$se["occ_k27me3"]
      row_ols$SE_Block2Mb_K27m3   <- res_ols$se_block["occ_k27me3"]
      row_ols$Std_Beta_K27me3 <- res_ols$std_b["occ_k27me3"]
      row_ols$VIF_K27me3      <- res_ols$vif["occ_k27me3"]
      row_ols$GAM_edf_K27me3  <- edf_k27m3
      row_ols$Unique_sr2_K27me3 <- sr2_k27m3_unique
    }
    
    row_ivw <- data.frame(
      Target = t_name, Weighting = "Response-only IVW",
      Beta_K27me2 = res_ivw$coef["occ_k27me2"], SE_SpatialHAC_K27 = res_ivw$se["occ_k27me2"], SE_Block2Mb_K27 = res_ivw$se_block["occ_k27me2"], p_SpatialHAC_K27 = res_ivw$p["occ_k27me2"], Std_Beta_K27me2 = res_ivw$std_b["occ_k27me2"],
      Beta_EZH2   = res_ivw$coef["occ_ezh2"],   SE_SpatialHAC_EZH2 = res_ivw$se["occ_ezh2"],   SE_Block2Mb_EZH2 = res_ivw$se_block["occ_ezh2"],   p_SpatialHAC_EZH2 = res_ivw$p["occ_ezh2"],   Std_Beta_EZH2   = res_ivw$std_b["occ_ezh2"],
      VIF_K27me2  = res_ivw$vif["occ_k27me2"],  VIF_EZH2           = res_ivw$vif["occ_ezh2"],  Condition_Kappa = res_ivw$kappa,
      OLS_R2 = r2_ivw, GAM_Dev_Expl = dev_exp_gam, GAM_CV_Dev_Expl = cv_dev_exp, GAM_edf_K27 = edf_k27, GAM_edf_EZH2 = edf_ezh,
      Unique_sr2_K27me2 = sr2_k27_ivw, Unique_sr2_EZH2 = sr2_ezh_ivw, Shared_sr2 = sr2_shr_ivw,
      stringsAsFactors = FALSE
    )
    if (has_k27me3) {
      row_ivw$Beta_K27me3     <- res_ivw$coef["occ_k27me3"]
      row_ivw$SE_SpatialHAC_K27m3 <- res_ivw$se["occ_k27me3"]
      row_ivw$SE_Block2Mb_K27m3   <- res_ivw$se_block["occ_k27me3"]
      row_ivw$Std_Beta_K27me3 <- res_ivw$std_b["occ_k27me3"]
      row_ivw$VIF_K27me3      <- res_ivw$vif["occ_k27me3"]
      row_ivw$GAM_edf_K27me3  <- edf_k27m3
      row_ivw$Unique_sr2_K27me3 <- sr2_k27m3_ivw
    }
    
    if (nrow(indep_metrics) == 0) {
      indep_metrics <- rbind(row_ols, row_ivw)
    } else {
      # Ensure matching columns before rbind
      for (col in union(names(indep_metrics), names(row_ols))) {
        if (!col %in% names(indep_metrics)) indep_metrics[[col]] <- NA
        if (!col %in% names(row_ols)) row_ols[[col]] <- NA
        if (!col %in% names(row_ivw)) row_ivw[[col]] <- NA
      }
      indep_metrics <- rbind(indep_metrics, row_ols, row_ivw)
    }
    
    decomp_long <- rbind(decomp_long,
      data.frame(Target = t_name, Weighting = "Unweighted OLS", Component = "Unique H3K27me2 (sr²)", Variance_Explained = sr2_k27_unique, stringsAsFactors = FALSE),
      data.frame(Target = t_name, Weighting = "Unweighted OLS", Component = "Unique EZH2 (sr²)",     Variance_Explained = sr2_ezh_unique, stringsAsFactors = FALSE)
    )
    if (has_k27me3) {
      decomp_long <- rbind(decomp_long,
        data.frame(Target = t_name, Weighting = "Unweighted OLS", Component = "Unique H3K27me3 (sr²)", Variance_Explained = sr2_k27m3_unique, stringsAsFactors = FALSE)
      )
    }
    decomp_long <- rbind(decomp_long,
      data.frame(Target = t_name, Weighting = "Unweighted OLS", Component = "Shared / Overlap",      Variance_Explained = sr2_shared,     stringsAsFactors = FALSE),
      data.frame(Target = t_name, Weighting = "Response-only IVW", Component = "Unique H3K27me2 (sr²)", Variance_Explained = sr2_k27_ivw, stringsAsFactors = FALSE),
      data.frame(Target = t_name, Weighting = "Response-only IVW", Component = "Unique EZH2 (sr²)",     Variance_Explained = sr2_ezh_ivw, stringsAsFactors = FALSE)
    )
    if (has_k27me3) {
      decomp_long <- rbind(decomp_long,
        data.frame(Target = t_name, Weighting = "Response-only IVW", Component = "Unique H3K27me3 (sr²)", Variance_Explained = sr2_k27m3_ivw, stringsAsFactors = FALSE)
      )
    }
    decomp_long <- rbind(decomp_long,
      data.frame(Target = t_name, Weighting = "Response-only IVW", Component = "Shared / Overlap",      Variance_Explained = sr2_shr_ivw,     stringsAsFactors = FALSE)
    )
  }
  
  tsv_metrics_path <- file.path(out_dir, sprintf("occupancy_colocalization_independence_metrics_%dbp%s.tsv", bs, sfx))
  write.table(indep_metrics, tsv_metrics_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nSaved occupancy metrics table to:", tsv_metrics_path, "\n")
  
  # Plot variance decomposition
  comp_levels <- if (has_k27me3) c("Shared / Overlap", "Unique H3K27me3 (sr²)", "Unique EZH2 (sr²)", "Unique H3K27me2 (sr²)") else c("Shared / Overlap", "Unique EZH2 (sr²)", "Unique H3K27me2 (sr²)")
  decomp_long$Component <- factor(decomp_long$Component, levels = comp_levels)
  decomp_long$Target    <- factor(decomp_long$Target, levels = c("WT CBX2", "GST-CBX2", "GST-CBX7"))
  decomp_long$Weighting <- factor(decomp_long$Weighting, levels = c("Unweighted OLS", "Response-only IVW"))
  
  pal_cols <- c("Shared / Overlap" = "#d9d9d9", "Unique EZH2 (sr²)" = "#3182bd", "Unique H3K27me2 (sr²)" = "#e6550d")
  if (has_k27me3) pal_cols["Unique H3K27me3 (sr²)"] <- "#7570b3"
  
  p_decomp <- ggplot(decomp_long, aes(x = Target, y = Variance_Explained, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black", linewidth = 0.5) +
    facet_wrap(~ Weighting) +
    scale_fill_manual(values = pal_cols) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = "Target", y = "Proportion of Occupancy Variance Explained (sr²)", fill = "Variance Component",
         title = sprintf("Occupancy Variance Decomposition (%d-bp bins)", bs)) +
    my_theme
  
  save_svg(p_decomp, file.path(out_dir, sprintf("occupancy_variance_decomposition_sr2_%dbp%s.svg", bs, sfx)), width = 10, height = 5.5)
  
  # ---------------------------------------------------------------------------
  # 3. Steiger's Z-Tests for Dependent Overlapping Correlations
  # ---------------------------------------------------------------------------
  cat("\n--- Steiger's Z-Tests: Comparing Correlation with K27me2 Occupancy ---\n")
  r_wt_k27  <- cor(df_merged$occ_cbx2,    df_merged$occ_k27me2, method = "pearson")
  r_gst2_k2 <- cor(df_merged$occ_gstcbx2, df_merged$occ_k27me2, method = "pearson")
  r_gst7_k2 <- cor(df_merged$occ_gstcbx7, df_merged$occ_k27me2, method = "pearson")
  r_wt_gst2 <- cor(df_merged$occ_cbx2,    df_merged$occ_gstcbx2, method = "pearson")
  r_wt_gst7 <- cor(df_merged$occ_cbx2,    df_merged$occ_gstcbx7, method = "pearson")
  r_gst2_7  <- cor(df_merged$occ_gstcbx2, df_merged$occ_gstcbx7, method = "pearson")
  
  cat("  -> Comparison 1: WT CBX2 vs GST-CBX2 over H3K27me2 Occupancy\n")
  cc_1 <- cocor.dep.groups.overlap(r.jk = r_wt_k27, r.jh = r_gst2_k2, r.kh = r_wt_gst2, n = N_use)
  print(cc_1)
  
  cat("\n  -> Comparison 2: WT CBX2 vs GST-CBX7 over H3K27me2 Occupancy\n")
  cc_2 <- cocor.dep.groups.overlap(r.jk = r_wt_k27, r.jh = r_gst7_k2, r.kh = r_wt_gst7, n = N_use)
  print(cc_2)
  
  cat("\n  -> Comparison 3: GST-CBX2 vs GST-CBX7 over H3K27me2 Occupancy\n")
  cc_3 <- cocor.dep.groups.overlap(r.jk = r_gst2_k2, r.jh = r_gst7_k2, r.kh = r_gst2_7, n = N_use)
  print(cc_3)
  
  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # 4. Stratification by EZH2 Occupancy Tertiles (Testing Independence in Low EZH2 Regions)
  # ---------------------------------------------------------------------------
  cat("\n--- EZH2 Occupancy Tertile Stratification ---\n")
  q_ezh <- quantile(df_merged$occ_ezh2, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  df_merged$ezh2_tertile <- cut(df_merged$occ_ezh2, breaks = q_ezh, include.lowest = TRUE,
                                labels = c("Low EZH2 / Absent", "Moderate EZH2", "High EZH2"))
  
  strat_df <- data.frame()
  for (ter in levels(df_merged$ezh2_tertile)) {
    sub_df <- df_merged[df_merged$ezh2_tertile == ter, ]
    n_sub  <- nrow(sub_df)
    
    for (t_pair in list(c("occ_cbx2", "WT CBX2"), c("occ_gstcbx2", "GST-CBX2"), c("occ_gstcbx7", "GST-CBX7"))) {
      t_col  <- t_pair[1]
      t_name <- t_pair[2]
      
      # Fit OLS bivariate inside tertile
      m_fit <- lm(as.formula(paste(t_col, "~ occ_k27me2")), data = sub_df)
      hac_sub <- calc_spatial_hac_vcov(m_fit, seqnames = sub_df$seqnames, pos = sub_df$start, max_lag_dist = 500000)
      blk_sub <- calc_block2mb_se(m_fit, sub_df)
      
      b_val   <- coef(m_fit)["occ_k27me2"]
      sd_y_sub <- sd(sub_df[[t_col]], na.rm = TRUE)
      sd_x_sub <- sd(sub_df$occ_k27me2, na.rm = TRUE)
      std_b_sub <- b_val * (sd_x_sub / sd_y_sub)
      
      # Fit GAM inside tertile
      g_fit <- tryCatch(gam(as.formula(paste(t_col, "~ s(occ_k27me2, bs='cs')")), data = sub_df), error = function(e) NULL)
      gam_dev <- if (!is.null(g_fit)) summary(g_fit)$dev.expl else NA_real_
      cv_dev  <- gam_cv_dev_expl(paste(t_col, "~ s(occ_k27me2, bs='cs')"), sub_df, k = 5)
      
      strat_df <- rbind(strat_df, data.frame(
        ezh2_tertile = ter, Target = t_name, N = n_sub,
        Spearman_rho = cor(sub_df[[t_col]], sub_df$occ_k27me2, method = "spearman"),
        Beta = b_val, Std_Beta = std_b_sub,
        SE_SpatialHAC = hac_sub$se["occ_k27me2"], CI_SpatialHAC_low = hac_sub$ci_low["occ_k27me2"], CI_SpatialHAC_high = hac_sub$ci_high["occ_k27me2"],
        SE_Block2Mb   = blk_sub$se["occ_k27me2"], CI_Block2Mb_low   = blk_sub$ci_low["occ_k27me2"], CI_Block2Mb_high   = blk_sub$ci_high["occ_k27me2"],
        OLS_R2 = summary(m_fit)$r.squared, GAM_Dev_Expl = gam_dev, GAM_CV_Dev_Expl = cv_dev,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  print(strat_df)
  tsv_strat_path <- file.path(out_dir, sprintf("occupancy_ezh2_tertile_stratification_%dbp%s.tsv", bs, sfx))
  write.table(strat_df, tsv_strat_path, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ---------------------------------------------------------------------------
  # 5. Stratification by H3K27me3 Occupancy Tertiles (Item 6)
  # ---------------------------------------------------------------------------
  if (has_k27me3) {
    cat("\n--- H3K27me3 Occupancy Tertile Stratification ---\n")
    q_k27m3 <- quantile(df_merged$occ_k27me3, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
    df_merged$k27me3_tertile <- cut(df_merged$occ_k27me3, breaks = q_k27m3, include.lowest = TRUE,
                                    labels = c("Low H3K27me3 / Absent", "Moderate H3K27me3", "High H3K27me3"))
    
    strat_m3_df <- data.frame()
    for (ter in levels(df_merged$k27me3_tertile)) {
      sub_df <- df_merged[df_merged$k27me3_tertile == ter, ]
      n_sub  <- nrow(sub_df)
      
      for (t_pair in list(c("occ_cbx2", "WT CBX2"), c("occ_gstcbx2", "GST-CBX2"), c("occ_gstcbx7", "GST-CBX7"))) {
        t_col  <- t_pair[1]
        t_name <- t_pair[2]
        
        m_fit <- lm(as.formula(paste(t_col, "~ occ_k27me2")), data = sub_df)
        hac_sub <- calc_spatial_hac_vcov(m_fit, seqnames = sub_df$seqnames, pos = sub_df$start, max_lag_dist = 500000)
        blk_sub <- calc_block2mb_se(m_fit, sub_df)
        
        b_val   <- coef(m_fit)["occ_k27me2"]
        sd_y_sub <- sd(sub_df[[t_col]], na.rm = TRUE)
        sd_x_sub <- sd(sub_df$occ_k27me2, na.rm = TRUE)
        std_b_sub <- b_val * (sd_x_sub / sd_y_sub)
        
        g_fit <- tryCatch(gam(as.formula(paste(t_col, "~ s(occ_k27me2, bs='cs')")), data = sub_df), error = function(e) NULL)
        gam_dev <- if (!is.null(g_fit)) summary(g_fit)$dev.expl else NA_real_
        cv_dev  <- gam_cv_dev_expl(paste(t_col, "~ s(occ_k27me2, bs='cs')"), sub_df, k = 5)
        
        strat_m3_df <- rbind(strat_m3_df, data.frame(
          k27me3_tertile = ter, Target = t_name, N = n_sub,
          Spearman_rho = cor(sub_df[[t_col]], sub_df$occ_k27me2, method = "spearman"),
          Beta = b_val, Std_Beta = std_b_sub,
          SE_SpatialHAC = hac_sub$se["occ_k27me2"], CI_SpatialHAC_low = hac_sub$ci_low["occ_k27me2"], CI_SpatialHAC_high = hac_sub$ci_high["occ_k27me2"],
          SE_Block2Mb   = blk_sub$se["occ_k27me2"], CI_Block2Mb_low   = blk_sub$ci_low["occ_k27me2"], CI_Block2Mb_high   = blk_sub$ci_high["occ_k27me2"],
          OLS_R2 = summary(m_fit)$r.squared, GAM_Dev_Expl = gam_dev, GAM_CV_Dev_Expl = cv_dev,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    print(strat_m3_df)
    tsv_strat_m3_path <- file.path(out_dir, sprintf("occupancy_k27me3_tertile_stratification_%dbp%s.tsv", bs, sfx))
    write.table(strat_m3_df, tsv_strat_m3_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  make_strat_bar <- function(strat_data, tertile_col, title_str, filename) {
    long_strat <- rbind(
      data.frame(Tertile = strat_data[[tertile_col]], Target = strat_data$Target, Metric = "Spearman \u03c1", Value = strat_data$Spearman_rho, stringsAsFactors = FALSE),
      data.frame(Tertile = strat_data[[tertile_col]], Target = strat_data$Target, Metric = "GAM Dev. Expl.", Value = strat_data$GAM_Dev_Expl, stringsAsFactors = FALSE)
    )
    
    p <- ggplot(long_strat, aes(x = Tertile, y = Value, fill = Metric)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, color = "black") +
      facet_wrap(~ Target) +
      geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(width = 0.7), vjust = -0.4, size = 3.5, family = PLOT_FONT) +
      scale_fill_manual(values = c("Spearman \u03c1" = "#7570b3", "GAM Dev. Expl." = "#d95f02")) +
      scale_y_continuous(limits = c(-0.1, 1.05)) +
      labs(x = "Occupancy Tertile", y = "Correlation (\u03c1) / Deviance Explained",
           title = title_str) +
      my_theme + theme(axis.text.x = element_text(angle = 20, hjust = 1))
    save_svg(p, file.path(out_dir, sprintf("%s_%dbp%s.svg", filename, bs, sfx)), width = 11, height = 5.5)
  }
  
  make_strat_bar(strat_df, "ezh2_tertile", "WT CBX2, GST-CBX2 & GST-CBX7 vs K27me2 across EZH2 Tertiles", "occupancy_all_targets_ezh2_strat")
  if (has_k27me3) {
    make_strat_bar(strat_m3_df, "k27me3_tertile", "WT CBX2, GST-CBX2 & GST-CBX7 vs K27me2 across H3K27me3 Tertiles", "occupancy_all_targets_k27me3_strat")
  }
  
  cat(sprintf("\n=== Occupancy analysis for %d-bp complete ===\n", bs))
}

run_analysis(bin_size, label = "")

if (bin_size_sensitivity > 0) {
  rds_sens <- file.path(lfc_dir, sprintf("lfc_bins_k27me2_%dbp.rds", bin_size_sensitivity))
  if (file.exists(rds_sens)) {
    run_analysis(bin_size_sensitivity, label = "sensitivity")
  } else {
    cat(sprintf("\nNOTE: %d-bp sensitivity RDS not found at %s (skipping sensitivity run).\n", bin_size_sensitivity, rds_sens))
  }
}

cat("\n=== correlate_occupancy.R finished successfully ===\n")
