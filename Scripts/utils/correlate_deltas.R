#!/usr/bin/env Rscript

# =============================================================================
# correlate_deltas.R
#
# Delta-vs-Delta (ΔlogFC vs ΔlogFC) genome-wide correlation analysis.
# Requires pre-computed LFC bin RDS files from compute_lfc_bins.R (step_7i).
#
# Analyses performed
# ------------------
# Core plots (my_theme, base R svg):
#   A. ΔCBX2 vs ΔH3K27me2   — hexbin scatterplot + Spearman ρ, Pearson r
#   B. ΔCBX7 vs ΔH3K27me2   — hexbin scatterplot + Spearman ρ, Pearson r
#
# Extra 1: Partial correlation controlling for ΔEZH2 (ppcor::pcor.test)
#   - ΔCBX2 | ΔH3K27me2, holding ΔEZH2 constant
#   - ΔCBX7 | ΔH3K27me2, holding ΔEZH2 constant
#   + Weighted regression: lm(delta_CBX ~ delta_K27me2 + delta_EZH2,
#                              weights = 1/se_cbx^2)
#
# Extra 2: Steiger's Z test (cocor::cocor.dep.groups.overlap)
#   Tests whether r(ΔCBX2, ΔK27me2) ≠ r(ΔCBX7, ΔK27me2)
#
# Extra 3: ΔEZH2 tertile stratification
#   Spearman ρ for CBX2/CBX7 vs H3K27me2 within each EZH2 tertile
#
# Extra 4: Sensitivity check at 2kb bins
#   Repeats core + partial correlation at smaller spatial scale
# =============================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(ggplot2)
  library(dplyr)
  library(ppcor)
  library(cocor)
  library(showtext)
})

# ---------------------------------------------------------------------------
# Font setup (Helvetica / Arial fallback)
# ---------------------------------------------------------------------------
font_families_available <- font_families()
if (!"Helvetica" %in% font_families_available) {
  hv_candidates <- c(
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/type1/gsfonts/n019003l.pfb",
    "/usr/share/fonts/truetype/freefont/FreeSans.ttf"
  )
  hv_found <- hv_candidates[file.exists(hv_candidates)]
  if (length(hv_found) > 0) {
    font_add("Helvetica", hv_found[1])
    cat("Font mapped: Helvetica →", hv_found[1], "\n")
  } else {
    cat("NOTE: No Helvetica/Arial-like font file found; using ggplot2 default.\n")
  }
}
showtext_auto()

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args_cli <- commandArgs(trailingOnly = TRUE)

get_flag <- function(flag, default = NULL) {
  idx <- which(args_cli == flag)
  if (length(idx) == 0) return(default)
  if (idx[1] + 1 > length(args_cli)) return(default)
  args_cli[idx[1] + 1]
}

bin_size         <- as.integer(get_flag("--bin_size",              default = "10000"))
bin_size_sens    <- as.integer(get_flag("--bin_size_sensitivity",  default = "2000"))
out_dir          <- get_flag("--out_dir", default = "Analysis_Data/delta_vs_delta")
lfc_dir          <- get_flag("--lfc_dir", default = out_dir)
lfc_type         <- get_flag("--lfc_type", default = "shrunk")
lfc_col          <- if (lfc_type == "unshrunk") "logFC_raw" else "logFC_shrunk"

cat(sprintf("\n=== correlate_deltas.R ===\n"))
cat(sprintf("  Primary bin size    : %d bp\n", bin_size))
cat(sprintf("  Sensitivity bin size: %d bp\n", bin_size_sens))
cat(sprintf("  LFC type            : %s (%s)\n", lfc_type, lfc_col))
cat(sprintf("  LFC source dir      : %s\n", lfc_dir))
cat(sprintf("  Output dir          : %s\n\n", out_dir))

# ---------------------------------------------------------------------------
# my_theme (sourced from inhibitor_resistance_analysis.R conventions)
# ---------------------------------------------------------------------------
my_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid    = element_blank(),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 1.2),
    strip.background = element_blank(),
    strip.text    = element_text(size = 17, family = "Helvetica"),
    plot.title    = element_text(face = "bold", family = "Helvetica"),
    legend.position  = "right",
    legend.title  = element_text(size = 13, family = "Helvetica"),
    legend.text   = element_text(size = 12, family = "Helvetica"),
    axis.text     = element_text(color = "black", size = 13, family = "Helvetica"),
    axis.title    = element_text(color = "black", size = 15, family = "Helvetica"),
    axis.line     = element_blank()
  )

# ---------------------------------------------------------------------------
# Helper: load LFC bin RDS for a target at a given bin size
# ---------------------------------------------------------------------------
load_lfc <- function(target, bs, od = lfc_dir) {
  target_safe <- tolower(gsub("-", "_", target))
  rds <- file.path(od, sprintf("lfc_bins_%s_%dbp.rds", target_safe, bs))
  if (!file.exists(rds)) stop("LFC RDS not found: ", rds)
  cat("Loading:", rds, "\n")
  readRDS(rds)
}

# ---------------------------------------------------------------------------
# Helper: inner join across 4 targets on bin_id
# ---------------------------------------------------------------------------
merge_targets <- function(bs) {
  k27   <- load_lfc("K27me2",  bs) %>% dplyr::rename(lfc_k27me2  = .data[[lfc_col]], se_k27me2  = lfcSE, alc_k27me2  = aveLogCPM)
  cbx2  <- load_lfc("CBX2",    bs) %>% dplyr::rename(lfc_cbx2    = .data[[lfc_col]], se_cbx2    = lfcSE, alc_cbx2    = aveLogCPM)
  cbx7  <- load_lfc("GSTCBX7", bs) %>% dplyr::rename(lfc_cbx7    = .data[[lfc_col]], se_cbx7    = lfcSE, alc_cbx7    = aveLogCPM)
  ezh2  <- load_lfc("EZH2",    bs) %>% dplyr::rename(lfc_ezh2    = .data[[lfc_col]], se_ezh2    = lfcSE, alc_ezh2    = aveLogCPM)
  gcbx2 <- load_lfc("GSTCBX2", bs) %>% dplyr::rename(lfc_gstcbx2 = .data[[lfc_col]], se_gstcbx2 = lfcSE, alc_gstcbx2 = aveLogCPM)

  # Inner join on bin_id — common bins across all 5 targets
  df <- k27  %>% dplyr::select(bin_id, seqnames, start, end, lfc_k27me2, se_k27me2, alc_k27me2) %>%
    dplyr::inner_join(cbx2  %>% dplyr::select(bin_id, lfc_cbx2, se_cbx2, alc_cbx2),             by = "bin_id") %>%
    dplyr::inner_join(gcbx2 %>% dplyr::select(bin_id, lfc_gstcbx2, se_gstcbx2, alc_gstcbx2),    by = "bin_id") %>%
    dplyr::inner_join(cbx7  %>% dplyr::select(bin_id, lfc_cbx7, se_cbx7, alc_cbx7),             by = "bin_id") %>%
    dplyr::inner_join(ezh2  %>% dplyr::select(bin_id, lfc_ezh2, se_ezh2, alc_ezh2),             by = "bin_id") %>%
    dplyr::filter(complete.cases(.))

  has_k27me3 <- file.exists(file.path(lfc_dir, sprintf("lfc_bins_k27me3_%dbp.rds", bs))) && !grepl("lost_k27me3", lfc_dir)
  if (has_k27me3) {
    k27m3 <- load_lfc("K27me3", bs, od = lfc_dir) %>% dplyr::rename(lfc_k27me3 = .data[[lfc_col]], se_k27me3 = lfcSE, alc_k27me3 = aveLogCPM)
    df <- df %>% dplyr::inner_join(k27m3 %>% dplyr::select(bin_id, lfc_k27me3, se_k27me3, alc_k27me3), by = "bin_id")
  }

  cat(sprintf("\nCommon bins across all targets at %d bp: N = %d (has_k27me3 = %s)\n", bs, nrow(df), has_k27me3))
  df
}

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Helper: Cluster-Robust Standard Errors and Confidence Intervals via sandwich package
# ---------------------------------------------------------------------------
calc_cluster_vcov <- function(mod, cluster) {
  cluster <- as.factor(cluster)
  G <- nlevels(cluster)
  
  # Use industry-standard sandwich::vcovCL (HC1 adjustment)
  V <- sandwich::vcovCL(mod, cluster = cluster, type = "HC1")
  se <- sqrt(pmax(0, diag(V)))
  names(se) <- names(coef(mod))
  
  t_crit <- qt(0.975, df = max(1, G - 1))
  ci_low  <- coef(mod) - t_crit * se
  ci_high <- coef(mod) + t_crit * se
  
  list(se = se, ci_low = ci_low, ci_high = ci_high, df = G - 1, vcov = V)
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
  c_factor <- N / (N - k)
  V <- c_factor * (B %*% M_hac %*% B)
  colnames(V) <- names(coef(mod))
  rownames(V) <- names(coef(mod))
  se <- sqrt(pmax(0, diag(V)))
  names(se) <- names(coef(mod))
  
  t_crit <- qt(0.975, df = max(1, N - k))
  ci_low  <- coef(mod) - t_crit * se
  ci_high <- coef(mod) + t_crit * se
  
  list(se = se, ci_low = ci_low, ci_high = ci_high, df = N - k, vcov = V)
}

# ---------------------------------------------------------------------------
# Helper: Block2Mb Standard Errors
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Helper: 5-fold CV deviance explained for GAMs
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Helper: exact VIF from linear regression model
# ---------------------------------------------------------------------------
calc_vif_from_mod <- function(mod, df_data) {
  X_mat <- model.matrix(mod)[, -1, drop = FALSE]
  if (ncol(X_mat) < 2) return(setNames(1.0, colnames(X_mat)))
  w <- weights(mod)
  if (is.null(w)) w <- rep(1, nrow(X_mat))
  cov_w <- cov.wt(X_mat, wt = w)$cov
  cor_w <- cov2cor(cov_w)
  vifs <- tryCatch(diag(solve(cor_w)), error = function(e) rep(NA_real_, ncol(X_mat)))
  names(vifs) <- colnames(X_mat)
  vifs
}

# ---------------------------------------------------------------------------
# Helper: annotate plot with correlation stats
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Helper: annotate plot with correlation stats and exact beta slope
# ---------------------------------------------------------------------------
corr_label <- function(x, y, N, w = NULL) {
  r_spear <- cor(x, y, method = "spearman")
  if (!is.null(w)) {
    b_val <- coef(lm(y ~ x, weights = w))[2]
    sprintf("Spearman \u03c1 = %.3f\nIVW \u03b2 = %.3f\nN = %s",
            r_spear, b_val, formatC(N, format = "d", big.mark = ","))
  } else {
    b_val <- coef(lm(y ~ x))[2]
    sprintf("Spearman \u03c1 = %.3f\nOLS \u03b2 = %.3f\nN = %s",
            r_spear, b_val, formatC(N, format = "d", big.mark = ","))
  }
}

# ---------------------------------------------------------------------------
# Helper: hexbin scatter with correlation stats and beta slope
# ---------------------------------------------------------------------------
make_scatter <- function(df, xcol, ycol, xlab, ylab, title_str, wcol = NULL) {
  x <- df[[xcol]]; y <- df[[ycol]]
  w <- if (!is.null(wcol)) df[[wcol]] else NULL
  lbl <- corr_label(x, y, nrow(df), w)

  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    stat_bin_hex(bins = 80, aes(fill = after_stat(log10(count)))) +
    scale_fill_gradientn(
      colours = c("#f0f0f0", "#4575b4", "#d73027"),
      name    = expression(log[10](N))
    )
  if (!is.null(wcol)) {
    p <- p + geom_smooth(aes(weight = .data[[wcol]]), method = "lm", se = FALSE, color = "black", linewidth = 0.9, linetype = "dashed")
  } else {
    p <- p + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8, linetype = "dashed")
  }
  p + annotate("text", x = Inf, y = -Inf, label = lbl,
               hjust = 1.05, vjust = -0.3, size = 3.8, family = "Helvetica") +
    labs(x = xlab, y = ylab, title = title_str) +
    my_theme
}

# ---------------------------------------------------------------------------
# Helper: save to SVG via base R device
# ---------------------------------------------------------------------------
save_svg <- function(p, path, width = 5.5, height = 5) {
  svg(path, width = width, height = height)
  print(p)
  invisible(dev.off())
  cat("Saved:", path, "\n")
}

# ---------------------------------------------------------------------------
# Helper: print correlation + p-value table row
# ---------------------------------------------------------------------------
print_corr_row <- function(label, x, y) {
  n  <- length(x)
  rs <- cor.test(x, y, method = "spearman", exact = FALSE)
  rp <- cor.test(x, y, method = "pearson")
  cat(sprintf("  %-36s  Spearman ρ=%+.4f (p=%.2e)  Pearson r=%+.4f (p=%.2e)  N=%d\n",
              label, rs$estimate, rs$p.value, rp$estimate, rp$p.value, n))
}

# ---------------------------------------------------------------------------
# MAIN ANALYSIS FUNCTION (parameterised by bin size)
# ---------------------------------------------------------------------------
run_analysis <- function(bs, label = "", save_core_plots = TRUE) {

  cat(sprintf("\n%s\n### ANALYSIS — %d-bp bins %s###\n%s\n",
              strrep("=", 60), bs, if (nzchar(label)) paste0("[", label, "] ") else "", strrep("=", 60)))

  df <- merge_targets(bs)
  N  <- nrow(df)

  # Precision proxy from aveLogCPM (read count abundance ~ 2^aveLogCPM) instead of derived lfcSE
  df$w_cbx2    <- 2^df$alc_cbx2
  df$w_gstcbx2 <- 2^df$alc_gstcbx2
  df$w_cbx7    <- 2^df$alc_cbx7

  # Combined IVW precision across response and predictors (1 / (sigma_y^2 + sigma_x1^2 + sigma_x2^2))
  df$w_comb_cbx2    <- 1 / (2^(-df$alc_cbx2)    + 2^(-df$alc_k27me2) + 2^(-df$alc_ezh2))
  df$w_comb_gstcbx2 <- 1 / (2^(-df$alc_gstcbx2) + 2^(-df$alc_k27me2) + 2^(-df$alc_ezh2))
  df$w_comb_cbx7    <- 1 / (2^(-df$alc_cbx7)    + 2^(-df$alc_k27me2) + 2^(-df$alc_ezh2))

  # ---- Core plots ----------------------------------------------------------
  if (save_core_plots) {
    pA <- make_scatter(df,
                       xcol = "lfc_k27me2", ycol = "lfc_cbx2",
                       xlab = expression(Delta ~ "H3K27me2  " * (log[2]*FC)),
                       ylab = expression(Delta ~ "WT CBX2  " * (log[2]*FC)),
                       title_str = "ΔCBX2 (WT) vs ΔH3K27me2 [OLS]")

    pB <- make_scatter(df,
                       xcol = "lfc_k27me2", ycol = "lfc_cbx7",
                       xlab = expression(Delta ~ "H3K27me2  " * (log[2]*FC)),
                       ylab = expression(Delta ~ "GST-CBX7  " * (log[2]*FC)),
                       title_str = "ΔCBX7 vs ΔH3K27me2 [OLS]")

    pC <- make_scatter(df,
                       xcol = "lfc_k27me2", ycol = "lfc_gstcbx2",
                       xlab = expression(Delta ~ "H3K27me2  " * (log[2]*FC)),
                       ylab = expression(Delta ~ "GST-CBX2  " * (log[2]*FC)),
                       title_str = "ΔCBX2 (GST) vs ΔH3K27me2 [OLS]")

    # Also make target vs EZH2 OLS scatter plots with annotated beta values
    pA_e <- make_scatter(df, "lfc_ezh2", "lfc_cbx2",
                         expression(Delta ~ "EZH2  " * (log[2]*FC)),
                         expression(Delta ~ "WT CBX2  " * (log[2]*FC)),
                         "ΔCBX2 (WT) vs ΔEZH2 [OLS]")
    pB_e <- make_scatter(df, "lfc_ezh2", "lfc_cbx7",
                         expression(Delta ~ "EZH2  " * (log[2]*FC)),
                         expression(Delta ~ "GST-CBX7  " * (log[2]*FC)),
                         "ΔCBX7 vs ΔEZH2 [OLS]")
    pC_e <- make_scatter(df, "lfc_ezh2", "lfc_gstcbx2",
                         expression(Delta ~ "EZH2  " * (log[2]*FC)),
                         expression(Delta ~ "GST-CBX2  " * (log[2]*FC)),
                         "ΔCBX2 (GST) vs ΔEZH2 [OLS]")

    sfx <- if (nzchar(label)) paste0("_", tolower(gsub(" ", "_", label))) else ""
    save_svg(pA, file.path(out_dir, sprintf("delta_vs_delta_cbx2_vs_k27me2_%dbp%s.svg", bs, sfx)))
    save_svg(pB, file.path(out_dir, sprintf("delta_vs_delta_cbx7_vs_k27me2_%dbp%s.svg", bs, sfx)))
    save_svg(pC, file.path(out_dir, sprintf("delta_vs_delta_gstcbx2_vs_k27me2_%dbp%s.svg", bs, sfx)))
    save_svg(pA_e, file.path(out_dir, sprintf("delta_vs_delta_cbx2_vs_ezh2_%dbp%s.svg", bs, sfx)))
    save_svg(pB_e, file.path(out_dir, sprintf("delta_vs_delta_cbx7_vs_ezh2_%dbp%s.svg", bs, sfx)))
    save_svg(pC_e, file.path(out_dir, sprintf("delta_vs_delta_gstcbx2_vs_ezh2_%dbp%s.svg", bs, sfx)))

    pA_w <- make_scatter(df, "lfc_k27me2", "lfc_cbx2",
                         expression(Delta ~ "H3K27me2  " * (log[2]*FC)),
                         expression(Delta ~ "WT CBX2  " * (log[2]*FC)),
                         "ΔCBX2 (WT) vs ΔH3K27me2 [Combined IVW]", wcol = "w_comb_cbx2")
    pB_w <- make_scatter(df, "lfc_k27me2", "lfc_cbx7",
                         expression(Delta ~ "H3K27me2  " * (log[2]*FC)),
                         expression(Delta ~ "GST-CBX7  " * (log[2]*FC)),
                         "ΔCBX7 vs ΔH3K27me2 [Combined IVW]", wcol = "w_comb_cbx7")
    pC_w <- make_scatter(df, "lfc_k27me2", "lfc_gstcbx2",
                         expression(Delta ~ "H3K27me2  " * (log[2]*FC)),
                         expression(Delta ~ "GST-CBX2  " * (log[2]*FC)),
                         "ΔCBX2 (GST) vs ΔH3K27me2 [Combined IVW]", wcol = "w_comb_gstcbx2")

    pA_ew <- make_scatter(df, "lfc_ezh2", "lfc_cbx2",
                          expression(Delta ~ "EZH2  " * (log[2]*FC)),
                          expression(Delta ~ "WT CBX2  " * (log[2]*FC)),
                          "ΔCBX2 (WT) vs ΔEZH2 [Combined IVW]", wcol = "w_comb_cbx2")
    pB_ew <- make_scatter(df, "lfc_ezh2", "lfc_cbx7",
                          expression(Delta ~ "EZH2  " * (log[2]*FC)),
                          expression(Delta ~ "GST-CBX7  " * (log[2]*FC)),
                          "ΔCBX7 vs ΔEZH2 [Combined IVW]", wcol = "w_comb_cbx7")
    pC_ew <- make_scatter(df, "lfc_ezh2", "lfc_gstcbx2",
                          expression(Delta ~ "EZH2  " * (log[2]*FC)),
                          expression(Delta ~ "GST-CBX2  " * (log[2]*FC)),
                          "ΔCBX2 (GST) vs ΔEZH2 [Combined IVW]", wcol = "w_comb_gstcbx2")

    save_svg(pA_w, file.path(out_dir, sprintf("delta_vs_delta_cbx2_vs_k27me2_combined_ivw_%dbp%s.svg", bs, sfx)))
    save_svg(pB_w, file.path(out_dir, sprintf("delta_vs_delta_cbx7_vs_k27me2_combined_ivw_%dbp%s.svg", bs, sfx)))
    save_svg(pC_w, file.path(out_dir, sprintf("delta_vs_delta_gstcbx2_vs_k27me2_combined_ivw_%dbp%s.svg", bs, sfx)))
    save_svg(pA_ew, file.path(out_dir, sprintf("delta_vs_delta_cbx2_vs_ezh2_combined_ivw_%dbp%s.svg", bs, sfx)))
    save_svg(pB_ew, file.path(out_dir, sprintf("delta_vs_delta_cbx7_vs_ezh2_combined_ivw_%dbp%s.svg", bs, sfx)))
    save_svg(pC_ew, file.path(out_dir, sprintf("delta_vs_delta_gstcbx2_vs_ezh2_combined_ivw_%dbp%s.svg", bs, sfx)))
  }

  # ---- Raw correlation summary --------------------------------------------
  cat("\n--- Marginal Correlations ---\n")
  print_corr_row("WT CBX2 ~ K27me2",  df$lfc_k27me2, df$lfc_cbx2)
  print_corr_row("GST-CBX2 ~ K27me2", df$lfc_k27me2, df$lfc_gstcbx2)
  print_corr_row("GST-CBX7 ~ K27me2", df$lfc_k27me2, df$lfc_cbx7)
  print_corr_row("WT CBX2 ~ EZH2",    df$lfc_ezh2,   df$lfc_cbx2)
  print_corr_row("GST-CBX2 ~ EZH2",   df$lfc_ezh2,   df$lfc_gstcbx2)
  print_corr_row("GST-CBX7 ~ EZH2",   df$lfc_ezh2,   df$lfc_cbx7)
  print_corr_row("K27me2 ~ EZH2",     df$lfc_ezh2,   df$lfc_k27me2)

  # ---- Extra 1: Partial correlation (ppcor) --------------------------------
  cat("\n--- Extra 1: Partial Correlation (controlling for ΔEZH2) ---\n")
  pc_cbx2    <- pcor.test(df$lfc_cbx2,    df$lfc_k27me2, df$lfc_ezh2, method = "spearman")
  pc_gstcbx2 <- pcor.test(df$lfc_gstcbx2, df$lfc_k27me2, df$lfc_ezh2, method = "spearman")
  pc_cbx7    <- pcor.test(df$lfc_cbx7,    df$lfc_k27me2, df$lfc_ezh2, method = "spearman")
  cat(sprintf("  ΔCBX2 (WT)  ~ ΔK27me2 | ΔEZH2:  ρ_partial=%+.4f  p=%.2e  df=%d\n",
              pc_cbx2$estimate, pc_cbx2$p.value, pc_cbx2$n - pc_cbx2$gp - 2))
  cat(sprintf("  ΔCBX2 (GST) ~ ΔK27me2 | ΔEZH2:  ρ_partial=%+.4f  p=%.2e  df=%d\n",
              pc_gstcbx2$estimate, pc_gstcbx2$p.value, pc_gstcbx2$n - pc_gstcbx2$gp - 2))
  cat(sprintf("  ΔCBX7 (GST) ~ ΔK27me2 | ΔEZH2:  ρ_partial=%+.4f  p=%.2e  df=%d\n",
              pc_cbx7$estimate, pc_cbx7$p.value, pc_cbx7$n - pc_cbx7$gp - 2))

  has_k27me3 <- "lfc_k27me3" %in% names(df)
  if (has_k27me3) {
    df$w_comb_cbx2    <- 1 / (2^(-df$alc_cbx2)    + 2^(-df$alc_k27me2) + 2^(-df$alc_ezh2) + 2^(-df$alc_k27me3))
    df$w_comb_gstcbx2 <- 1 / (2^(-df$alc_gstcbx2) + 2^(-df$alc_k27me2) + 2^(-df$alc_ezh2) + 2^(-df$alc_k27me3))
    df$w_comb_cbx7    <- 1 / (2^(-df$alc_cbx7)    + 2^(-df$alc_k27me2) + 2^(-df$alc_ezh2) + 2^(-df$alc_k27me3))
  }
  
  pred_cols   <- if (has_k27me3) c("lfc_k27me2", "lfc_ezh2", "lfc_k27me3") else c("lfc_k27me2", "lfc_ezh2")
  formula_rhs <- paste(pred_cols, collapse = " + ")

  cat("\n--- Extra 1b: Weighted Linear Regression Sensitivity Analysis ---\n")
  lm_ols_cbx2    <- lm(as.formula(paste("lfc_cbx2 ~", formula_rhs)), data = df)
  lm_ols_gstcbx2 <- lm(as.formula(paste("lfc_gstcbx2 ~", formula_rhs)), data = df)
  lm_ols_cbx7    <- lm(as.formula(paste("lfc_cbx7 ~", formula_rhs)), data = df)

  lm_resp_cbx2    <- lm(as.formula(paste("lfc_cbx2 ~", formula_rhs)), data = df, weights = w_cbx2)
  lm_resp_gstcbx2 <- lm(as.formula(paste("lfc_gstcbx2 ~", formula_rhs)), data = df, weights = w_gstcbx2)
  lm_resp_cbx7    <- lm(as.formula(paste("lfc_cbx7 ~", formula_rhs)), data = df, weights = w_cbx7)

  lm_comb_cbx2    <- lm(as.formula(paste("lfc_cbx2 ~", formula_rhs)), data = df, weights = w_comb_cbx2)
  lm_comb_gstcbx2 <- lm(as.formula(paste("lfc_gstcbx2 ~", formula_rhs)), data = df, weights = w_comb_gstcbx2)
  lm_comb_cbx7    <- lm(as.formula(paste("lfc_cbx7 ~", formula_rhs)), data = df, weights = w_comb_cbx7)

  # Spatial block ID: 2-Mb blocks along each chromosome
  df$block_2mb <- paste0(df$seqnames, "_", floor(df$start / 2e6))

  models_list <- list(
    "WT CBX2"  = list("Unweighted OLS" = lm_ols_cbx2, "Response-only IVW" = lm_resp_cbx2, "Combined IVW" = lm_comb_cbx2),
    "GST-CBX2" = list("Unweighted OLS" = lm_ols_gstcbx2, "Response-only IVW" = lm_resp_gstcbx2, "Combined IVW" = lm_comb_gstcbx2),
    "GST-CBX7" = list("Unweighted OLS" = lm_ols_cbx7, "Response-only IVW" = lm_resp_cbx7, "Combined IVW" = lm_comb_cbx7)
  )

  spatial_ci_table <- data.frame()

  for (tname in names(models_list)) {
    cat(sprintf("\n=== Target: %s ===\n", tname))
    for (mname in names(models_list[[tname]])) {
      mod <- models_list[[tname]][[mname]]
      s <- summary(mod)
      cf <- s$coefficients
      vif_vals <- calc_vif_from_mod(mod, df)

      for (pred_name in pred_cols) {
        b_p <- cf[pred_name, "Estimate"]; se_naive_p <- cf[pred_name, "Std. Error"]
        ci_naive_p_low <- b_p - 1.96 * se_naive_p; ci_naive_p_high <- b_p + 1.96 * se_naive_p

        res_chr <- calc_cluster_vcov(mod, df$seqnames)
        se_chr_p <- res_chr$se[pred_name]; ci_chr_p_low <- res_chr$ci_low[pred_name]; ci_chr_p_high <- res_chr$ci_high[pred_name]

        res_blk <- calc_cluster_vcov(mod, df$block_2mb)
        se_blk_p <- res_blk$se[pred_name]; ci_blk_p_low <- res_blk$ci_low[pred_name]; ci_blk_p_high <- res_blk$ci_high[pred_name]

        res_hac <- calc_spatial_hac_vcov(mod, df$seqnames, df$start, max_lag_dist = 500000)
        se_hac_p <- res_hac$se[pred_name]; ci_hac_p_low <- res_hac$ci_low[pred_name]; ci_hac_p_high <- res_hac$ci_high[pred_name]

        t_hac_p <- b_p / pmax(se_hac_p, 1e-6); p_hac_p <- 2 * (1 - pnorm(abs(t_hac_p)))

        spatial_ci_table <- rbind(spatial_ci_table, data.frame(
          Target = tname, Weighting = mname, Predictor = pred_name, Estimate = b_p, VIF = vif_vals[pred_name],
          SE_Naive = se_naive_p, CI_Naive_Low = ci_naive_p_low, CI_Naive_High = ci_naive_p_high,
          SE_Chr = se_chr_p, CI_Chr_Low = ci_chr_p_low, CI_Chr_High = ci_chr_p_high,
          SE_Block2Mb = se_blk_p, CI_Block2Mb_Low = ci_blk_p_low, CI_Block2Mb_High = ci_blk_p_high,
          SE_SpatialHAC = se_hac_p, CI_SpatialHAC_Low = ci_hac_p_low, CI_SpatialHAC_High = ci_hac_p_high,
          t_SpatialHAC = t_hac_p, p_SpatialHAC = p_hac_p,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  sfx <- if (nzchar(label)) paste0("_", tolower(gsub(" ", "_", label))) else ""
  tsv_spatial_path <- file.path(out_dir, sprintf("regression_summary_with_spatial_CIs_%dbp%s.tsv", bs, sfx))
  write.table(spatial_ci_table, tsv_spatial_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nSaved spatial CI regression table:", tsv_spatial_path, "\n")

  # Coefficient sensitivity plots (plotting estimates with Spatial HAC and 2Mb Block CI bars)
  coef_k27 <- spatial_ci_table[spatial_ci_table$Predictor == "lfc_k27me2", ]
  coef_k27$Target <- factor(coef_k27$Target, levels = c("WT CBX2", "GST-CBX2", "GST-CBX7"))
  coef_k27$Weighting <- factor(coef_k27$Weighting, levels = c("Unweighted OLS", "Response-only IVW", "Combined IVW"))

  p_coef_k27 <- ggplot(coef_k27, aes(x = Target, y = Estimate, color = Weighting, shape = Weighting)) +
    geom_pointrange(aes(ymin = CI_Block2Mb_Low, ymax = CI_Block2Mb_High), position = position_dodge(width = 0.5), size = 0.8) +
    scale_color_manual(values = c("Unweighted OLS" = "#7fcdbb", "Response-only IVW" = "#2c7fb8", "Combined IVW" = "#d95f02")) +
    labs(x = "Target Protein",
         y = expression(beta["ΔH3K27me2"] ~ "Coefficient (controlling for ΔEZH2)"),
         title = "Sensitivity of ΔH3K27me2 Coefficient Across Weighting Schemes",
         subtitle = "Error bars = 95% Confidence Intervals accounting for 2-Mb Spatial Autocorrelation") +
    my_theme +
    theme(legend.position = "top")

  save_svg(p_coef_k27, file.path(out_dir, sprintf("delta_vs_delta_k27me2_coef_sensitivity_%dbp%s.svg", bs, sfx)), width = 8.5, height = 5.5)

  coef_ezh2 <- spatial_ci_table[spatial_ci_table$Predictor == "lfc_ezh2", ]
  coef_ezh2$Target <- factor(coef_ezh2$Target, levels = c("WT CBX2", "GST-CBX2", "GST-CBX7"))
  coef_ezh2$Weighting <- factor(coef_ezh2$Weighting, levels = c("Unweighted OLS", "Response-only IVW", "Combined IVW"))

  p_coef_ezh2 <- ggplot(coef_ezh2, aes(x = Target, y = Estimate, color = Weighting, shape = Weighting)) +
    geom_pointrange(aes(ymin = CI_Block2Mb_Low, ymax = CI_Block2Mb_High), position = position_dodge(width = 0.5), size = 0.8) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    scale_color_manual(values = c("Unweighted OLS" = "#7fcdbb", "Response-only IVW" = "#2c7fb8", "Combined IVW" = "#d95f02")) +
    labs(x = "Target Protein",
         y = expression(beta["ΔEZH2"] ~ "Coefficient (controlling for ΔK27me2)"),
         title = "Sensitivity of ΔEZH2 Coefficient Across Weighting Schemes",
         subtitle = "Error bars = 95% Confidence Intervals accounting for 2-Mb Spatial Autocorrelation") +
    my_theme +
    theme(legend.position = "top")

  save_svg(p_coef_ezh2, file.path(out_dir, sprintf("delta_vs_delta_ezh2_coef_sensitivity_%dbp%s.svg", bs, sfx)), width = 8.5, height = 5.5)

  # Multivariable Regression Forest Plot comparing K27me2 vs EZH2 side by side
  spatial_ci_table$Predictor_Label <- ifelse(spatial_ci_table$Predictor == "lfc_k27me2", "ΔH3K27me2 (β_K27me2)", "ΔEZH2 (β_EZH2)")
  spatial_ci_table$Target <- factor(spatial_ci_table$Target, levels = c("WT CBX2", "GST-CBX2", "GST-CBX7"))
  spatial_ci_table$Weighting <- factor(spatial_ci_table$Weighting, levels = c("Unweighted OLS", "Response-only IVW", "Combined IVW"))

  p_forest <- ggplot(spatial_ci_table, aes(x = Predictor_Label, y = Estimate, color = Weighting, shape = Weighting)) +
    geom_pointrange(aes(ymin = CI_SpatialHAC_Low, ymax = CI_SpatialHAC_High), position = position_dodge(width = 0.6), size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.6) +
    facet_wrap(~ Target, nrow = 1) +
    scale_color_manual(values = c("Unweighted OLS" = "#7fcdbb", "Response-only IVW" = "#2c7fb8", "Combined IVW" = "#d95f02")) +
    labs(x = "Predictor in Multivariable Regression",
         y = "Partial Regression Coefficient (β ± 95% Spatial HAC CI)",
         title = "Multivariable Regression Forest Plot: Independent Contributions to ΔCBX",
         subtitle = "Error bars = 95% CIs accounting for 1D Spatial HAC Autocorrelation") +
    my_theme +
    theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))

  save_svg(p_forest, file.path(out_dir, sprintf("coef_forest_plot_%dbp%s.svg", bs, sfx)), width = 10, height = 5.5)

  # ---- Extra 1c: Colocalization Independence & Variance Decomposition (sr²) ----
  cat("\n--- Extra 1c: Colocalization Independence & Variance Decomposition (sr²) ---\n")
  wtd_var <- function(x, w = NULL) {
    if (is.null(w)) return(var(x, na.rm = TRUE))
    w <- w / sum(w)
    m <- sum(w * x)
    sum(w * (x - m)^2) / (1 - sum(w^2))
  }

  indep_metrics <- data.frame()
  decomp_long   <- data.frame()

  for (tname in names(models_list)) {
    yvar <- switch(tname, "WT CBX2" = "lfc_cbx2", "GST-CBX2" = "lfc_gstcbx2", "GST-CBX7" = "lfc_cbx7")
    wvar <- switch(tname, "WT CBX2" = "w_comb_cbx2", "GST-CBX2" = "w_comb_gstcbx2", "GST-CBX7" = "w_comb_cbx7")

    for (mname in names(models_list[[tname]])) {
      mod_full <- models_list[[tname]][[mname]]
      w_vec    <- if (mname == "Unweighted OLS") NULL else if (mname == "Combined IVW") df[[wvar]] else switch(tname, "WT CBX2"=df$w_cbx2, "GST-CBX2"=df$w_gstcbx2, "GST-CBX7"=df$w_cbx7)
      r2_full  <- summary(mod_full)$r.squared

      for (pred_name in pred_cols) {
        r2_no_pred <- if (!is.null(w_vec)) summary(lm(as.formula(paste(yvar, "~", paste(setdiff(pred_cols, pred_name), collapse = " + "))), data = df, weights = w_vec))$r.squared else summary(lm(as.formula(paste(yvar, "~", paste(setdiff(pred_cols, pred_name), collapse = " + "))), data = df))$r.squared
        sr2_unique <- pmax(0, r2_full - r2_no_pred)
        
        sd_y <- sqrt(wtd_var(df[[yvar]], w_vec))
        sd_p <- sqrt(wtd_var(df[[pred_name]], w_vec))
        b_p  <- coef(mod_full)[pred_name]
        std_b_p <- b_p * (sd_p / sd_y)
        
        row_p <- spatial_ci_table[spatial_ci_table$Target == tname & spatial_ci_table$Weighting == mname & spatial_ci_table$Predictor == pred_name, ]
        
        indep_metrics <- rbind(indep_metrics, data.frame(
          Target = tname, Weighting = mname, Predictor = pred_name,
          Beta = b_p, SE_SpatialHAC = row_p$SE_SpatialHAC, p_SpatialHAC = row_p$p_SpatialHAC, Std_Beta = std_b_p,
          R2_Full = r2_full, Unique_sr2 = sr2_unique,
          stringsAsFactors = FALSE
        ))
        decomp_long <- rbind(decomp_long, data.frame(
          Target = tname, Weighting = mname, Component = sprintf("Unique %s (sr²)", pred_name), Variance_Explained = sr2_unique, stringsAsFactors = FALSE
        ))
      }
      sum_unique <- sum(decomp_long$Variance_Explained[decomp_long$Target == tname & decomp_long$Weighting == mname])
      sr2_shared <- pmax(0, r2_full - sum_unique)
      decomp_long <- rbind(decomp_long, data.frame(
        Target = tname, Weighting = mname, Component = "Shared / Overlap", Variance_Explained = sr2_shared, stringsAsFactors = FALSE
      ))
    }
  }

  tsv_indep_path <- file.path(out_dir, sprintf("colocalization_independence_metrics_%dbp%s.tsv", bs, sfx))
  write.table(indep_metrics, tsv_indep_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Saved colocalization independence metrics table:", tsv_indep_path, "\n")

  decomp_long$Component <- factor(decomp_long$Component, levels = c("Shared / Overlap", "Unique ΔEZH2 (sr²)", "Unique ΔH3K27me2 (sr²)"))
  decomp_long$Target    <- factor(decomp_long$Target, levels = c("WT CBX2", "GST-CBX2", "GST-CBX7"))
  decomp_long$Weighting <- factor(decomp_long$Weighting, levels = c("Unweighted OLS", "Response-only IVW", "Combined IVW"))

  p_decomp <- ggplot(decomp_long, aes(x = Target, y = Variance_Explained, fill = Component)) +
    geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black") +
    facet_wrap(~ Weighting, nrow = 1) +
    scale_fill_manual(values = c("Unique ΔH3K27me2 (sr²)" = "#1b9e77", "Unique ΔEZH2 (sr²)" = "#d95f02", "Shared / Overlap" = "#7570b3")) +
    labs(x = "Target Protein", y = "Proportion of Total Variance (R²) Explained",
         title = "Variance Decomposition (sr²): Independent vs Shared Contributions to ΔCBX",
         subtitle = "Unique sr² = Variance explained exclusively by predictor holding the other constant") +
    my_theme + theme(legend.position = "top")

  save_svg(p_decomp, file.path(out_dir, sprintf("variance_decomposition_sr2_%dbp%s.svg", bs, sfx)), width = 10, height = 5.5)

  # ---- Extra 2: Steiger's Z test (cocor) -----------------------------------
  cat("\n--- Extra 2: Steiger's Z Tests (Comparing r(target, K27me2)) ---\n")
  r_wt_k27  <- cor(df$lfc_cbx2,    df$lfc_k27me2, method = "pearson")
  r_gst_k27 <- cor(df$lfc_gstcbx2, df$lfc_k27me2, method = "pearson")
  r_cbx7_k27<- cor(df$lfc_cbx7,    df$lfc_k27me2, method = "pearson")

  r_wt_gst  <- cor(df$lfc_cbx2,    df$lfc_gstcbx2, method = "pearson")
  r_wt_cbx7 <- cor(df$lfc_cbx2,    df$lfc_cbx7,    method = "pearson")
  r_gst_cbx7<- cor(df$lfc_gstcbx2, df$lfc_cbx7,    method = "pearson")

  cat(sprintf("  Pearson r(WT CBX2, K27me2)  = %.4f\n  Pearson r(GST-CBX2, K27me2) = %.4f\n  Pearson r(GST-CBX7, K27me2) = %.4f\n\n",
              r_wt_k27, r_gst_k27, r_cbx7_k27))

  cat("  -> Comparison 1: WT CBX2 vs GST-CBX2 (Does recombinant/N-term GST-CBX2 correlate better than WT CBX2?)\n")
  cc_wt_gst <- cocor.dep.groups.overlap(r.jk = r_wt_k27, r.jh = r_gst_k27, r.kh = r_wt_gst, n = N)
  print(cc_wt_gst)

  cat("\n  -> Comparison 2: GST-CBX2 vs GST-CBX7 (Comparing the two GST-tagged constructs over shared K27me2)\n")
  cc_gst2_7 <- cocor.dep.groups.overlap(r.jk = r_gst_k27, r.jh = r_cbx7_k27, r.kh = r_gst_cbx7, n = N)
  print(cc_gst2_7)

  cat("\n  -> Comparison 3: WT CBX2 vs GST-CBX7 (Baseline comparison)\n")
  cc_wt_7 <- cocor.dep.groups.overlap(r.jk = r_wt_k27, r.jh = r_cbx7_k27, r.kh = r_wt_cbx7, n = N)
  print(cc_wt_7)

  # ---- Extra 3: EZH2 & K27me3 tertile stratifications ---------------------
  cat("\n--- Extra 3: Stratification by ΔEZH2 Tertile ---\n")
  df$ezh2_tertile <- cut(df$lfc_ezh2,
                          breaks = quantile(df$lfc_ezh2, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                          labels = c("Stable EZH2", "Moderate EZH2 loss", "High EZH2 loss"),
                          include.lowest = TRUE)

  strat_summary <- data.frame()
  for (ter in levels(df$ezh2_tertile)) {
    sub_df <- df[df$ezh2_tertile == ter, ]
    n_sub  <- nrow(sub_df)
    
    for (t_pair in list(c("lfc_cbx2", "WT CBX2"), c("lfc_gstcbx2", "GST-CBX2"), c("lfc_cbx7", "GST-CBX7"))) {
      t_col  <- t_pair[1]
      t_name <- t_pair[2]
      
      m_fit <- lm(as.formula(paste(t_col, "~ lfc_k27me2")), data = sub_df)
      hac_sub <- calc_spatial_hac_vcov(m_fit, seqnames = sub_df$seqnames, pos = sub_df$start, max_lag_dist = 500000)
      blk_sub <- calc_block2mb_se(m_fit, sub_df)
      
      b_val   <- coef(m_fit)["lfc_k27me2"]
      sd_y_sub <- sd(sub_df[[t_col]], na.rm = TRUE)
      sd_x_sub <- sd(sub_df$lfc_k27me2, na.rm = TRUE)
      std_b_sub <- b_val * (sd_x_sub / sd_y_sub)
      
      g_fit <- tryCatch(gam(as.formula(paste(t_col, "~ s(lfc_k27me2, bs='cs')")), data = sub_df), error = function(e) NULL)
      gam_dev <- if (!is.null(g_fit)) summary(g_fit)$dev.expl else NA_real_
      cv_dev  <- gam_cv_dev_expl(paste(t_col, "~ s(lfc_k27me2, bs='cs')"), sub_df, k = 5)
      
      strat_summary <- rbind(strat_summary, data.frame(
        ezh2_tertile = ter, Target = t_name, N = n_sub,
        Spearman_rho = cor(sub_df[[t_col]], sub_df$lfc_k27me2, method = "spearman"),
        Beta = b_val, Std_Beta = std_b_sub,
        SE_SpatialHAC = hac_sub$se["lfc_k27me2"], CI_SpatialHAC_low = hac_sub$ci_low["lfc_k27me2"], CI_SpatialHAC_high = hac_sub$ci_high["lfc_k27me2"],
        SE_Block2Mb   = blk_sub$se["lfc_k27me2"], CI_Block2Mb_low   = blk_sub$ci_low["lfc_k27me2"], CI_Block2Mb_high   = blk_sub$ci_high["lfc_k27me2"],
        OLS_R2 = summary(m_fit)$r.squared, GAM_Dev_Expl = gam_dev, GAM_CV_Dev_Expl = cv_dev,
        stringsAsFactors = FALSE
      ))
    }
  }
  print(strat_summary)
  tsv_strat_path <- file.path(out_dir, sprintf("delta_vs_delta_ezh2_tertile_stratification_%dbp%s.tsv", bs, sfx))
  write.table(strat_summary, tsv_strat_path, sep = "\t", quote = FALSE, row.names = FALSE)

  if (has_k27me3) {
    cat("\n--- H3K27me3 Tertile Stratification ---\n")
    df$k27me3_tertile <- cut(df$lfc_k27me3,
                              breaks = quantile(df$lfc_k27me3, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                              labels = c("Low ΔH3K27me3", "Moderate ΔH3K27me3", "High ΔH3K27me3"),
                              include.lowest = TRUE)
    strat_m3_summary <- data.frame()
    for (ter in levels(df$k27me3_tertile)) {
      sub_df <- df[df$k27me3_tertile == ter, ]
      n_sub  <- nrow(sub_df)
      for (t_pair in list(c("lfc_cbx2", "WT CBX2"), c("lfc_gstcbx2", "GST-CBX2"), c("lfc_cbx7", "GST-CBX7"))) {
        t_col  <- t_pair[1]
        t_name <- t_pair[2]
        m_fit <- lm(as.formula(paste(t_col, "~ lfc_k27me2")), data = sub_df)
        hac_sub <- calc_spatial_hac_vcov(m_fit, seqnames = sub_df$seqnames, pos = sub_df$start, max_lag_dist = 500000)
        blk_sub <- calc_block2mb_se(m_fit, sub_df)
        b_val   <- coef(m_fit)["lfc_k27me2"]
        sd_y_sub <- sd(sub_df[[t_col]], na.rm = TRUE)
        sd_x_sub <- sd(sub_df$lfc_k27me2, na.rm = TRUE)
        std_b_sub <- b_val * (sd_x_sub / sd_y_sub)
        g_fit <- tryCatch(gam(as.formula(paste(t_col, "~ s(lfc_k27me2, bs='cs')")), data = sub_df), error = function(e) NULL)
        gam_dev <- if (!is.null(g_fit)) summary(g_fit)$dev.expl else NA_real_
        cv_dev  <- gam_cv_dev_expl(paste(t_col, "~ s(lfc_k27me2, bs='cs')"), sub_df, k = 5)
        strat_m3_summary <- rbind(strat_m3_summary, data.frame(
          k27me3_tertile = ter, Target = t_name, N = n_sub,
          Spearman_rho = cor(sub_df[[t_col]], sub_df$lfc_k27me2, method = "spearman"),
          Beta = b_val, Std_Beta = std_b_sub,
          SE_SpatialHAC = hac_sub$se["lfc_k27me2"], CI_SpatialHAC_low = hac_sub$ci_low["lfc_k27me2"], CI_SpatialHAC_high = hac_sub$ci_high["lfc_k27me2"],
          SE_Block2Mb   = blk_sub$se["lfc_k27me2"], CI_Block2Mb_low   = blk_sub$ci_low["lfc_k27me2"], CI_Block2Mb_high   = blk_sub$ci_high["lfc_k27me2"],
          OLS_R2 = summary(m_fit)$r.squared, GAM_Dev_Expl = gam_dev, GAM_CV_Dev_Expl = cv_dev,
          stringsAsFactors = FALSE
        ))
      }
    }
    print(strat_m3_summary)
    tsv_strat_m3_path <- file.path(out_dir, sprintf("delta_vs_delta_k27me3_tertile_stratification_%dbp%s.tsv", bs, sfx))
    write.table(strat_m3_summary, tsv_strat_m3_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  invisible(list(df = df, pc_cbx2 = pc_cbx2, pc_gstcbx2 = pc_gstcbx2, pc_cbx7 = pc_cbx7))
}

# ---------------------------------------------------------------------------
# Run at primary bin size
# ---------------------------------------------------------------------------
res_primary <- run_analysis(bin_size, label = "", save_core_plots = TRUE)

# ---------------------------------------------------------------------------
# Extra 4: Sensitivity check at 2kb bins
# ---------------------------------------------------------------------------
cat(sprintf("\n%s\n### Extra 4: Sensitivity Check at %d-bp bins ###\n%s\n",
            strrep("=", 60), bin_size_sens, strrep("=", 60)))

# Check whether 2kb LFC files exist; skip gracefully if step_7i hasn't been run at 2kb
tryCatch({
  run_analysis(bin_size_sens, label = "sensitivity", save_core_plots = TRUE)
}, error = function(e) {
  cat(sprintf("  NOTE: 2kb sensitivity files not found (%s). Run step_7i at --bin_size 2000 first.\n",
              conditionMessage(e)))
})

cat("\n=== correlate_deltas.R complete ===\n")
