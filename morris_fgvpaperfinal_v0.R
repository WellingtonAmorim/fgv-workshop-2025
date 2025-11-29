# =========================================================
# Global Sensitivity — Morris (Elementary Effects)
# SFC model (t=600) + numerical safeguards
# Saves under Desktop/.../Graficos FGV final/baseline
# =========================================================

suppressPackageStartupMessages({
  library(sensitivity)  # Morris
  library(ggplot2)      # Plots
  library(ggrepel)      # Non-overlapping labels
  library(dplyr)
  library(tibble)       # Printing tibble with n = ...
  library(tidyr)
})

set.seed(123)

# ---------- Clean gray theme ----------
theme_pub_gray <- function(base_size = 11, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      axis.title      = element_text(color = "grey20", margin = margin(t = 2)),
      axis.text       = element_text(color = "grey25"),
      panel.grid.minor= element_blank(),
      panel.grid.major= element_line(color = "grey85", linewidth = .3),
      axis.line       = element_line(color = "grey40", linewidth = .3),
      legend.position = "bottom",
      legend.title    = element_blank(),
      strip.background= element_rect(fill = "grey95", color = NA),
      strip.text      = element_text(face = "bold", color = "grey20")
    )
}

# ---------- Numeric helpers ----------
eps_num  <- 1e-9
clip <- function(x, lo, hi) pmax(lo, pmin(hi, x))
safe_div <- function(a, b, eps = eps_num) {
  b2 <- ifelse(is.finite(b) & abs(b) >= eps, b, ifelse(b >= 0, eps, -eps))
  a / b2
}

# =========================================================
# 1) Factors (parameters) and ranges
# =========================================================
factors <- c("delta", "sigma", "phi0", "phi1", "phi2", "phi3")
binf    <- c(-0.9,   0.10,   0.00,  0.01,  0.01,  0.01)
bsup    <- c( 0.2,   0.90,   0.05,  0.10,  0.10,  0.10)

# =========================================================
# 2) Outputs of interest (final values at t = 600)
# =========================================================
outputs <- c("yd", "ybp", "OSIN", "indebtY", "epsilon", "pi", "rdg")

# =========================================================
# 3) Model function: runs t=600 and returns final outputs
#    Variable parameters: delta (del), sigma, phi0..phi3
# =========================================================
run_model <- function(delta, sigma, phi0, phi1, phi2, phi3) {
  
  # ----- Fixed (baseline) parameters -----
  eps_G   <- 2.00
  eps_T   <- 0.78
  pi_G    <- 1.85
  pi_T    <- 1.00
  
  eta     <- -0.06
  psi     <- -0.57
  
  theta_max <- 0.5
  theta1    <- 1.0
  theta2    <- 10.0
  vi_max    <- 0.7
  vi1       <- 0.2
  vi2       <- 10.0
  
  lambda1 <- -0.015
  lambda2 <- -0.331
  lambda3 <- 0.113
  lambda4 <- -0.73
  lambda5 <- 0.30
  pd_bar  <- 0.0
  pf_bar  <- 0.000
  
  alpha1 <- 0.68
  alpha2 <- 0.27
  alpha3 <- 0.37
  alpha4 <- 0.16
  
  beta0 <- 0.95
  beta1 <- 4.0
  beta2 <- 0.01
  beta3 <- 0.03
  
  rho   <- 0.50
  tau1  <- 0.32
  tau2  <- 0.46
  z_ext <- 0.02
  
  # ----- Time dimension and allocation -----
  t <- 600
  
  et <- pdt <- pft <- epsilon <- pi <- theta <- vi <- rdg <- ye <- ydt <- ct <- numeric(t)
  PORT <- FDI <- port <- fdi <- embi <- OSIN <- osin <- gn <- numeric(t)
  it <- it_star <- gi <- gi_star <- numeric(t)
  Yt <- K <- Er <- Ex <- E <- R <- xt <- ybt <- numeric(t)
  ExDebtRes <- InDebt <- C <- numeric(t)
  
  KY_ratio <- InDebtY_ratio <- EY_ratio <- numeric(t)  # auxiliary ratios
  
  # ----- Initial conditions (t = 1) -----
  Yt[1]       <- 100.0
  K[1]        <- 200.0
  it_star[1]  <- 0.01
  it[1]       <- 0.013
  Ex[1]       <- 2.0
  Er[1]       <- 2.0
  FDI[1]      <- 2.0
  PORT[1]     <- 2.0
  ExDebtRes[1]<- 25.0
  InDebt[1]   <- 30.0
  E[1]        <- 15.0
  C[1]        <- FDI[1] + PORT[1]
  R[1]        <- E[1] + C[1]
  OSIN[1]     <- 0.6
  
  fdi[1]      <- 0.005
  gi[1]       <- 0
  gi_star[1]  <- -0.01
  rdg[1]      <- 0.0
  ybt[1]      <- 0.02
  epsilon[1]  <- 1.0
  pi[1]       <- 1.3
  pdt[1]      <- pd_bar
  et[1]       <- 0.0
  KY_ratio[1] <- K[1]/Yt[1]
  InDebtY_ratio[1] <- InDebt[1]/Yt[1]
  EY_ratio[1] <- E[1]/Yt[1]
  
  pft[] <- pf_bar
  
  # ----- Main loop -----
  for (i in 2:t) {
    
    # (A) Growth of i* and i (exogenous path for gi*)
    gi_star[i] <- -0.01
    gi[i]      <- delta * gi_star[i]
    
    # (B) Flows and R&D
    port[i] <- alpha1 * (gi[i-1] - gi_star[i-1] - embi[i-1])
    rdg[i]  <- beta0 * ye[i-1] + (1 - exp(-beta1 * rdg[i-1])) * fdi[i-1]
    fdi[i]  <- alpha2 * ybt[i-1] + alpha3 * rdg[i-1] - alpha4 * embi[i-1]
    
    # (C) Risk and domestic prices
    embi[i] <- phi0 + phi1 * osin[i-1] - phi2 * (epsilon[i-1] - pi[i-1]) + phi3 * gn[i-1]
    pdt[i]  <- pd_bar + lambda4 * gi_star[i-1] + lambda5 * et[i-1]
    
    # (D) Real exchange rate (lags)
    Er[i] <- (1 + et[i-1] + pft[i-1] - pdt[i-1]) * Er[i-1]
    
    # (E) Logistics (use rdg[i-1]*Er[i])
    vi[i]    <- vi_max    / (1 + vi1    * exp(-vi2    * (rdg[i-1] * Er[i])))
    theta[i] <- theta_max / (1 + theta1 * exp(-theta2 * (rdg[i-1] * Er[i])))
    
    # (F) Nominal exchange rate (lagged embi) and nominal exports
    et[i] <- lambda1 * (1 - theta[i]) * pdt[i] + lambda2 * port[i] + lambda3 * embi[i-1]
    et[i] <- clip(et[i], -0.95, 0.95)                 # avoid (1+et) <= 0
    Ex[i] <- (1 + et[i-1]) * Ex[i-1]
    
    # (G) Composite elasticities
    epsilon[i] <- theta[i] * eps_G + (1 - theta[i]) * eps_T
    pi[i]      <- (1 - vi[i]) * pi_G + vi[i] * pi_T
    
    # (H) Terms of trade and real exports (clamp 1+xt)
    xt[i] <- eta * (pdt[i] - pft[i] - et[i]) + epsilon[i] * z_ext
    f_E   <- 1 + xt[i]
    if (!is.finite(f_E) || f_E < eps_num) f_E <- eps_num
    E[i]  <- f_E * E[i-1]
    
    # (I) Financial stocks (before ybp) — safe weights
    FDI[i]  <- (1 + fdi[i])  * FDI[i-1]
    PORT[i] <- (1 + port[i]) * PORT[i-1]
    C[i]    <- FDI[i] + PORT[i]
    if (!is.finite(C[i]) || abs(C[i]) < eps_num) C[i] <- eps_num
    
    R[i]    <- E[i] + C[i]
    if (!is.finite(R[i]) || abs(R[i]) < eps_num) R[i] <- eps_num
    
    wF <- clip(FDI[i] / C[i], 0, 1)
    wP <- 1 - wF
    ct[i] <- wF * fdi[i] + wP * port[i]
    
    # (J) Thirlwall–Hussain (safe divisions)
    ybt[i] <- (
      (((E[i]/R[i]) * eta) + psi + 1) * (pdt[i] - et[i] - pft[i]) +
        (E[i]/R[i]) * epsilon[i] * z_ext +
        (C[i]/R[i]) * (ct[i] - pdt[i])
    ) / pi[i]
    
    # (K) Domestic and effective demand
    term_prod <- safe_div(rho * Yt[i-1], K[i-1])
    ydt[i] <- beta2 + beta3 * ( term_prod - it_star[i-1] * OSIN[i-1] - it[i-1] * safe_div(InDebt[i-1], Yt[i-1]) )
    ye[i]  <- min(ybt[i], ydt[i])
    
    # (L) Green structural change
    gn[i] <- tau1 * ybt[i-1] - tau2 * rdg[i-1]
    
    # (M) Debts and external vulnerability (uses 'sigma')
    ExDebtRes[i] <- ExDebtRes[i-1] * (1 + it_star[i-1]) + sigma * PORT[i]
    InDebt[i]    <- (1 + it[i-1]) * InDebt[i-1] + (1 - sigma) * PORT[i]
    
    denom_OSIN <- safe_div(Yt[i-1], E[i-1])
    OSIN[i]    <- ExDebtRes[i-1] / denom_OSIN
    if (!is.finite(OSIN[i])) OSIN[i] <- OSIN[i-1]
    osin[i]    <- safe_div(OSIN[i] - OSIN[i-1], OSIN[i-1])
    
    # (N) Output, capital and interest
    Yt[i]      <- (1 + ye[i])  * Yt[i-1]
    K[i]       <- (1 + ydt[i]) * K[i-1]
    it[i]      <- (1 + gi[i-1])     * it[i-1]
    it_star[i] <- (1 + gi_star[i-1]) * it_star[i-1]
    
    # Auxiliary ratios
    KY_ratio[i]       <- safe_div(K[i], Yt[i])
    InDebtY_ratio[i]  <- safe_div(InDebt[i], Yt[i])
    EY_ratio[i]       <- safe_div(E[i], Yt[i])
  }
  
  # Final outputs (t = 600)
  c(
    yd       = ydt[t],
    ybp      = ybt[t],
    OSIN     = OSIN[t],
    indebtY  = safe_div(InDebt[t], Yt[t]),
    epsilon  = epsilon[t],
    pi       = pi[t],
    rdg      = rdg[t]
  )
}

# =========================================================
# 4) Morris design and model runs
# =========================================================
morris_design <- morris(
  model   = NULL,
  factors = factors,
  r       = 10,                             # number of trajectories
  design  = list(type = "oat", levels = 4, grid.jump = 2),
  binf    = binf,
  bsup    = bsup
)

X <- data.matrix(morris_design$X)
n_runs <- nrow(X)
Y <- matrix(NA_real_, n_runs, length(outputs))
colnames(Y) <- outputs

# Run model for each X row
for (i in seq_len(n_runs)) {
  p <- X[i, ]
  yi <- run_model(delta = p[1], sigma = p[2], phi0 = p[3], phi1 = p[4], phi2 = p[5], phi3 = p[6])
  Y[i, ] <- yi[outputs]
}

# =========================================================
# 5) μ* (mu_star) and σ per output × factor
# =========================================================
sens_results <- data.frame()

for (j in seq_along(outputs)) {
  yj <- Y[, j]
  # robust imputation (median) for non-finite values
  if (any(!is.finite(yj))) {
    fin <- is.finite(yj)
    if (any(fin)) {
      med <- median(yj[fin]); yj[!fin] <- med
    } else {
      yj[] <- 0
    }
  }
  morris_output <- morris_design
  morris_output <- tell(morris_output, yj)
  EE <- morris_output$ee
  
  mu_star <- apply(EE, 2, function(x) mean(abs(x)))
  sigma_v <- apply(EE, 2, sd)
  
  sens_results <- rbind(
    sens_results,
    data.frame(
      output    = outputs[j],
      parameter = factors,
      mu_star   = as.numeric(mu_star),
      sigma     = as.numeric(sigma_v)
    )
  )
}

# =========================================================
# 6) Visualizations — titles removed; axes/legends in English
# =========================================================
theme_set(theme_pub_gray())

# (A) μ* vs σ
plot_data <- sens_results %>% filter(is.finite(mu_star), is.finite(sigma))

p <- ggplot(plot_data, aes(x = mu_star, y = sigma, label = parameter)) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, min.segment.length = 0) +
  facet_wrap(~ output, scales = "free") +
  labs(
    title = NULL,
    x = expression(paste(mu^"*", " (mean absolute EE)")),
    y = expression(paste(sigma, " (std. dev. of EE)"))
  )
print(p)

# Rebuild signed EEs (tidy)
ee_list <- vector("list", length(outputs)); names(ee_list) <- outputs
for (j in seq_along(outputs)) {
  yj <- Y[, j]
  if (any(!is.finite(yj))) {
    fin <- is.finite(yj); yj[!fin] <- if (any(fin)) median(yj[fin]) else 0
  }
  mo <- morris_design
  mo <- tell(mo, yj)
  EE <- as.data.frame(mo$ee); colnames(EE) <- factors
  ee_list[[j]] <- EE |>
    mutate(traj = dplyr::row_number(), output = outputs[j]) |>
    tidyr::pivot_longer(all_of(factors), names_to = "parameter", values_to = "ee")
}
ee_long <- dplyr::bind_rows(ee_list)

# Directional stats
dir_stats <- ee_long |>
  group_by(output, parameter) |>
  summarise(
    mu       = mean(ee, na.rm = TRUE),
    mu_star  = mean(abs(ee), na.rm = TRUE),
    sigma    = sd(ee, na.rm = TRUE),
    q05      = quantile(ee, 0.05, na.rm = TRUE),
    q50      = quantile(ee, 0.50, na.rm = TRUE),
    q95      = quantile(ee, 0.95, na.rm = TRUE),
    pos_share = mean(ee > 0, na.rm = TRUE),
    dir_consistency = ifelse(mu_star > 0, abs(mu)/mu_star, NA_real_),
    .groups = "drop"
  )

# (B) Signed bars (μ) + q05–q95
p_signed <- dir_stats |>
  ggplot(aes(x = reorder(parameter, mu), y = mu)) +
  geom_hline(yintercept = 0, linewidth = .3, color = "grey60") +
  geom_col(fill = "grey25") +
  geom_errorbar(aes(ymin = q05, ymax = q95), width = 0.2, linewidth = .4, color = "grey35") +
  coord_flip() +
  facet_wrap(~ output, scales = "free_x") +
  labs(title = NULL, x = NULL, y = expression(mu~"(signed mean of EE)"))
print(p_signed)

# (C) EE distributions (violin + box)
p_violins <- ee_long |>
  ggplot(aes(x = parameter, y = ee)) +
  geom_hline(yintercept = 0, linewidth = .3, color = "grey70") +
  geom_violin(trim = TRUE, scale = "width", fill = "grey85", color = "grey55") +
  geom_boxplot(width = 0.12, outlier.size = 0.6, fill = "grey35", color = "grey20") +
  coord_flip() +
  facet_wrap(~ output, scales = "free_x") +
  labs(title = NULL, x = NULL, y = "Elementary effects (EE)")
print(p_violins)

# (D) μ vs μ*
p_mu_mustar <- dir_stats |>
  ggplot(aes(x = mu, y = mu_star, label = parameter)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = .3) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = .3) +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 3, min.segment.length = 0) +
  facet_wrap(~ output, scales = "free") +
  labs(title = NULL,
       x = expression(mu~"(signed mean)"),
       y = expression(mu^"*"))
print(p_mu_mustar)

# (E) Heatmap of normalized μ
dir_stats_norm <- dir_stats |>
  group_by(output) |>
  mutate(mu_norm = mu / max(abs(mu), na.rm = TRUE)) |>
  ungroup()

p_heat <- dir_stats_norm |>
  ggplot(aes(x = parameter, y = output, fill = mu_norm)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish, name = "Normalized mu") +
  labs(title = NULL, x = NULL, y = NULL)
print(p_heat)

# (Optional) Ranking by |μ|
dir_ranking <- dir_stats |>
  mutate(abs_mu = abs(mu)) |>
  arrange(output, desc(abs_mu)) |>
  select(output, parameter, mu, mu_star, sigma, pos_share, dir_consistency) |>
  as_tibble()
print(dir_ranking, n = nrow(dir_ranking))

# =========================================================
# SAVE all plots to Desktop / "Graficos FGV final/baseline"
# =========================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a) || is.na(a) || length(a) == 0) b else a

get_desktop_path <- function(){
  candidatos <- unique(c(
    file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Área de Trabalho"),
    file.path(Sys.getenv("USERPROFILE"), "OneDrive", "Desktop"),
    file.path(Sys.getenv("USERPROFILE"), "Área de Trabalho"),
    file.path(Sys.getenv("USERPROFILE"), "Desktop"),
    path.expand("~/Desktop")
  ))
  alvo <- candidatos[dir.exists(candidatos)][1] %||% getwd()
  return(alvo)
}

OUT_DIR <- file.path(get_desktop_path(), "Graficos FGV final")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot_obj, nome, w = 12, h = 8, dpi = 300){
  stopifnot(inherits(plot_obj, "ggplot"))
  ggsave(filename = file.path(OUT_DIR, paste0(nome, ".png")),
         plot = plot_obj, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(filename = file.path(OUT_DIR, paste0(nome, ".pdf")),
         plot = plot_obj, width = w, height = h, device = "pdf", bg = "white")
}

save_plot(p,            "morris_mu_star_sigma")
save_plot(p_signed,     "morris_signed_mu")
save_plot(p_violins,    "morris_violins")
save_plot(p_mu_mustar,  "morris_mu_vs_mustar")
save_plot(p_heat,       "morris_heat_mu_norm")

message(sprintf("Plots saved to: %s", OUT_DIR))
