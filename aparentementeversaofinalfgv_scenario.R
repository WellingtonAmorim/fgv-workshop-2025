# =========================================================
# SCENARIO VERSION
# Pure model – corrected equations, ready to run
# =========================================================
rm(list = ls())

# ---------------------------
# 1) Parameters (constants)
# ---------------------------
eps_G   <- 2.00
eps_T   <- 0.78
pi_G    <- 1.85
pi_T    <- 1.00

eta     <- -0.06      # price elasticity of exports
psi     <- -0.57      # price elasticity of imports

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
pd_bar  <- 0.0    # intercept of pdt
pf_bar  <- 0.000  # constant pft

alpha1 <- 0.68
alpha2 <- 0.27
alpha3 <- 0.37
alpha4 <- 0.16

phi0 <- 0.01
phi1 <- 0.05
phi2 <- 0.02
phi3 <- 0.01

beta0 <- 0.95
beta1 <- 4.0
beta2 <- 0.01
beta3 <- 0.03

rho   <- 0.50
sigma <- 0.50
tau1  <- 0.32
tau2  <- 0.46
z_ext <- 0.02
del   <- 0

# ---------------------------
# 2) Allocation and initials
# ---------------------------
t <- 600

et <- prices <- pdt <- pft <- epsilon <- pi <- theta <- vi <- rdg <- ye <- ydt <- ct <-
  effect_one <- effect_two <- effect_three <- numeric(t)

PORT <- FDI <- port <- fdi <- embi <- OSIN <- osin <- gn <- numeric(t)
it <- it_star <- gi <- gi_star <- numeric(t)
Yt <- K <- Er <- Ex <- E <- R <- xt <- ybt <- numeric(t)
ExDebtRes <- InDebt <- C <- numeric(t)
KY_ratio <- numeric(t)
InDebtY_ratio <- numeric(t)
EY_ratio <- numeric(t)

# Initial conditions (t = 1)
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

effect_one[1]   <- 0
effect_two[1]   <- 0
effect_three[1] <- 0

# constant pft (vector)
pft[] <- pf_bar

# ---------------------------
# 3) Main loop
# ---------------------------
for (i in 2:t) {
  
  gi_star[i] <- -0.01
  gi[i]      <- del * gi_star[i]
  
  # (A) Flows and R&D
  port[i] <- alpha1 * (gi[i-1] - gi_star[i-1] - embi[i-1])
  rdg[i]  <- beta0 * ye[i-1] + (1 - exp(-beta1 * rdg[i-1])) * fdi[i-1]
  fdi[i]  <- alpha2 * ybt[i-1] + alpha3 * rdg[i-1] - alpha4 * embi[i-1]
  
  # (B) Risk and domestic prices
  embi[i] <- phi0 + phi1 * osin[i-1] - phi2 * (epsilon[i-1] - pi[i-1]) + phi3 * gn[i-1]
  pdt[i]  <- pd_bar + lambda4 * gi_star[i-1] + lambda5 * et[i-1]
  
  # (C) Real exchange rate (lags)
  Er[i] <- (1 + et[i-1] + pft[i-1] - pdt[i-1]) * Er[i-1]
  
  # (D) Logistic functions (use rdg[i-1]*Er[i])
  vi[i]    <- vi_max    / (1 + vi1    * exp(-vi2    * (rdg[i-1] * Er[i])))
  theta[i] <- theta_max / (1 + theta1 * exp(-theta2 * (rdg[i-1] * Er[i])))
  
  # (E) Nominal exchange rate (lagged embi)
  et[i] <- lambda1 * (1 - theta[i]) * pdt[i] + lambda2 * port[i] + lambda3 * embi[i-1]
  Ex[i] <- (1 + et[i-1]) * Ex[i-1]
  
  # (F) Composite elasticities
  epsilon[i] <- theta[i] * eps_G + (1 - theta[i]) * eps_T
  pi[i]      <- (1 - vi[i]) * pi_G + vi[i] * pi_T
  
  # (G) Terms of trade and real exports
  xt[i] <- eta * (pdt[i] - pft[i] - et[i]) + epsilon[i] * z_ext
  E[i]  <- (1 + xt[i]) * E[i-1]
  
  # (H) Financial stocks (before ybp)
  FDI[i]  <- (1 + fdi[i])  * FDI[i-1]
  PORT[i] <- (1 + port[i]) * PORT[i-1]
  C[i]    <- FDI[i] + PORT[i]
  R[i]    <- E[i] + C[i]
  ct[i]   <- (FDI[i]/C[i]) * fdi[i] + (PORT[i]/C[i]) * port[i]
  
  # (I) Thirlwall–Hussain (corrected form)
  ybt[i] <- (
    (((E[i]/R[i]) * eta) + psi + 1) * (pdt[i] - et[i] - pft[i]) +
      (E[i]/R[i]) * epsilon[i] * z_ext +
      (C[i]/R[i]) * (ct[i] - pdt[i])
  ) / pi[i]
  
  # (J) Domestic and effective demand
  ydt[i] <- beta2 + beta3 * (
    (rho * Yt[i-1]) / K[i-1] -
      it_star[i-1] * OSIN[i-1] -
      it[i-1] * (InDebt[i-1] / Yt[i-1])
  )
  ye[i] <- min(ybt[i], ydt[i])
  
  # (K) Green structural change
  gn[i] <- tau1 * ybt[i-1] - tau2 * rdg[i-1]
  
  # (L) Debt and external vulnerability
  ExDebtRes[i] <- ExDebtRes[i-1] * (1 + it_star[i-1]) + sigma * PORT[i]
  InDebt[i]    <- (1 + it[i-1]) * InDebt[i-1] + (1 - sigma) * PORT[i]
  OSIN[i]      <- ExDebtRes[i-1] / (Yt[i-1] / Ex[i-1])
  osin[i]      <- (OSIN[i] - OSIN[i-1]) / OSIN[i-1]
  
  # (M) Output, capital and interest rates
  Yt[i]      <- (1 + ye[i])  * Yt[i-1]
  K[i]       <- (1 + ydt[i]) * K[i-1]
  it[i]      <- (1 + gi[i-1])      * it[i-1]
  it_star[i] <- (1 + gi_star[i-1]) * it_star[i-1]
  
  KY_ratio[i]      <- K[i]/Yt[i]
  InDebtY_ratio[i] <- InDebt[i]/Yt[i]
  EY_ratio[i]      <- E[i]/Yt[i]
  prices[i]        <- pdt[i] + pft[i] + et[i]
  
  # (N) Decomposition of the BOP-constrained growth rate
  effect_one[i]   <- ((E[i]/R[i]) * epsilon[i] * z_ext) / pi[i]
  effect_two[i]   <- ((C[i]/R[i]) * (ct[i] - pdt[i])) / pi[i]
  effect_three[i] <- ((((E[i]/R[i]) * eta) + psi + 1) * (pdt[i] - et[i] - pft[i])) / pi[i]
}

# ---------------------------
# 4) Organized output
# ---------------------------
resultado <- data.frame(
  t = 0:(t-1),
  ye = ye, yd = ydt, ybp = ybt, xt = xt,
  effect_one   = effect_one,
  effect_two   = effect_two,
  effect_three = effect_three,
  eps = epsilon, pi = pi,
  theta = theta, vi = vi,
  R = R, E = E, C = C,
  FDI = FDI, PORT = PORT,
  OSIN = OSIN, ExDebtRes = ExDebtRes, InDebt = InDebt,
  Y = Yt, K = K,
  i_star = it_star, i = it,
  Er = Er, Ex = Ex,
  pdt = pdt, et = et, pft = pft,
  prices = prices,
  ct = ct, fdi = fdi, port = port, embi = embi,
  gi = gi, gi_star = gi_star, gn = gn, rdg = rdg, osin = osin,
  ky_ratio = KY_ratio, IndebtY_ratio = InDebtY_ratio, EY_ratio = EY_ratio
)

# ================================
# 5) PLOTS — single, dual and multi-series, neutral scientific theme
# ================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales)
})

# ---- time filter
start_period <- 1
df <- resultado |>
  mutate(time = t) |>
  filter(time >= start_period)

# ---- neutral scientific theme (white background, minimal color)
# larger base_size to make axis numbers readable in documents
theme_pub_sci <- function(base_size = 16, base_family = "sans"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.4),
      panel.grid.major = element_line(color = "grey90", linewidth = .3),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(color = "grey20",
                                      size  = base_size,
                                      margin = margin(t = 4)),
      axis.text        = element_text(color = "grey25",
                                      size  = base_size * 0.9),
      axis.ticks       = element_line(color = "grey50", linewidth = .3),
      plot.title       = element_text(face = "bold", color = "grey15",
                                      size  = base_size * 1.1,
                                      margin = margin(b = 6), hjust = 0),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = base_size * 0.9)
    )
}

# ---- scales helpers
xb <- scale_x_continuous(
  breaks = pretty_breaks(8),
  expand = expansion(mult = c(.01, .02))
)
yg <- scale_y_continuous(labels = label_number())

# ---- color palette (neutral, publication-ready)
col_single <- "#1F3B4D"  # dark blue
col_a      <- col_single
col_b      <- "#8C8C8C"  # neutral grey
col_c      <- "#595959"  # slightly lighter grey

# -----------------------------
# 5.1 Generic plotting functions
# -----------------------------

# single time series (with legend, first letter capitalized)
plot_single <- function(var_name, label = NULL){
  if (is.null(label)) label <- var_name
  
  # capitalize first letter of legend text
  label <- paste0(
    toupper(substr(label, 1, 1)),
    substr(label, 2, nchar(label))
  )
  
  ggplot(df, aes(x = time)) +
    geom_line(aes(y = .data[[var_name]], color = label),
              linewidth = 1.4, lineend = "round") +
    scale_color_manual(
      values = setNames(col_single, label)
    ) +
    xb + yg +
    labs(
      title = NULL,
      x     = "Time (periods)",
      y     = NULL,
      color = NULL
    ) +
    theme_pub_sci()
}

# two time series in the same panel (with legend)
plot_dual_ts <- function(var1, var2,
                         lab1 = NULL, lab2 = NULL){
  if (is.null(lab1)) lab1 <- var1
  if (is.null(lab2)) lab2 <- var2
  
  ggplot(df, aes(x = time)) +
    geom_line(aes(y = .data[[var1]], color = lab1),
              linewidth = 1.4, lineend = "round") +
    geom_line(aes(y = .data[[var2]], color = lab2),
              linewidth = 1.4, lineend = "round") +
    scale_color_manual(
      values = setNames(c(col_a, col_b), c(lab1, lab2))
    ) +
    xb + yg +
    labs(
      title = NULL,
      x     = "Time (periods)",
      y     = NULL,
      color = NULL
    ) +
    theme_pub_sci()
}

# -----------------------------
# 5.2 Single-series plots for all variables
# -----------------------------
vars_to_plot <- setdiff(names(resultado), "t")

plots <- lapply(vars_to_plot, function(v) plot_single(v))
names(plots) <- vars_to_plot

# -----------------------------
# 5.3 Specific dual plots (with descriptive labels)
# -----------------------------

# BOP-constrained vs effective growth
p_ybp_ye <- plot_dual_ts(
  var1 = "ybp",
  var2 = "ye",
  lab1 = "BOP-constrained growth (ybp)",
  lab2 = "Effective demand growth (ye)"
)

# NEW: BOP-constrained vs domestic demand growth
p_ybp_yd <- plot_dual_ts(
  var1 = "ybp",
  var2 = "yd",
  lab1 = "BOP-constrained growth (ybp)",
  lab2 = "Domestic demand growth (yd)"
)

# Real vs nominal exchange rate
p_Er_Ex <- plot_dual_ts(
  var1 = "Er",
  var2 = "Ex",
  lab1 = "Real exchange rate (Er)",
  lab2 = "Nominal exchange rate (Ex)"
)

# Elasticities epsilon and pi
p_eps_pi <- plot_dual_ts(
  var1 = "eps",
  var2 = "pi",
  lab1 = "Elasticity epsilon (eps)",
  lab2 = "Elasticity pi (pi)"
)

# -----------------------------
# 5.4 Three-effects decomposition in a single panel
# -----------------------------

p_3effects <- ggplot(df, aes(x = time)) +
  geom_line(
    aes(y = effect_one,   color = "Effect 1 – elasticity term"),
    linewidth = 1.2, lineend = "round"
  ) +
  geom_line(
    aes(y = effect_two,   color = "Effect 2 – composition term"),
    linewidth = 1.2, lineend = "round"
  ) +
  geom_line(
    aes(y = effect_three, color = "Effect 3 – relative-price term"),
    linewidth = 1.2, lineend = "round"
  ) +
  scale_color_manual(
    values = c(
      "Effect 1 – elasticity term"      = col_a,
      "Effect 2 – composition term"     = col_b,
      "Effect 3 – relative-price term"  = col_c
    )
  ) +
  xb + yg +
  labs(
    title = NULL,
    x     = "Time (periods)",
    y     = NULL,
    color = NULL
  ) +
  theme_pub_sci()

# -----------------------------
# 5.5 Save all plots to 4K folder "FGV_final_plots_4k"
# -----------------------------
plots_dir_4k <- "FGV_final_plots_4k"
dir.create(plots_dir_4k, showWarnings = FALSE)

## 5.5.1 Single-series plots (all variables from `resultado`)
for (i in seq_along(vars_to_plot)) {
  v <- vars_to_plot[i]
  
  # torna o nome "seguro" para arquivo
  file_safe <- gsub("[^A-Za-z0-9_]+", "_", v)
  
  # inclui um índice na frente para evitar conflito entre osin / OSIN em sistemas case-insensitive
  filename  <- sprintf("%03d_%s.png", i, file_safe)
  
  ggsave(
    filename = file.path(plots_dir_4k, filename),
    plot     = plots[[v]],
    width    = 3840, height = 2160, units = "px", dpi = 300
  )
}
## 5.5.2 Specific dual plots: ybp × ye, ybp × yd, Er × Ex, eps × pi
ggsave(
  filename = file.path(plots_dir_4k, "ybp_ye_dual.png"),
  plot     = p_ybp_ye,
  width    = 3840, height = 2160, units = "px", dpi = 300
)

ggsave(
  filename = file.path(plots_dir_4k, "ybp_yd_dual.png"),
  plot     = p_ybp_yd,
  width    = 3840, height = 2160, units = "px", dpi = 300
)

ggsave(
  filename = file.path(plots_dir_4k, "Er_Ex_dual.png"),
  plot     = p_Er_Ex,
  width    = 3840, height = 2160, units = "px", dpi = 300
)

ggsave(
  filename = file.path(plots_dir_4k, "eps_pi_dual.png"),
  plot     = p_eps_pi,
  width    = 3840, height = 2160, units = "px", dpi = 300
)

## 5.5.3 Three-effects decomposition plot
ggsave(
  filename = file.path(plots_dir_4k, "three_effects.png"),
  plot     = p_3effects,
  width    = 3840, height = 2160, units = "px", dpi = 300
)
