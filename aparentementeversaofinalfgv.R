# =========================================================
# Modelo (puro) – Equações corrigidas, prontos para rodar
# =========================================================
rm(list = ls())

# ---------------------------
# 1) Parâmetros (constantes)
# ---------------------------
eps_G   <- 2.00
eps_T   <- 0.78
pi_G    <- 1.85
pi_T    <- 1.00

eta     <- -0.06      # elasticidade-preço das exportações
psi     <- -0.57      # elasticidade-preço das importações

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
pd_bar  <- 0.0    # intercepto de pdt
pf_bar  <- 0.000      # pft constante

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
del <-0
# ---------------------------
# 2) Alocação e iniciais
# ---------------------------
t <- 600

et <- pdt <- pft <- epsilon <- pi <- theta <- vi <- rdg <- ye <- ydt <- ct <- numeric(t)
PORT <- FDI <- port <- fdi <- embi <- OSIN <- osin <- gn <- numeric(t)
it <- it_star <- gi <- gi_star <- numeric(t)
Yt <- K <- Er <- Ex <- E <- R <- xt <- ybt <- numeric(t)
ExDebtRes <- InDebt <- C <- numeric(t)
KY_ratio<- numeric(t)
InDebtY_ratio<- numeric(t)
EY_ratio<- numeric(t)

# Iniciais (t = 1)
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
gi[1]       <- 0.0
gi_star[1]  <- 0.0
rdg[1]      <- 0.0
ybt[1]      <- 0.02
epsilon[1]  <- 1.0
pi[1]       <- 1.3
pdt[1]      <- pd_bar
et[1]       <- 0.0
KY_ratio[1]<-K[1]/Yt[1]
InDebtY_ratio[1]<-InDebt[1]/Yt[1]
EY_ratio[1]<-E[1]/Yt[1]

# pft constante (vetor)
pft[] <- pf_bar

# ---------------------------
# 3) Loop principal
# ---------------------------
for (i in 2:t) {
  
  gi_star[i]<--0.00
  gi[i]<-del*gi_star[i]
  # (A) Fluxos e I&D
  port[i] <- alpha1 * (gi[i-1] - gi_star[i-1] - embi[i-1])
  rdg[i]  <- beta0 * ye[i-1] + (1 - exp(-beta1 * rdg[i-1])) * fdi[i-1]
  fdi[i]  <- alpha2 * ybt[i-1] + alpha3 * rdg[i-1] - alpha4 * embi[i-1]
  
  # (B) Risco e preços internos
  embi[i] <- phi0 + phi1 * osin[i-1] - phi2 * (epsilon[i-1] - pi[i-1]) + phi3 * gn[i-1]
  pdt[i]  <- pd_bar + lambda4 * gi_star[i-1] + lambda5 * et[i-1]
  
  # (C) Câmbio real (lags)
  Er[i] <- (1 + et[i-1] + pft[i-1] - pdt[i-1]) * Er[i-1]
  
  # (D) Logísticas (usam rdg[i-1]*Er[i])
  vi[i]    <- vi_max    / (1 + vi1    * exp(-vi2    * (rdg[i-1] * Er[i])))
  theta[i] <- theta_max / (1 + theta1 * exp(-theta2 * (rdg[i-1] * Er[i])))
  
  # (E) Câmbio nominal (embi defasado)
  et[i] <- lambda1 * (1 - theta[i]) * pdt[i] + lambda2 * port[i] + lambda3 * embi[i-1]
  Ex[i] <- (1 + et[i-1]) * Ex[i-1]
  
  # (F) Elasticidades compostas
  epsilon[i] <- theta[i] * eps_G + (1 - theta[i]) * eps_T
  pi[i]      <- (1 - vi[i]) * pi_G + vi[i] * pi_T
  
  # (G) Termos de troca e exportações reais
  xt[i] <- eta * (pdt[i] - pft[i] - et[i]) + epsilon[i] * z_ext
  E[i]  <- (1 + xt[i]) * E[i-1]
  
  # (H) Estoques financeiros (antes de ybp)
  FDI[i]  <- (1 + fdi[i])  * FDI[i-1]
  PORT[i] <- (1 + port[i]) * PORT[i-1]
  C[i]    <- FDI[i] + PORT[i]
  R[i]    <- E[i] + C[i]
  ct[i]   <- (FDI[i]/C[i]) * fdi[i] + (PORT[i]/C[i]) * port[i]
  
  # (I) Thirlwall–Hussain (forma corrigida)
  ybt[i] <- (
    (((E[i]/R[i]) * eta) + psi + 1) * (pdt[i] - et[i] - pft[i]) +
      (E[i]/R[i]) * epsilon[i] * z_ext +
      (C[i]/R[i]) * (ct[i] - pdt[i])
  ) / pi[i]
  
  # (J) Demanda doméstica e efetiva
  ydt[i] <- beta2 + beta3 * (
    (rho * Yt[i-1]) / K[i-1] -
      it_star[i-1] * OSIN[i-1] -
      it[i-1] * (InDebt[i-1] / Yt[i-1])
  )
  ye[i] <- min(ybt[i], ydt[i])
  
  # (K) Mudança estrutural verde (sinal correto)
  gn[i] <- tau1 * ybt[i-1] - tau2 * rdg[i-1]
  
  # (L) Dívidas e vulnerabilidade externa
  ExDebtRes[i] <- ExDebtRes[i-1] * (1 + it_star[i-1]) + sigma * PORT[i]
  InDebt[i]    <- (1 + it[i-1]) * InDebt[i-1] + (1 - sigma) * PORT[i]
  OSIN[i]      <- ExDebtRes[i-1] / (Yt[i-1] / Ex[i-1])
  osin[i]      <- (OSIN[i] - OSIN[i-1]) / OSIN[i-1]
  
  # (M) Produto, capital e juros
  Yt[i]      <- (1 + ye[i])  * Yt[i-1]
  K[i]       <- (1 + ydt[i]) * K[i-1]
  it[i]      <- (1 + gi[i-1])     * it[i-1]
  it_star[i] <- (1 + gi_star[i-1]) * it_star[i-1]
  
  
  KY_ratio[i]<-K[i]/Yt[i]
  InDebtY_ratio[i]<-InDebt[i]/Yt[i]
  EY_ratio[i]<-E[i]/Yt[i]
  
}

# ---------------------------
# 4) Saída organizada
# ---------------------------
resultado <- data.frame(
  t = 0:(t-1),
  ye = ye, yd = ydt, ybp = ybt, xt = xt,
  eps = epsilon, pi = pi,
  theta = theta, vi = vi,
  R = R, E = E, C = C,
  FDI = FDI, PORT = PORT,
  OSIN = OSIN, ExDebtRes = ExDebtRes, InDebt = InDebt,
  Y = Yt, K = K,
  i_star = it_star, i = it,
  Er = Er, Ex = Ex,
  pdt = pdt, et = et, pft = pft,
  ct = ct, fdi = fdi, port = port, embi = embi,
  gi = gi, gi_star = gi_star, gn = gn, rdg=rdg, osin=osin, ky_ratio=KY_ratio, IndebtY_ratio=InDebtY_ratio, EY_ratio=EY_ratio
)

# ================================
# 5) GRÁFICOS — CINZA, LIMPOS, t>=3
# ================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(scales); library(patchwork)
})

# ---- filtro temporal
start_period <- 1
df <- resultado |>
  mutate(time = t) |>
  filter(time >= start_period)

# ---- tema de publicação em cinza (inspirado em minimal-ink + Nature)
theme_pub_gray <- function(base_size = 11, base_family = "sans"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(face = "bold", color = "grey10", margin = margin(b = 5)),
      plot.subtitle   = element_text(color = "grey35"),
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

# ---- helpers
xb <- scale_x_continuous(breaks = pretty_breaks(8), expand = expansion(mult = c(.01,.02)))
yg <- scale_y_continuous(labels = label_number())

# Linhas sempre em cinza-escuro; séries múltiplas com dois tons
col_single <- "grey20"
col_a      <- "grey20"
col_b      <- "grey55"

plot_single <- function(var, title){
  ggplot(df, aes(time, .data[[var]])) +
    geom_line(linewidth = 1.15, lineend = "round", color = col_single) +
    xb + yg + labs(title = title, x = NULL, y = NULL) +
    theme_pub_gray()
}

plot_dual <- function(var1, var2, lab1, lab2, title = NULL){
  dd <- df |>
    transmute(time, !!lab1 := .data[[var1]], !!lab2 := .data[[var2]]) |>
    pivot_longer(-time, names_to = "serie", values_to = "valor")
  ggplot(dd, aes(time, valor, color = serie)) +
    geom_line(linewidth = 1.15, lineend = "round") +
    scale_color_manual(values = setNames(c(col_a, col_b), c(lab1, lab2))) +
    xb + yg + labs(title = title, x = NULL, y = NULL) +
    theme_pub_gray()
}

# -------- gráficos (mesma estrutura dos seus prints)
p_ye       <- plot_single("ye", "ye")
p_y_vs_ybp <- plot_dual("ye", "ybp", "ye", "ybp")

# epsilon t (eps) x pi t (pi)
p_eps_pi <- df |>
  transmute(time, `epsilon t` = eps, `pit` = pi) |>
  pivot_longer(-time, names_to = "serie", values_to = "valor") |>
  ggplot(aes(time, valor, color = serie)) +
  geom_line(linewidth = 1.15, lineend = "round") +
  scale_color_manual(values = c("epsilon t" = col_a, "pit" = col_b)) +
  xb + yg + labs(title = NULL, x = NULL, y = NULL) +
  theme_pub_gray()

# K/Y
p_KY   <- plot_single("ky_ratio", "K/Y")

# gn (com linha-zero discreta)
p_gn <- ggplot(df, aes(time, gn)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = .35) +
  geom_line(linewidth = 1.15, lineend = "round", color = col_single) +
  xb + yg + labs(title = "gn", x = NULL, y = NULL) +
  theme_pub_gray()

# rdg
p_rdg  <- plot_single("rdg", "rdg")

# OSIN
p_OSIN <- plot_single("OSIN", "OSIN")

# indebt/Y  (atenção: nome no data.frame = IndebtY_ratio)
p_InDebtY <- plot_single("IndebtY_ratio", "indebt/Y")

# -------- painel final (proporções equilibradas)
painel <- (p_ye | p_y_vs_ybp) /
  (p_eps_pi | p_KY) /
  (p_gn | p_rdg) /
  (p_OSIN | p_InDebtY)

painel

# (opcional) exporte em alta
# dir.create("plots", showWarnings = FALSE)
# ggsave("plots/painel_cinza.png", painel, width = 14, height = 12, dpi = 300)
# ggsave("plots/ye_cinza.png", p_ye, width = 7, height = 4.5, dpi = 300)
