# ===========================================================
# flint.R
# Flexible interactions with learned powers — comparing:
#  (A) binary switches (zeta) vs (B) horseshoe on kappa
# Datasets: Sim A/B/C/D, ToothGrowth (Gaussian), warpbreaks (NegBin),
#           PlantGrowth (Gaussian), InsectSprays (Poisson)
# Notes:
#  - Interactions use powered covariates: f(x) = sign(x)*|x|^power
#  - No p-values; only PPCs, RMSE/MAE, Rhat, ESS, Heidelberger
# ===========================================================
suppressPackageStartupMessages({
  library(nimble)
  library(coda)
})
set.seed(123)

# ------------------------- helpers -------------------------
zscore <- function(x) as.numeric(scale(x))
spow   <- function(x, p, eps=1e-6) sign(x) * (abs(x)+eps)^p

build_pairs <- function(P) {
  if (P < 2) return(list(pair1=integer(0), pair2=integer(0)))
  cmb <- t(combn(P, 2))
  list(pair1 = cmb[,1], pair2 = cmb[,2])
}

ppc_gaussian <- function(S, y, prefix="y_rep[") {
  yc <- S[, startsWith(colnames(S), prefix), drop=FALSE]
  yhat <- rowMeans(yc)
  c(obs_mean = mean(y), rep_mean = mean(yhat),
    obs_sd   = sd(y),   rep_sd   = mean(apply(yc, 1, sd)),
    RMSE = sqrt(mean((y - yhat)^2)), MAE = mean(abs(y - yhat)))
}
ppc_count <- function(S, y, prefix="y_rep[") {
  yc <- S[, startsWith(colnames(S), prefix), drop=FALSE]
  yhat <- rowMeans(yc)
  c(obs_mean = mean(y), rep_mean = mean(yhat),
    obs_var  = var(y),  rep_var  = mean(apply(yc, 1, var)),
    RMSE = sqrt(mean((y - yhat)^2)), MAE = mean(abs(y - yhat)))
}

diag_summary <- function(mcmc_list) {
  gd <- try(gelman.diag(mcmc_list, autoburnin = FALSE), silent=TRUE)
  rhat_max <- if (inherits(gd, "try-error")) NA_real_ else suppressWarnings(max(gd$psrf[,1], na.rm=TRUE))
  ess <- try(effectiveSize(mcmc_list), silent=TRUE)
  ess_min <- if (inherits(ess, "try-error")) NA_real_ else suppressWarnings(min(ess, na.rm=TRUE))
  heidel_pass <- try({
    chsum <- lapply(mcmc_list, function(ch) {
      hd <- heidel.diag(ch)
      mean(hd$stest[,"stationary"] == 1, na.rm = TRUE)
    })
    mean(unlist(chsum))
  }, silent=TRUE)
  heidel_pass <- if (inherits(heidel_pass, "try-error")) NA_real_ else heidel_pass
  c(Rhat_max = rhat_max, ESS_min = ess_min, Heidel_pass = heidel_pass)
}

top_kappa_table <- function(S, kappa_prefix="kappa[", top=6) {
  keep <- grepl("^kappa\\[[0-9]+\\]$", colnames(S))
  if (!any(keep)) return(NULL)
  K  <- S[, keep, drop=FALSE]
  md <- apply(abs(K), 2, median)
  ord <- order(md, decreasing = TRUE)[seq_len(min(top, length(md)))]
  out <- cbind(kappa = names(md)[ord],
               median = round(md[ord], 3),
               lo     = round(apply(K[,ord,drop=FALSE], 2, function(v) quantile(v, .025)), 3),
               hi     = round(apply(K[,ord,drop=FALSE], 2, function(v) quantile(v, .975)), 3))
  rownames(out) <- NULL
  out
}

summ_ci <- function(x) c(mean=mean(x), lo=quantile(x,0.025), hi=quantile(x,0.975))

# ----------------- nimbleCode generator --------------------
# Interactions are products of powered covariates:
#  int_term[i,m] = MULT[m] * kappa[m] * (sign(x_j1)*|x_j1|^psi1[m]) * (sign(x_j2)*|x_j2|^psi2[m])
# j1/j2 are constant index maps (passed in via constants), so no node redefinitions.
code_flex <- function(N, P, Q, family = c("poisson","gaussian","nbinom"),
                      prior = c("switch","horseshoe")) {
  family <- match.arg(family)
  prior  <- match.arg(prior)
  nimbleCode({
    for (i in 1:N) {
      # mains: signed-power contributions
      for (j in 1:P) {
        s[i,j] <- 2*step(xz[i,j]) - 1
        a[i,j] <- abs(xz[i,j]) + eps
        main_contrib[i,j] <- beta[j] * s[i,j] * pow(a[i,j], gamma[j])
      }
      if (P > 1) { lin_main[i] <- sum(main_contrib[i, 1:P]) } else { lin_main[i] <- main_contrib[i,1] }

      # interactions (all pairs) using constant index maps j1/j2
      for (m in 1:Q) {
        f1[i,m] <- s[i, j1[m]] * pow(a[i, j1[m]], psi1[m])
        f2[i,m] <- s[i, j2[m]] * pow(a[i, j2[m]], psi2[m])
        int_term[i,m] <- MULT[m] * kappa[m] * f1[i,m] * f2[i,m]
      }
      if (Q > 1) { lin_int[i] <- sum(int_term[i, 1:Q]) } else { lin_int[i] <- int_term[i,1] }

      # likelihood per family
      if (FAM == 1) {
        log(mu[i]) <- alpha + lin_main[i] + lin_int[i]
        y[i] ~ dpois(mu[i]); y_rep[i] ~ dpois(mu[i])
      }
      if (FAM == 2) {
        mu_g[i] <- alpha + lin_main[i] + lin_int[i]
        y[i] ~ dnorm(mu_g[i], prec); y_rep[i] ~ dnorm(mu_g[i], prec)
      }
      if (FAM == 3) {
        log(mu_nb[i]) <- alpha + lin_main[i] + lin_int[i]
        y[i] ~ dnbinom(size=delta, prob = delta / (delta + mu_nb[i]))
        y_rep[i] ~ dnbinom(size=delta, prob = delta / (delta + mu_nb[i]))
      }
    }

    # priors: mains & powers
    alpha ~ dnorm(0, 1.0E-4)
    for (j in 1:P) { beta[j] ~ dnorm(0, prec_beta); gamma[j] ~ dunif(0.25, 3) }
    for (m in 1:Q) { psi1[m] ~ dunif(0.25, 3); psi2[m] ~ dunif(0.25, 3) }

    # switch or horseshoe on kappa
    if (USE_SWITCH == 1) {
      for (m in 1:Q) {
        zeta[m] ~ dbern(pi_int);  MULT[m] <- zeta[m]
        kappa[m] ~ dnorm(0, prec_kappa)
      }
      pi_int ~ dbeta(1, 9);  prec_kappa ~ dgamma(0.5, 0.5)
    }
    if (USE_HS == 1) {
      tau2 ~ dinvgamma(0.5, 1 / xi);  xi ~ dinvgamma(0.5, 1)
      for (m in 1:Q) {
        lambda2[m] ~ dinvgamma(0.5, 1 / nu[m]);  nu[m] ~ dinvgamma(0.5, 1)
        prec_kappa_m[m] <- 1 / (tau2 * lambda2[m])
        kappa[m] ~ dnorm(0, prec_kappa_m[m])
        MULT[m] <- 1
      }
    }
    prec_beta ~ dgamma(0.5, 0.5)
    if (FAM == 2) prec ~ dgamma(1,1)
    if (FAM == 3) delta ~ dgamma(0.5,0.5)

    # overall interaction magnitude
    kappa_mag <- sum(pow(kappa[1:Q], 2))
  })
}

# --------------- generic fit wrapper -----------------------
run_flex <- function(y, XZ, family = c("poisson","gaussian","nbinom"),
                     prior = c("switch","horseshoe"),
                     niter=6000, nburn=3000, thin=2, nchains=2) {
  family <- match.arg(family); prior <- match.arg(prior)
  fam_code <- switch(family, poisson=1L, gaussian=2L, nbinom=3L)
  N <- nrow(XZ); P <- ncol(XZ)
  pairs <- build_pairs(P); Q <- length(pairs$pair1)  # P>=2 here, so Q>=1

  code <- code_flex(N,P,Q,family,prior)
  consts <- list(N=N, P=P, Q=Q, eps=1e-6,
                 j1 = as.integer(pairs$pair1),
                 j2 = as.integer(pairs$pair2),
                 FAM = fam_code,
                 USE_SWITCH = as.integer(prior=="switch"),
                 USE_HS     = as.integer(prior=="horseshoe"))
  data <- list(y = y, xz = as.matrix(XZ))

  init_fun <- function() {
    ini <- list(alpha = 0,
                beta  = rnorm(P,0,0.5),
                gamma = runif(P,0.8,1.2),
                psi1  = runif(Q,0.8,1.2),
                psi2  = runif(Q,0.8,1.2),
                prec_beta = 1)
    if (prior == "switch") {
      ini$pi_int <- 0.15;  ini$zeta <- rbinom(Q,1,0.1)
      ini$kappa <- rnorm(Q,0,0.1); ini$prec_kappa <- 1
    } else {
      ini$tau2 <- 1; ini$xi <- 1; ini$lambda2 <- rep(1, Q); ini$nu <- rep(1, Q)
      ini$kappa <- rnorm(Q,0,0.1)
    }
    if (family == "gaussian") ini$prec <- 1/var(y)
    if (family == "nbinom")   ini$delta <- 1
    ini
  }

  model <- nimbleModel(code, constants=consts, data=data, inits=init_fun())
  monitors <- c("alpha",
                paste0("beta[",1:P,"]"),
                paste0("gamma[",1:P,"]"),
                paste0("psi1[",1:Q,"]"),
                paste0("psi2[",1:Q,"]"),
                paste0("kappa[",1:Q,"]"),
                "kappa_mag",
                paste0("y_rep[",1:N,"]"))
  if (prior == "switch")    monitors <- c(monitors, "pi_int")
  if (prior == "horseshoe") monitors <- c(monitors, "tau2")

  conf <- configureMCMC(model, monitors = monitors)
  # slice sampling for power parameters
  for (nm in c(paste0("gamma[",1:P,"]"), paste0("psi1[",1:Q,"]"), paste0("psi2[",1:Q,"]"))) {
    if (nm %in% model$getNodeNames()) { conf$removeSamplers(nm); conf$addSampler(nm, type="slice") }
  }
  mcmc   <- buildMCMC(conf)
  cmodel <- compileNimble(model)
  cmcmc  <- compileNimble(mcmc, project=model)

  tm <- system.time({
    samps <- runMCMC(cmcmc, niter=niter, nburnin=nburn, thin=thin,
                     nchains=nchains, samplesAsCodaMCMC=TRUE, progressBar=FALSE)
  })
  S <- as.matrix(do.call(rbind, samps))
  list(samples = samps, Smat = S, time = unname(tm["elapsed"]),
       P=P, Q=Q, prior=prior, family=family)
}

report_fit <- function(name, y, fit) {
  cat("\n--- ", name, " [", fit$prior, "] ---\n", sep="")
  fam <- fit$family; S <- fit$Smat
  ppc <- if (fam == "gaussian") ppc_gaussian(S, y) else ppc_count(S, y)
  di  <- diag_summary(fit$samples)
  cat(sprintf("Runtime (s): %.2f | Rhat_max: %.3f | ESS_min: %.0f | Heidel_pass: %.2f\n",
              fit$time, di["Rhat_max"], di["ESS_min"], di["Heidel_pass"]))
  if (fam == "gaussian") {
    cat(sprintf("PPC: obs_mean=%.3f rep_mean=%.3f | obs_sd=%.3f rep_sd=%.3f | RMSE=%.3f MAE=%.3f\n",
                ppc["obs_mean"], ppc["rep_mean"], ppc["obs_sd"], ppc["rep_sd"], ppc["RMSE"], ppc["MAE"]))
  } else {
    cat(sprintf("PPC: obs_mean=%.3f rep_mean=%.3f | obs_var=%.3f rep_var=%.3f | RMSE=%.3f MAE=%.3f\n",
                ppc["obs_mean"], ppc["rep_mean"], ppc["obs_var"], ppc["rep_var"], ppc["RMSE"], ppc["MAE"]))
  }
  if ("kappa_mag" %in% colnames(S)) {
    km <- S[,"kappa_mag"]; cat("||kappa||^2 ~ ", paste(round(summ_ci(km),3), collapse=" "), "\n", sep="")
  }
  tk <- top_kappa_table(S); if (!is.null(tk)) { cat("Top |kappa| (median [2.5%,97.5%]):\n"); print(tk, row.names=FALSE) }
}

# ------------------- DATA PREPARATION ----------------------
# Simulations
# Sim A: 3 cont, no interaction (Poisson)
N3 <- 600
X3 <- cbind(rnorm(N3,0.3,1.0), rnorm(N3,-0.4,1.1), rnorm(N3,0.0,1.2))
XZ3 <- apply(X3, 2, zscore)
alpha_A <- 1.0; beta_A <- c(0.9,-0.8,0.0); gamma_A <- c(1.0,1.5,1.2)
eta_A <- alpha_A + beta_A[1]*spow(XZ3[,1], gamma_A[1]) + beta_A[2]*spow(XZ3[,2], gamma_A[2])
y_A <- rpois(N3, exp(eta_A))

# Sim B: 3 cont, x1:x2 interaction
kappa_B <- 0.6; psi_i1_B <- 1.0; psi_i2_B <- 0.9
eta_B <- alpha_A +
  beta_A[1]*spow(XZ3[,1], gamma_A[1]) +
  beta_A[2]*spow(XZ3[,2], gamma_A[2]) +
  kappa_B*spow(XZ3[,1], psi_i1_B)*spow(XZ3[,2], psi_i2_B)
y_B <- rpois(N3, exp(eta_B))

# Sim C: 5 cont, mains with powers 1/3 & 3, plus 3 interactions among first 3
N5 <- 700
X5 <- cbind(rnorm(N5), rnorm(N5), rnorm(N5), rnorm(N5), rnorm(N5))
XZ5 <- apply(X5, 2, zscore)
alpha_C <- 0.8
beta_C  <- c(0.8, -0.7, 0.0, 0.5, 0.3)
gamma_C <- c(1.0, 1.5, 1.2, 1/3, 3.0)
kappa_C <- c(0.5, 0.3, -0.4) # pairs (1,2),(1,3),(2,3)
eta_C <- alpha_C +
  beta_C[1]*spow(XZ5[,1], gamma_C[1]) +
  beta_C[2]*spow(XZ5[,2], gamma_C[2]) +
  beta_C[3]*spow(XZ5[,3], gamma_C[3]) +
  beta_C[4]*spow(XZ5[,4], gamma_C[4]) +
  beta_C[5]*spow(XZ5[,5], gamma_C[5]) +
  kappa_C[1]*spow(XZ5[,1],1)*spow(XZ5[,2],1) +
  kappa_C[2]*spow(XZ5[,1],1)*spow(XZ5[,3],1) +
  kappa_C[3]*spow(XZ5[,2],1)*spow(XZ5[,3],1)
y_C <- rpois(N5, exp(eta_C))

# Sim D: factor F (+/-1) + 3 cont; interactions x1:x2 and F:x3 (Poisson)
ND <- 650
x1 <- rnorm(ND, 0.2, 1.0); x2 <- rnorm(ND, -0.5, 1.0); x3 <- rnorm(ND, 0.0, 1.0)
XZ_D <- cbind(zscore(x1), zscore(x2), zscore(x3))
F <- factor(sample(c("A","B"), ND, replace=TRUE)); F_sign <- ifelse(F=="A", 1, -1)
XZ_D <- cbind(XZ_D, zscore(F_sign))  # treat as fourth predictor
alpha_D <- 0.7; beta_D <- c(0.7,-0.6,0.4, 0.3); gamma_D <- c(1.0,1.5,0.7,1.0)
kappa12 <- 0.4; phi_F3 <- 0.5; psi_F3 <- 1.2
eta_D <- alpha_D +
  beta_D[1]*spow(XZ_D[,1], gamma_D[1]) +
  beta_D[2]*spow(XZ_D[,2], gamma_D[2]) +
  beta_D[3]*spow(XZ_D[,3], gamma_D[3]) +
  beta_D[4]*spow(XZ_D[,4], gamma_D[4]) +
  kappa12  * spow(XZ_D[,1],1) * spow(XZ_D[,2],1) +
  phi_F3   * spow(XZ_D[,3], psi_F3) * spow(XZ_D[,4], 1)
y_D <- rpois(ND, exp(eta_D))

# Built-ins
data(ToothGrowth, package="datasets")
TG <- ToothGrowth
TG$dose_z <- zscore(TG$dose); TG$supp_s <- ifelse(TG$supp==levels(TG$supp)[1], 1, -1)
XZ_TG <- cbind(TG$dose_z, zscore(TG$supp_s)) # P=2; interaction is dose×supp
y_TG  <- TG$len

data(warpbreaks, package="datasets")
WB <- warpbreaks
w_s   <- ifelse(WB$wool==levels(WB$wool)[1], 1, -1)
Tmat  <- contr.sum(nlevels(WB$tension))
TX    <- Tmat[as.integer(WB$tension), , drop=FALSE] # two columns
XZ_WB <- cbind(zscore(w_s), scale(TX, center=TRUE, scale=TRUE)) # P=3
y_WB  <- WB$breaks

data(PlantGrowth, package="datasets")
PG <- PlantGrowth
PgMat <- contr.sum(nlevels(PG$group))
PgX   <- PgMat[as.integer(PG$group), , drop=FALSE] # 2 columns
XZ_PG <- scale(PgX, center=TRUE, scale=TRUE)       # P=2
y_PG  <- PG$weight

data(InsectSprays, package="datasets")
ISp <- InsectSprays
SpMat <- contr.sum(nlevels(ISp$spray))
SpX   <- SpMat[as.integer(ISp$spray), , drop=FALSE] # 5 columns
XZ_IS <- scale(SpX, center=TRUE, scale=TRUE)        # P=5
y_IS  <- ISp$count

# ------------------- RUN & REPORT --------------------------
run_both <- function(name, y, XZ, family) {
  cat("\n== ", name, " ==\n", sep="")
  fit_sw <- run_flex(y, XZ, family=family, prior="switch")
  report_fit(name, y, fit_sw)
  fit_hs <- run_flex(y, XZ, family=family, prior="horseshoe")
  report_fit(name, y, fit_hs)
  invisible(list(switch=fit_sw, horseshoe=fit_hs))
}

cat("=== SIMULATIONS ===\n")
resA <- run_both("Sim A (Poisson, 3 vars, no interaction)", y_A, XZ3, "poisson")
resB <- run_both("Sim B (Poisson, 3 vars, x1:x2)",           y_B, XZ3, "poisson")
resC <- run_both("Sim C (Poisson, 5 vars, powers 1/3 & 3, 3 ints)", y_C, XZ5, "poisson")
resD <- run_both("Sim D (Poisson, x1,x2,x3 + factor)",        y_D, XZ_D, "poisson")

cat("\n=== BUILT-INS ===\n")
resTG <- run_both("ToothGrowth (Gaussian; dose × supp)",  y_TG, XZ_TG, "gaussian")
resWB <- run_both("warpbreaks (NegBin; wool × tension)",  y_WB, XZ_WB, "nbinom")
resPG <- run_both("PlantGrowth (Gaussian; groups)",       y_PG, XZ_PG, "gaussian")
resIS <- run_both("InsectSprays (Poisson; spray)",        y_IS, XZ_IS, "poisson")

# ---------------- Greek mapping (ASCII) --------------------
cat("\n--- Greek mapping (ASCII) ---\n")
cat("alpha = intercept\n",
    "beta[j] = main coefficient for predictor j\n",
    "gamma[j] = learned power for predictor j (main term)\n",
    "kappa[m] = interaction coefficient for pair m = (j1[m], j2[m])\n",
    "psi1[m], psi2[m] = learned powers for the two legs of pair m\n",
    "zeta[m] = binary switch for interaction m (switch model only)\n",
    "tau2, lambda2[m] = horseshoe global/local scales (horseshoe only)\n",
    "prec (Gaussian) = residual precision; delta (NegBin) = size parameter\n", sep="")
