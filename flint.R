#!/usr/bin/env Rscript
# ===========================================================
# flint_batch_psis_pairs.R  (PAIRWISE ONLY, NO K-FOLD HERE)
# - Flexible powers on mains & pairwise interactions
# - Datasets: Sims, Built-ins, NY Tick, GBIF×WorldClim
# - Priors:   horseshoe, switch
# - Samplers: RW, slice, ESS, NUTS (if nimbleHMC present)
# - Powers bounded: Uniform(-3, 3)
# - Computes pointwise log-lik inside model → PSIS-LOO
# - Saves frequent checkpoints to flintplus.Rdata
# ===========================================================
suppressPackageStartupMessages({
  library(nimble)
  library(nimbleHMC)
  library(coda)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(tibble)
  library(loo)
})

set.seed(123)

# ---------------------- SAVE POINTS ------------------------
do_save <- TRUE
save_point <- function(tag = NULL) {
  if (!do_save) return(invisible(NULL))
  try(save.image("flintplus.Rdata", version = 2), silent = TRUE)
  if (!is.null(tag)) message(sprintf(" [save_point: %s -> flintplus.Rdata]", tag))
}

# ---------------------- CONFIG -----------------------------
cfg <- list(
  priors_to_run      = c("horseshoe","switch"),
  sampler_strategies = c("rw","slice","ess","nuts"),

  # MCMC knobs (bump later for long runs)
  niter = 5e5, nburn = 3e5, thin = 5, nchains = 3,

  # GBIF×WorldClim aggregation
  gbif_csv   = "data/gbif_occ/bradypus_variegatus_gbif.csv",
  wc_tif     = "data/worldclim_bio/worldclim_bio_10min.tif",
  gbif_zero_multiplier = 3,
  gbif_topK_predictors = 6,

  ny_tick_dir = "data/ny_tick",
  results_dir = "results",
  collin_dir  = "results/_collinearity_info"
)

dir.create(cfg$results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$collin_dir,  showWarnings = FALSE, recursive = TRUE)

# ---------------------- HELPERS ----------------------------
zscore <- function(x) as.numeric(scale(x))
spow   <- function(x, p, eps=1e-6) sign(x) * (abs(x)+eps)^p
tagify <- function(s) gsub("^_+|_+$","", gsub("[^A-Za-z0-9]+","_", s))

ppc_gaussian <- function(S, y, prefix="y_rep[") {
  yc <- S[, startsWith(colnames(S), prefix), drop=FALSE]
  yhat <- colMeans(yc)
  c(
    obs_mean = mean(y), rep_mean = mean(yhat),
    obs_sd = sd(y),     rep_sd  = sqrt(mean(apply(yc,2,var))),
    RMSE = sqrt(mean((y - yhat)^2)),
    MAE  = mean(abs(y - yhat))
  )
}
ppc_count <- function(S, y, prefix="y_rep[") {
  yc <- S[, startsWith(colnames(S), prefix), drop=FALSE]
  yhat <- colMeans(yc)
  c(
    obs_mean = mean(y), rep_mean = mean(yhat),
    obs_var  = var(y),  rep_var  = mean(apply(yc,2,var)),
    RMSE = sqrt(mean((y - yhat)^2)),
    MAE  = mean(abs(y - yhat))
  )
}

diag_summary <- function(mcmc_list) {
  gd <- try(gelman.diag(mcmc_list, autoburnin = FALSE), silent = TRUE)
  rhat_max <- if (inherits(gd, "try-error")) NA_real_ else {
    ps <- gd$psrf
    col <- if (!is.null(colnames(ps)) && "Point est." %in% colnames(ps)) "Point est." else 1
    suppressWarnings(max(ps[, col], na.rm = TRUE))
  }
  ess <- try(effectiveSize(mcmc_list), silent = TRUE)
  ess_min <- if (inherits(ess, "try-error")) NA_real_ else suppressWarnings(min(ess, na.rm = TRUE))
  hw <- try({
    per_chain <- lapply(mcmc_list, function(ch) {
      hd <- heidel.diag(ch)
      cols <- colnames(hd)
      stest_ok <- if ("stest" %in% cols) mean(hd[, "stest"] == 1, na.rm = TRUE) else NA_real_
      htest_ok <- if ("htest" %in% cols) mean(hd[, "htest"] == 1, na.rm = TRUE) else NA_real_
      start_md <- if ("start" %in% cols) median(hd[, "start"], na.rm = TRUE) else NA_real_
      pval_md  <- if ("pvalue" %in% cols) median(hd[, "pvalue"], na.rm = TRUE) else NA_real_
      c(stest_frac = stest_ok, htest_frac = htest_ok, start_med = start_md, pvalue_med = pval_md)
    })
    mat <- do.call(rbind, per_chain)
    c(Heidel_stable_frac   = mean(mat[, "stest_frac"],  na.rm = TRUE),
      Heidel_halfwidth_frac= mean(mat[, "htest_frac"],  na.rm = TRUE),
      Heidel_start_median  = median(mat[, "start_med"],  na.rm = TRUE),
      Heidel_pvalue_median = median(mat[, "pvalue_med"], na.rm = TRUE))
  }, silent = TRUE)
  if (inherits(hw, "try-error")) {
    c(Rhat_max = rhat_max, ESS_min = ess_min,
      Heidel_stable_frac = NA_real_, Heidel_halfwidth_frac = NA_real_,
      Heidel_start_median = NA_real_, Heidel_pvalue_median = NA_real_)
  } else c(Rhat_max = rhat_max, ESS_min = ess_min, hw)
}

build_pairs <- function(P) {
  if (P < 2) return(list(p1 = integer(0), p2 = integer(0)))
  cmb <- combn(P, 2); list(p1 = as.integer(cmb[1, ]), p2 = as.integer(cmb[2, ]))
}

# ---------- COLLINEARITY (info only; no dropping) ----------
collin_report <- function(name, X, path_txt) {
  cat(sprintf("\n[Collinearity] %s: P=%d\n", name, ncol(X)))
  con <- file(path_txt, open="wt"); on.exit(close(con), add=TRUE)
  writeLines(sprintf("# Collinearity report: %s", name), con)
  writeLines(sprintf("P = %d", ncol(X)), con)
  if (ncol(X) < 2) { writeLines("Too few predictors.", con); return(invisible(NULL)) }

  Xc <- scale(X, TRUE, TRUE)
  C  <- try(stats::cor(Xc), silent=TRUE)
  if (!inherits(C, "try-error") && all(is.finite(C))) {
    ev <- eigen(C, symmetric=TRUE, only.values=TRUE)$values
    kappa_cond <- sqrt(max(ev, na.rm=TRUE)/min(ev, na.rm=TRUE))
    writeLines(sprintf("Condition number (corr): %.2f", kappa_cond), con)
    writeLines(sprintf("Max |corr| off-diagonal: %.3f", max(abs(C[row(C)!=col(C)]), na.rm=TRUE)), con)
  } else writeLines("Correlation matrix unstable/NA.", con)

  if (ncol(X) <= 30) {
    vifs <- rep(NA_real_, ncol(X))
    for (j in seq_len(ncol(X))) {
      df <- as.data.frame(Xc)
      fit <- try(lm(df[[j]] ~ . , data = df[ , -j, drop=FALSE]), silent=TRUE)
      if (!inherits(fit, "try-error")) {
        R2 <- summary(fit)$r.squared
        vifs[j] <- 1/(1 - max(R2, 0))
      }
    }
    writeLines(sprintf("VIF: median=%.2f | max=%.2f", median(vifs, na.rm=TRUE), max(vifs, na.rm=TRUE)), con)
  } else writeLines("VIF skipped (P>30).", con)
}

# ----------------- SAMPLER ASSIGNMENT ----------------------
assign_samplers <- function(conf, model, P, Q,
                            strategy = c("auto","slice","rw","ess","nuts"),
                            monitors = NULL) {
  strategy <- match.arg(strategy)
  betas   <- paste0("beta[",  1:P, "]")
  gammas  <- paste0("gamma[", 1:P, "]")
  psi1s   <- if (Q>0) paste0("psi1[",  1:Q, "]") else character()
  psi2s   <- if (Q>0) paste0("psi2[",  1:Q, "]") else character()
  kappas  <- if (Q>0) paste0("kappa[", 1:Q, "]") else character()

  rm_add <- function(nodes, type) {
    nodes <- nodes[nodes %in% model$getNodeNames()]
    if (!length(nodes)) return()
    for (nm in nodes) {
      conf$removeSamplers(nm)
      conf$addSampler(target = nm, type = type, control = list())  # explicit target + control
    }
  }

  if (strategy %in% c("auto","slice")) { rm_add(c(betas, gammas, psi1s, psi2s, kappas), "slice"); return(conf) }
  if (strategy == "rw")                { rm_add(c(betas, gammas, psi1s, psi2s, kappas), "RW");    return(conf) }
  if (strategy == "ess") {
    rm_add(betas, "ess"); rm_add(kappas, "ess")
    rm_add(c(gammas, psi1s, psi2s), "slice")
    return(conf)
  }

  # ---- NUTS / HMC ----
  if (strategy == "nuts") {
    if (!requireNamespace("nimbleHMC", quietly = TRUE)) {
      message("nimbleHMC not installed; using slice for 'nuts'.")
      rm_add(c(betas, gammas, psi1s, psi2s, kappas), "slice")
      return(conf)
    }
    # Prefer configureHMC() which builds a NEW config (keeps monitors if we pass them)
    if ("configureHMC" %in% getNamespaceExports("nimbleHMC")) {
      conf_hmc <- nimbleHMC::configureHMC(model, type = "NUTS",
                                          monitors = monitors, print = FALSE)
      return(conf_hmc)
    }
    # Fallback: retrofit an existing conf (requires replace=TRUE)
    if ("addHMC" %in% getNamespaceExports("nimbleHMC")) {
      nimbleHMC::addHMC(conf, type = "NUTS", replace = TRUE)
      return(conf)
    }
    message("nimbleHMC present but no usable HMC config fn; using slice.")
    rm_add(c(betas, gammas, psi1s, psi2s, kappas), "slice")
    return(conf)
  }

  return(conf)
}


# -------- nimbleCode generator (pairs + ll[i]) --------------
code_pairs_Qaware <- function(N, P, Q, FAM = c(1,2,3), prior = c("switch","horseshoe")) {
  prior <- match.arg(prior)
  FAM <- match.arg(as.character(FAM), choices = c("1","2","3"))
  FAM <- as.integer(FAM)

  nimbleCode({
    # mains
    for (i in 1:N) {
      for (j in 1:P) {
        s[i,j] <- 2*step(xz[i,j]) - 1
        a[i,j] <- abs(xz[i,j]) + eps
        main_contrib[i,j] <- beta[j] * s[i,j] * pow(a[i,j], gamma[j])
      }
      if (P > 1) { lin_main[i] <- sum(main_contrib[i,1:P]) } else { lin_main[i] <- main_contrib[i,1] }
    }

    # pairs
    if (Q > 0) {
      for (i in 1:N) {
        for (m in 1:Q) {
          f1p[i,m] <- s[i, j1[m]] * pow(a[i, j1[m]], psi1[m])
          f2p[i,m] <- s[i, j2[m]] * pow(a[i, j2[m]], psi2[m])
          int_pair[i,m] <- MULT2[m] * kappa[m] * f1p[i,m] * f2p[i,m]
        }
        lin_pairs[i] <- sum(int_pair[i,1:Q])
      }
    } else {
      for (i in 1:N) lin_pairs[i] <- 0
    }

    # likelihood + pointwise log-lik
    if (FAM == 1) {  # Poisson
      for (i in 1:N) {
        log(mu[i]) <- alpha + lin_main[i] + lin_pairs[i]
        y[i] ~ dpois(mu[i]); y_rep[i] ~ dpois(mu[i])
        ll[i] <- dpois(y[i], mu[i], log=1)
      }
    }
    if (FAM == 2) {  # Gaussian
      for (i in 1:N) {
        mu_g[i] <- alpha + lin_main[i] + lin_pairs[i]
        y[i] ~ dnorm(mu_g[i], prec); y_rep[i] ~ dnorm(mu_g[i], prec)
        ll[i] <- dnorm(y[i], mean=mu_g[i], prec=prec, log=1)
      }
      prec ~ dgamma(1,1)
    }
    if (FAM == 3) {  # NegBin
      for (i in 1:N) {
        log(mu_nb[i]) <- alpha + lin_main[i] + lin_pairs[i]
        y[i] ~ dnbinom(size=delta, prob = delta / (delta + mu_nb[i]))
        y_rep[i] ~ dnbinom(size=delta, prob = delta / (delta + mu_nb[i]))
        ll[i] <- dnbinom(y[i], size=delta, prob = delta / (delta + mu_nb[i]), log=1)
      }
      delta ~ dgamma(0.5,0.5)
    }

    # priors (use compile-time flags!)
    alpha ~ dnorm(0, 1.0E-4)
    for (j in 1:P) { beta[j] ~ dnorm(0, prec_beta); gamma[j] ~ dunif(-3, 3) }
    prec_beta ~ dgamma(0.5, 0.5)

    if (USE_SWITCH == 1) {
      if (Q > 0) {
        for (m in 1:Q) {
          psi1[m] ~ dunif(-3, 3); psi2[m] ~ dunif(-3, 3)
          zeta2[m] ~ dbern(pi2); MULT2[m] <- zeta2[m]
          kappa[m] ~ dnorm(0, prec_kappa2)
        }
        pi2 ~ dbeta(1,9); prec_kappa2 ~ dgamma(0.5,0.5)
        kappa_mag2 <- sum(pow(kappa[1:Q],2))
      } else kappa_mag2 <- 0
    }

    if (USE_HS == 1) {
      if (Q > 0) {
        tau2_k ~ dinvgamma(0.5, 1 / xi_k);  xi_k ~ dinvgamma(0.5, 1)
        for (m in 1:Q) {
          lambda2_k[m] ~ dinvgamma(0.5, 1 / nu_k[m]);  nu_k[m] ~ dinvgamma(0.5, 1)
          prec_kappa_m[m] <- 1 / (tau2_k * lambda2_k[m])
          kappa[m] ~ dnorm(0, prec_kappa_m[m]); MULT2[m] <- 1
          psi1[m] ~ dunif(-3, 3); psi2[m] ~ dunif(-3, 3)
        }
        kappa_mag2 <- sum(pow(kappa[1:Q],2))
      } else { tau2_k <- 1; kappa_mag2 <- 0 }
    }
  })
}


# --------------- LOG-LIK EXTRACTOR FOR LOO -----------------
extract_loglik <- function(mcmc_list) {
  loglik_chain <- list(); chain_id <- integer(0)
  for (ch in seq_along(mcmc_list)) {
    M <- as.matrix(mcmc_list[[ch]])
    keep <- grepl("^ll\\[[0-9]+\\]$", colnames(M))
    if (!any(keep)) stop("No ll[i] columns found in MCMC output.")
    Mll <- M[, keep, drop=FALSE]
    idx <- as.integer(gsub("^ll\\[|\\]$", "", colnames(Mll)))
    Mll <- Mll[, order(idx), drop=FALSE]
    loglik_chain[[ch]] <- Mll
    chain_id <- c(chain_id, rep(ch, nrow(Mll)))
  }
  list(log_lik = do.call(rbind, loglik_chain), chain_id = chain_id)
}

compute_psis_loo <- function(mcmc_list) {
  ex <- extract_loglik(mcmc_list)
  r_eff <- try(loo::relative_eff(exp(ex$log_lik), chain_id = ex$chain_id), silent = TRUE)
  if (inherits(r_eff, "try-error")) r_eff <- NULL
  loo_out <- loo::loo(ex$log_lik, r_eff = r_eff, cores = 1)
  kvec <- try(as.numeric(loo_out$diagnostics$pareto_k), silent = TRUE)
  if (inherits(kvec, "try-error") || is.null(kvec)) kvec <- rep(NA_real_, ncol(ex$log_lik))
  list(
    loo = loo_out,
    k_max = suppressWarnings(max(kvec, na.rm = TRUE)),
    k_frac_gt_05 = mean(kvec > 0.5, na.rm = TRUE),
    k_frac_gt_07 = mean(kvec > 0.7, na.rm = TRUE),
    k_frac_gt_10 = mean(kvec > 1.0, na.rm = TRUE),
    k_values = kvec
  )
}

# ------------------- single fit ----------------------------
run_flex <- function(y, XZ, family = c("poisson","gaussian","nbinom"),
                     prior = c("switch","horseshoe"),
                     sampler_strategy = c("auto","slice","rw","ess","nuts"),
                     niter=cfg$niter, nburn=cfg$nburn, thin=cfg$thin, nchains=cfg$nchains) {

  family <- match.arg(family); prior <- match.arg(prior)
  sampler_strategy <- match.arg(sampler_strategy)

  N <- nrow(XZ); P <- ncol(XZ)
  pr <- build_pairs(P); Q <- length(pr$p1)

  code <- code_pairs_Qaware(N,P,Q, FAM = switch(family, poisson=1L, gaussian=2L, nbinom=3L), prior=prior)

  # ALWAYS pass j1/j2 as vectors (length Q), even if Q==1; pass integer(0) if Q==0
consts <- list(
  N = N, P = P, Q = Q, eps = 1e-6,
  FAM = switch(family, poisson = 1L, gaussian = 2L, nbinom = 3L),
  USE_SWITCH = as.integer(prior == "switch"),
  USE_HS     = as.integer(prior == "horseshoe"),
  j1 = if (Q > 0) as.integer(pr$p1) else integer(0),
  j2 = if (Q > 0) as.integer(pr$p2) else integer(0)
)
  data <- list(y = as.numeric(y), xz = as.matrix(XZ))

  init_fun <- function() {
    ini <- list(alpha = 0, beta = rnorm(P,0,0.5), gamma = runif(P,-0.2,0.2), prec_beta = 1)
    if (Q > 0) { ini$psi1 <- runif(Q,-0.2,0.2); ini$psi2 <- runif(Q,-0.2,0.2); ini$kappa <- rnorm(Q,0,0.1) }
    if (prior == "switch" && Q>0) { ini$pi2 <- 0.2; ini$zeta2 <- rbinom(Q,1,0.1); ini$prec_kappa2 <- 1 }
    if (prior == "horseshoe" && Q>0) { ini$tau2_k <- 1; ini$xi_k <- 1; ini$lambda2_k <- rep(1,Q); ini$nu_k <- rep(1,Q) }
    if (family == "gaussian") ini$prec <- 1/var(y)
    if (family == "nbinom")   ini$delta <- 1
    ini
  }

  cat("Defining model\n")
  model <- nimbleModel(code, constants=consts, data=data, inits=init_fun(), buildDerivs = TRUE)

  monitors <- c("alpha",
                if (P>0) c(paste0("beta[",1:P,"]"), paste0("gamma[",1:P,"]")) else NULL,
                if (Q>0) c(paste0("psi1[",1:Q,"]"), paste0("psi2[",1:Q,"]"), paste0("kappa[",1:Q,"]"), "kappa_mag2") else NULL,
                paste0("y_rep[",1:N,"]"),
                paste0("ll[",1:N,"]"),
                if (family=="gaussian") "prec" else NULL,
                if (family=="nbinom")  "delta" else NULL,
                if (prior=="horseshoe" && Q>0) "tau2_k" else NULL,
                if (prior=="switch"    && Q>0) "pi2"    else NULL)
  
  conf <- configureMCMC(model, monitors = monitors)
  conf <- assign_samplers(conf, model, P, Q, strategy = sampler_strategy, monitors = monitors)

conf$printSamplers(byType = TRUE)

  mcmc   <- buildMCMC(conf)
  cmodel <- compileNimble(model)
  cmcmc  <- compileNimble(mcmc, project=model)

  tm <- system.time({
    samps <- runMCMC(cmcmc, niter=niter, nburnin=nburn, thin=thin, nchains=nchains, samplesAsCodaMCMC=TRUE)
  })

  ps <- compute_psis_loo(samps)
  S <- as.matrix(do.call(rbind, samps))
  list(samples=samps, Smat=S, time=unname(tm["elapsed"]),
       family=family, P=P, Q=Q, prior=prior, strategy=sampler_strategy,
       loo=ps$loo, k_max=ps$k_max, k_frac_gt_05=ps$k_frac_gt_05, k_frac_gt_07=ps$k_frac_gt_07, k_frac_gt_10=ps$k_frac_gt_10,
       k_values=ps$k_values)
}

# --------------- report row + print ------------------------
tidy_report <- function(name, y, fit) {
  fam <- fit$family; S <- fit$Smat
  ppc <- if (fam=="gaussian") ppc_gaussian(S,y) else ppc_count(S,y)
  di <- diag_summary(fit$samples)
  lo <- fit$loo
  row <- tibble(
    dataset = name, prior = fit$prior, sampler = fit$strategy,
    runtime_sec = round(fit$time,2),
    Rhat_max = as.numeric(di["Rhat_max"]),
    ESS_min  = as.numeric(di["ESS_min"]),
    Heidel_stable_frac    = as.numeric(di["Heidel_stable_frac"]),
    Heidel_halfwidth_frac = as.numeric(di["Heidel_halfwidth_frac"]),
    RMSE = as.numeric(ppc["RMSE"]), MAE  = as.numeric(ppc["MAE"]),
    LOO_elpd   = as.numeric(lo$estimates["elpd_loo","Estimate"]),
    LOO_elpd_se= as.numeric(lo$estimates["elpd_loo","SE"]),
    LOO_p_loo  = as.numeric(lo$estimates["p_loo","Estimate"]),
    LOOIC      = as.numeric(lo$estimates["looic","Estimate"]),
    Pareto_k_max = fit$k_max,
    Pareto_k_gt0.5 = fit$k_frac_gt_05,
    Pareto_k_gt0.7 = fit$k_frac_gt_07,
    Pareto_k_gt1.0 = fit$k_frac_gt_10,
    fam = fit$family, P = fit$P, Q_pairs = fit$Q
  )
  print(row)
  row
}

# --------------- full suite per dataset --------------------
run_suite <- function(name, y, XZ, family) {
  # collinearity info (on full XZ)
  collin_path <- file.path(cfg$collin_dir, paste0(tagify(name), ".txt"))
  try(collin_report(name, XZ, collin_path), silent=TRUE)

  rows <- list()
  for (pr in cfg$priors_to_run) {
    for (st in cfg$sampler_strategies) {
      cat(sprintf("\n== %s | prior=%s | sampler=%s ==\n", name, pr, toupper(st)))
      fit <- run_flex(y, XZ, family=family, prior=pr, sampler_strategy=st)
      # save Pareto-k vector
      pk_path <- file.path(cfg$results_dir, paste0("pareto_k_", tagify(name), "_", pr, "_", toupper(st), ".csv"))
      readr::write_csv(tibble(obs = seq_along(fit$k_values), pareto_k = fit$k_values), pk_path)

      rows[[length(rows)+1]] <- tidy_report(name, y, fit)
      if (is.finite(fit$k_max) && fit$k_max > 0.7) {
        message(sprintf("  [PSIS warning] max Pareto-k=%.3f (>0.7). Consider k-fold later.", fit$k_max))
      }
      save_point(sprintf("%s_%s_%s", tagify(name), pr, st))
    }
  }
  bind_rows(rows)
}

# ================== DATA PREP ==============================
mk_simulations <- function() {
  # Sim A
  N3 <- 600; X3 <- cbind(rnorm(N3,0.3,1.0), rnorm(N3,-0.4,1.1), rnorm(N3,0.0,1.2))
  XZ3 <- apply(X3, 2, zscore)
  alpha_A <- 1.0; beta_A <- c(0.9,-0.8,0.0); gamma_A <- c(1.0,1.5,1.2)
  eta_A <- alpha_A + beta_A[1]*spow(XZ3[,1], gamma_A[1]) + beta_A[2]*spow(XZ3[,2], gamma_A[2])
  y_A <- rpois(N3, exp(eta_A))

  # Sim B
  kappa_B <- 0.6; psi_i1_B <- 1.0; psi_i2_B <- 0.9
  eta_B <- alpha_A + beta_A[1]*spow(XZ3[,1], gamma_A[1]) + beta_A[2]*spow(XZ3[,2], gamma_A[2]) +
           kappa_B*spow(XZ3[,1], psi_i1_B)*spow(XZ3[,2], psi_i2_B)
  y_B <- rpois(N3, exp(eta_B))

  # Sim C
  N5 <- 700; X5 <- cbind(rnorm(N5), rnorm(N5), rnorm(N5), rnorm(N5), rnorm(N5))
  XZ5 <- apply(X5, 2, zscore)
  alpha_C <- 0.6
  beta_C  <- c(0.5, -0.4, 0.0, 0.3, 0.15)
  gamma_C <- c(1.0, 1.5, 1.2, 1/3, 3.0)
  kappa_C <- c(0.3, 0.2, -0.2)
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

  list(
    list(name="Sim A (Pois, 3, no int)", y=y_A, XZ=unname(XZ3), fam="poisson"),
    list(name="Sim B (Pois, 3, x1:x2)", y=y_B, XZ=unname(XZ3), fam="poisson"),
    list(name="Sim C (Pois, 5, 3 pairs)", y=y_C, XZ=unname(XZ5), fam="poisson")
  )
}

mk_builtins <- function() {
  data(ToothGrowth, package="datasets")
  TG <- ToothGrowth
  TG$dose_z <- zscore(TG$dose); TG$supp_s <- ifelse(TG$supp==levels(TG$supp)[1], 1, -1)
  XZ_TG <- cbind(TG$dose_z, zscore(TG$supp_s)); y_TG <- TG$len

  data(warpbreaks, package="datasets")
  WB <- warpbreaks
  w_s <- ifelse(WB$wool==levels(WB$wool)[1], 1, -1)
  Tmat <- contr.sum(nlevels(WB$tension))
  TX   <- Tmat[as.integer(WB$tension), , drop=FALSE]
  XZ_WB <- cbind(zscore(w_s), scale(TX, TRUE, TRUE)); y_WB <- WB$breaks

  data(PlantGrowth, package="datasets")
  PG <- PlantGrowth
  PgMat <- contr.sum(nlevels(PG$group))
  PgX   <- PgMat[as.integer(PG$group), , drop=FALSE]
  XZ_PG <- scale(PgX, TRUE, TRUE); y_PG <- PG$weight

  data(InsectSprays, package="datasets")
  ISp <- InsectSprays
  SpMat <- contr.sum(nlevels(ISp$spray))
  SpX   <- SpMat[as.integer(ISp$spray), , drop=FALSE]
  XZ_IS <- scale(SpX, TRUE, TRUE); y_IS <- ISp$count

  list(
    list(name="ToothGrowth (Gauss; dose×supp)", y=y_TG, XZ=unname(XZ_TG), fam="gaussian"),
    list(name="warpbreaks (NegBin; wool×tension)", y=y_WB, XZ=unname(XZ_WB), fam="nbinom"),
    list(name="PlantGrowth (Gauss; groups)", y=y_PG, XZ=unname(XZ_PG), fam="gaussian"),
    list(name="InsectSprays (Pois; spray)", y=y_IS, XZ=unname(XZ_IS), fam="poisson")
  )
}

mk_ny_tick <- function(dir=cfg$ny_tick_dir) {
  pick_response <- function(df) {
    nm <- names(df); num <- nm[vapply(df, is.numeric, logical(1))]
    pref <- c("density","tickdensity","ticks","count","risk","encounter","infection","rate")
    candidates <- num[order(match(tolower(num), pref), na.last=NA)]
    drop <- nm[tolower(nm) %in% c("year","collectionyear","period","season")]
    candidates <- setdiff(candidates, drop)
    if (length(candidates)) candidates[1] else num[1]
  }
  prep_one <- function(path, label) {
    df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE, guess_max = 100000))
    for (j in names(df)) if (is.character(df[[j]]) && dplyr::n_distinct(df[[j]])<=12) df[[j]] <- factor(df[[j]])
    ycol <- pick_response(df); y <- df[[ycol]]

    fam <- if (all(is.finite(y)) && all(y >= 0) && all(abs(y - round(y)) < 1e-8)) "poisson" else "gaussian"

    keep <- setdiff(names(df), ycol)
    Xnum <- df[keep[vapply(df[keep], is.numeric, logical(1))]]
    Xfac <- df[keep[vapply(df[keep], is.factor,  logical(1))]]
    XZ <- NULL

    if (ncol(as.data.frame(Xnum))>0) XZ <- as.matrix(apply(Xnum,2,zscore))
    if (ncol(as.data.frame(Xfac))>0) {
      mm <- NULL
      for (cn in names(Xfac)) {
        M <- contr.sum(nlevels(Xfac[[cn]]))
        X <- M[as.integer(Xfac[[cn]]), , drop=FALSE]
        mm <- if (is.null(mm)) X else cbind(mm,X)
      }
      mm <- scale(mm, TRUE, TRUE)
      XZ <- if (is.null(XZ)) mm else cbind(XZ, mm)
    }
    XZ <- as.matrix(XZ)
    ok <- is.finite(y) & apply(XZ,1,function(r) all(is.finite(r)))
    y <- as.numeric(y[ok]); XZ <- XZ[ok,,drop=FALSE]
    if (ncol(XZ) < 2) return(NULL)
    list(name=paste0("NYTick ", label, " (", fam, ")"), y=y, XZ=unname(XZ), fam=fam)
  }

  items <- list()
  for (k in c("adults","nymphs","risk_adult","risk_nymph")) {
    pth <- file.path(dir, paste0(k, ".csv"))
    if (file.exists(pth)) {
      it <- try(prep_one(pth, k), silent=TRUE)
      if (!inherits(it, "try-error") && !is.null(it)) items[[length(items)+1]] <- it
    }
  }
  items
}

mk_gbif_worldclim <- function(gbif_csv=cfg$gbif_csv, wc_tif=cfg$wc_tif, topK=cfg$gbif_topK_predictors) {
  if (!file.exists(gbif_csv) || !file.exists(wc_tif)) return(list())

  occ <- suppressMessages(readr::read_csv(gbif_csv, show_col_types = FALSE, guess_max = 100000)) %>%
    filter(is.finite(decimalLatitude), is.finite(decimalLongitude))
  r <- terra::rast(wc_tif)

  pts <- terra::vect(occ[,c("decimalLongitude","decimalLatitude")], geom=c("decimalLongitude","decimalLatitude"), crs="EPSG:4326")
  cells <- terra::cellFromXY(r[[1]], terra::geom(pts))
  tab <- as.data.frame(table(cells), stringsAsFactors = FALSE)
  tab$cells <- as.integer(tab$cells); tab$y <- as.integer(tab$Freq); tab$Freq <- NULL

  xy <- terra::xyFromCell(r[[1]], tab$cells)
  cov <- as.data.frame(terra::extract(r, terra::vect(xy, geom=c("x","y"), crs="EPSG:4326")))
  colnames(cov) <- paste0("bio", sprintf("%02d", 1:19))
  pres <- cbind(tab["y"], cov) %>% drop_na()

  mask <- !is.na(r[[1]])
  n0 <- cfg$gbif_zero_multiplier * nrow(pres)
  set.seed(42)
  zcells <- sample(which(values(mask)[,1] == 1), size = min(n0, sum(values(mask)[,1]==1)))
  zcov <- as.data.frame(terra::extract(r, zcells)) %>% drop_na()
  zeros <- tibble(y = 0) %>% bind_cols(zcov)

  dat <- bind_rows(pres, zeros) %>% drop_na()
  X <- as.matrix(dat[ , -1, drop=FALSE]); y <- dat$y

  sds <- apply(X,2,sd)
  keep_by_var <- names(sort(sds, decreasing = TRUE))[1:min(topK*2, ncol(X))]
  X <- X[, keep_by_var, drop=FALSE]
  sc <- abs(sapply(seq_len(ncol(X)), function(j) cor(log1p(y), X[,j], use="pair")))
  ord <- order(sc, decreasing = TRUE)[1:min(topK, length(sc))]
  X <- X[, ord, drop=FALSE]

  XZ <- scale(X, TRUE, TRUE)
  list(list(name = sprintf("GBIF×WorldClim (Pois; K=%d)", ncol(XZ)), y = y,
            XZ = unname(XZ), fam = "poisson"))
}

# ================= MAIN (all datasets) =====================
datasets <- c(mk_simulations(), mk_builtins(), mk_ny_tick(), mk_gbif_worldclim())

all_rows <- list()
for (d in datasets) {
  if (is.null(d)) next
  message("\n=== DATASET: ", d$name, " | P=", ncol(d$XZ), " ===")
  # collinearity info
  collin_path <- file.path(cfg$collin_dir, paste0(tagify(d$name), ".txt"))
  try(collin_report(d$name, d$XZ, collin_path), silent=TRUE)

  res_list <- list()
  for (pr in cfg$priors_to_run) {
    for (st in cfg$sampler_strategies) {
      cat(sprintf("\n== %s | prior=%s | sampler=%s ==\n", d$name, pr, toupper(st)))
      fit <- run_flex(d$y, d$XZ, family=d$fam, prior=pr, sampler_strategy=st)
      # save Pareto-k vector
      pk_path <- file.path(cfg$results_dir, paste0("pareto_k_", tagify(d$name), "_", pr, "_", toupper(st), ".csv"))
      readr::write_csv(tibble(obs = seq_along(fit$k_values), pareto_k = fit$k_values), pk_path)

      res_list[[length(res_list)+1]] <- tidy_report(d$name, d$y, fit)
      if (is.finite(fit$k_max) && fit$k_max > 0.7) {
        message(sprintf("  [PSIS warning] max Pareto-k=%.3f (>0.7). Consider k-fold later.", fit$k_max))
      }
      save_point(sprintf("%s_%s_%s", tagify(d$name), pr, st))
    }
  }
  res <- bind_rows(res_list)
  out_path <- file.path(cfg$results_dir, paste0(tagify(d$name), ".csv"))
  readr::write_csv(res, out_path)
  all_rows[[length(all_rows)+1]] <- res
}
if (length(all_rows)) {
  summary_path <- file.path(cfg$results_dir, "_summary_all.csv")
  readr::write_csv(bind_rows(all_rows), summary_path)
}
save_point("after_main_runs")

cat("\nDone. Results in '", cfg$results_dir, "'. Frequent checkpoints in flintplus.Rdata\n", sep="")
