## HybridCR — cross-validation against reference R packages.
## Run:  Rscript compare/compare_R.R
## Requires: survival, cmprsk, etm   install.packages(c("survival","cmprsk","etm"))
## It simulates a two-cause exponential sample under the same hybrid Type-I scheme
## as HybridCR and compares (a) cause-specific MLEs and (b) the cumulative incidence
## functions against the Aalen–Johansen estimator (survival::survfit / etm).

suppressMessages({library(survival); library(cmprsk); library(etm)})
set.seed(42)

sim_hybrid <- function(n, a1, a2, beta){
  T1 <- rexp(n, a1); T2 <- rexp(n, a2)
  X  <- pmin(T1, T2); cause <- ifelse(T1 < T2, 1L, 2L)
  lam <- a1 + a2; T0 <- -log(1 - beta)/lam          # calendar horizon for fraction beta
  Z   <- pmin(X, T0); status <- ifelse(X <= T0, cause, 0L)  # 0 = censored
  data.frame(time = Z, status = status)
}

a1 <- 1.0; a2 <- 0.5; beta <- 0.75; n <- 5000
d  <- sim_hybrid(n, a1, a2, beta)

## (a) cause-specific exponential MLE  a_j = D_j / sum(Z)   (same estimator as HybridCR)
TAR <- sum(d$time); D1 <- sum(d$status==1); D2 <- sum(d$status==2)
cat(sprintf("Cause-specific MLE:  a1=%.4f (true %.2f)   a2=%.4f (true %.2f)\n",
            D1/TAR, a1, D2/TAR, a2))

## (b) Aalen–Johansen CIF via survival::survfit (multi-state), compared at a grid
fs  <- factor(d$status, levels=c(0,1,2), labels=c("censored","cause1","cause2"))
fit <- survfit(Surv(time, fs) ~ 1)
tg  <- quantile(d$time[d$status>0], c(.25,.5,.75))
cat("\nAalen-Johansen CIF (survival::survfit) at t =",
    paste(round(tg,2), collapse=", "), ":\n")
print(summary(fit, times = tg)$pstate)

## Parametric CIF for the fitted exponential model:  H_j(t)= (a_j/lam)(1-exp(-lam t))
lam_hat <- D1/TAR + D2/TAR
cif_par <- function(t, aj) (aj/lam_hat)*(1-exp(-lam_hat*t))
cat("\nParametric CIF (HybridCR model) at the same t:\n")
cat("  cause1:", round(cif_par(tg, D1/TAR),4), "\n")
cat("  cause2:", round(cif_par(tg, D2/TAR),4), "\n")

## cmprsk::cuminc as a second nonparametric reference
ci <- cuminc(d$time, d$status, cencode = 0)
cat("\n(cmprsk::cuminc computed; inspect ci for the two cause-specific CIFs.)\n")
cat("\nExpected: parametric and Aalen-Johansen CIFs agree to Monte-Carlo error\n",
    "under correct (exponential) specification.\n")
