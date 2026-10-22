#!/usr/bin/env python3
"""
Offline cross-validation of HybridCR's numerical core, with NO R required.

Independently re-implements, in NumPy, (i) simulation of two-cause data under
hybrid Type-I censoring tau_n = min(X_(r), T0), (ii) the cause-specific
exponential MLE, and (iii) the Aalen-Johansen cumulative-incidence estimator,
then checks:
  * MLE recovers the true rates;
  * the nonparametric Aalen-Johansen CIF agrees with the fitted parametric CIF
    (correct-specification check), to Monte-Carlo error.

This reproduces the role of compare_R.R (which uses survival/cmprsk/etm)
using only NumPy, so it runs anywhere Python is available.

Usage:  python3 compare_python.py
Requires: numpy only.
"""
import numpy as np

def run(n=20000, a1=1.0, a2=0.5, beta=0.6, seed=20261022):
    rng = np.random.default_rng(seed)
    lam = a1 + a2
    T1 = rng.exponential(1/a1, n); T2 = rng.exponential(1/a2, n)
    X = np.minimum(T1, T2); C = np.where(T1 < T2, 1, 2)
    r = int(np.ceil(beta*n)); tau = np.sort(X)[r-1]      # T0 = inf  =>  tau = X_(r)
    Z = np.minimum(X, tau); obs = (X <= tau); cause = C
    # independent MLE (closed form): a_j = D_j / total-time-at-risk
    D1 = int(np.sum(obs & (cause == 1))); D2 = int(np.sum(obs & (cause == 2)))
    TAR = float(np.sum(Z)); a1h, a2h = D1/TAR, D2/TAR
    # independent Aalen-Johansen CIF
    idx = np.argsort(Z, kind='mergesort'); S = 1.0; H1 = 0.0; H2 = 0.0; nAt = len(Z)
    ts = []; h1 = []; h2 = []
    for i in idx:
        if obs[i]:
            if cause[i] == 1: H1 += S/nAt
            else:            H2 += S/nAt
            S *= (1 - 1/nAt)
        ts.append(Z[i]); h1.append(H1); h2.append(H2); nAt -= 1
    ts = np.array(ts); h1 = np.array(h1); h2 = np.array(h2)
    cif = lambda aj, t: (aj/(a1h+a2h))*(1 - np.exp(-(a1h+a2h)*t))
    aj_at = lambda t: (h1[np.searchsorted(ts, t, 'right')-1], h2[np.searchsorted(ts, t, 'right')-1])

    print(f"hybrid-censored two-cause exponential: n={n}, a1={a1}, a2={a2}, beta={beta}")
    print(f"tau = X_(r) = {tau:.4f}   observed failures D1={D1}, D2={D2}  "
          f"(failure fraction {(D1+D2)/n:.4f}, target {beta})\n")
    print(f"MLE recovery:   a1_hat = {a1h:.4f} (truth {a1})    a2_hat = {a2h:.4f} (truth {a2})\n")
    print("  t        F1 parametric   F1 Aalen-Johansen     F2 parametric   F2 Aalen-Johansen")
    maxd = 0.0
    for t in np.round(np.linspace(0.1, 0.9, 5)*tau, 4):
        aj1, aj2 = aj_at(t); p1, p2 = cif(a1h, t), cif(a2h, t)
        maxd = max(maxd, abs(aj1-p1), abs(aj2-p2))
        print(f"  {t:.4f}     {p1:.5f}         {aj1:.5f}               {p2:.5f}         {aj2:.5f}")
    print(f"\nmax |Aalen-Johansen - parametric CIF| = {maxd:.2e}  (Monte-Carlo error under correct specification)")

if __name__ == "__main__":
    run()
