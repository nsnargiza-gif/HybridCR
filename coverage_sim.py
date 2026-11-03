#!/usr/bin/env python3
"""
Frequentist-coverage check of HybridCR's nominal-95% intervals, with no R required.

For the two-cause exponential model under hybrid Type-I censoring tau_n = min(X_(r), T0),
HybridCR reports the Wald interval  a_hat_j +/- 1.96 * se,  se = a_hat_j / sqrt(D_j),
with a_hat_j = D_j / TAR (D_j = cause-j failures, TAR = total time at risk). This script
simulates that interval over many replications and reports the empirical coverage of the
true rate, across censoring levels (beta) and sample sizes (n).

Usage:  python3 coverage_sim.py
Requires: numpy only.  Expected output: coverage near 0.95 throughout (Table 2 of the paper).
"""
import numpy as np

def run(a1=1.0, a2=0.5, R=5000, z=1.96, betas=(0.30, 0.50, 0.75, 1.00), ns=(50, 200), seed=424242):
    rng = np.random.default_rng(seed)
    print(f"Two-cause exponential (a1={a1}, a2={a2}); R={R} reps; 95% Wald interval, se=a_hat/sqrt(D).")
    print("Monte-Carlo s.e. of a coverage near 0.95 is ~%.3f\n" % np.sqrt(.95 * .05 / R))
    print(f"{'n':>4} {'beta':>5} {'cov a1':>8} {'cov a2':>8}")
    for n in ns:
        for beta in betas:
            r = int(np.ceil(beta * n))
            c1 = c2 = 0
            for _ in range(R):
                T1 = rng.exponential(1/a1, n); T2 = rng.exponential(1/a2, n)
                X = np.minimum(T1, T2); cause = np.where(T1 < T2, 1, 2)
                tau = np.sort(X)[r-1]                   # T0 = inf  =>  tau = X_(r)
                Z = np.minimum(X, tau); obs = X <= tau
                D1 = int(np.sum(obs & (cause == 1))); D2 = int(np.sum(obs & (cause == 2)))
                TAR = float(np.sum(Z))
                if D1 > 0:
                    ah = D1/TAR; se = ah/np.sqrt(D1)
                    c1 += (ah - z*se <= a1 <= ah + z*se)
                if D2 > 0:
                    ah = D2/TAR; se = ah/np.sqrt(D2)
                    c2 += (ah - z*se <= a2 <= ah + z*se)
            print(f"{n:>4} {beta:>5.2f} {c1/R:>8.3f} {c2/R:>8.3f}")

if __name__ == "__main__":
    run()
