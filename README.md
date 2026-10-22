# HybridCR

**A self-contained, browser-based tool for maximum-likelihood inference and study design in two-cause competing-risks models under hybrid Type-I censoring.**

HybridCR is a single HTML file (~110 kB) with no third-party libraries, no build step and no server. Open it in any modern browser — it runs offline, all computation is client-side, and no data leave your device. It is aimed at teaching, at pre-study design, and as an interactive companion to the methodological results in the references below.

Version **v1.0.0** · License **MIT** · Languages HTML5 / CSS3 / vanilla JavaScript (ES6), inline SVG.

---

## Quick start

- **Run:** double-click `HybridCR.html`, or open it in a browser, or serve the folder as a static site. Nothing to install.
- **Verified in:** Chromium 141 (the Chrome/Edge engine), on desktop and a mobile viewport, with no console errors. It uses only standard web APIs and is expected to run in current Firefox and Safari.

## Modules

1. **Simulation & estimation** — draw a two-cause sample under the hybrid scheme τₙ = min(X₍ᵣ₎, T₀), fit cause-specific exponential/Weibull intensities by Newton–Raphson, with cumulative incidence (parametric vs Aalen–Johansen), Cox–Snell diagnostics, and bootstrap CIs.
2. **Operational information & efficiency bound** — the operational Fisher information 𝓘ₒₚ vs the classical 𝓘_F; equals β·𝓘_F only for the exponential family, computed by quadrature for Weibull.
3. **Newton's method** — quadratic convergence of the MLE, error shown per iteration.
4. **Sample-size planner** — required n for a target standard error δ from the operational information.
5. **Risk dependence** — bias of the independence-based estimator under a Clayton copula, vs Kendall's τ, with Monte-Carlo uncertainty.
6. **Wald-test power** — exact non-centrality λ = n·hᵀ𝓘(θ₀)h for the q-dimensional hypothesis, with the direction-dependence range.

## Data input (CSV)

Upload your own data with three columns:

| column | meaning | values |
|--------|---------|--------|
| 1 — time | observed time (Z = min(X, T₀)) | number > 0 |
| 2 — cause | failure cause | 1 or 2 (any value ≠ 2 is read as 1) |
| 3 — status | event indicator | 1 = event, 0 = right-censored (optional; defaults to 1) |

Separators comma / semicolon / tab are accepted; a header row is auto-detected. Non-positive times and malformed rows are dropped. An example file and a template are available from the tool ("Download example CSV").

## Assumptions and scope

- Two competing causes; exponential or Weibull cause-specific intensities.
- The operational-information panels assume the horizon is set by the failure fraction: the fixed time **T₀ ≥ q_β** (the β-quantile), so stopping occurs at the r-th failure and τ\* = min(q_β, T₀) = q_β. If T₀ < q_β the study stops at T₀, the achieved failure fraction is below β, and those panels are outside their calibrated range.
- Latent cause times are assumed independent; Module 5 quantifies the bias when this is violated (exponential model, Clayton copula).
- Asymptotic (large-n) results; bootstrap CIs are provided as a small-sample cross-check and resample under the hybrid stopping rule.

## Tests

The core numerical routines are checked by an automated harness that **extracts the functions directly from `HybridCR.html`** (no duplicated source) and asserts closed-form reference values and edge cases:

```bash
node run_tests.cjs
```

It checks, among others: χ² quantiles (χ²₀.₉₅(1)=3.8415, χ²₀.₉₅(4)=9.4877), the central/non-central χ² identity, the closed-form operational variance, the exact Wald non-centrality and its eigenvalue range, positive-definiteness and β-monotonicity of 𝓘ₒₚ, MLE recovery of the truth on a seeded n=5000 sample (â₁≈1.01, â₂≈0.49), and edge cases (all-censored samples, malformed CSV rows, Clayton-copula bias zero at τ=0 and increasing in τ). Current status: **20/20 passing**. Stored reference values: `reference_values.json`.

## Cross-validation against R packages

Two cross-validation scripts reproduce HybridCR's MLE and cumulative incidence with independent code. `compare_python.py` (NumPy only, no R required) re-implements the simulation, the cause-specific MLE and the Aalen–Johansen estimator from scratch and checks them against HybridCR and against the fitted parametric CIF; `compare_R.R` does the same using the R packages `survival`, `cmprsk` and `etm`:

```bash
python3 compare_python.py     # no R needed
Rscript compare_R.R           # optional, if R is installed
```

## How to cite

If you use HybridCR, please cite the software (via its archived release — Mendeley Data DOI, to be minted on first release) together with the companion articles:

- N. S. Nurmukhamedova, G. S. Abdujalilova, M. Berdimuradov, et al. *A Bernstein–von Mises theorem for parametric competing risks under hybrid censoring.* **Mathematics** 14(18) (2026) 3277. https://doi.org/10.3390/math14183277
- N. S. Nurmukhamedova. *Simulation dataset on power and size of Wald, likelihood-ratio, and score test statistics in parametric competing risks models under hybrid censoring.* **Data in Brief** 68 (2026) 113136. https://doi.org/10.1016/j.dib.2026.113136

See `CITATION.cff` for machine-readable metadata.

## License

MIT — see `LICENSE`.
