# Power Analysis for Secondary Analysis in No-Report Paradigms

Analytic and Monte Carlo tools for planning the statistical power of the
**secondary analysis** in two-stage no-report designs, where momentary percepts
are first decoded from a physiological measure and a second neural variable is
then contrasted between the decoder-inferred percepts.

The code accompanies a methodological review of no-report paradigms in
multistable perception. It lets you convert an anticipated decoding accuracy and
a physical acquisition design into a required sample size, and it stress-tests
the analytic approximation against simulation when its assumptions are violated.

---

## The design being modelled

- **Stage 1 (primary analysis).** Each trial/sample is classified as percept A or
  B from a physiological measure (e.g. optokinetic nystagmus, pupil dilation,
  M/EEG, fMRI) with accuracy `acc`.
- **Stage 2 (secondary analysis).** A secondary measure (e.g. BOLD amplitude,
  band power) is compared between the classifier-labelled percepts. For each of
  `n` participants a mean difference `D_hat` is formed over `m` trials, and a
  group-level one-sample *t* test asks whether the mean `D_hat` differs from zero.

Power depends on the true effect, the classification accuracy, the number of
subjects and trials, between-subject heterogeneity, and — in continuous viewing —
on two distinct *effective* sample counts (see below).

---

## Requirements

- **R** (≥ 4.0 recommended); base R only for the statistics and plotting.
- Packages: **`dplyr`**, **`tidyr`** (used for the summary-table pipelines and
  for assembling p-value data frames). `ggplot2` is *not* required — all plots
  are base R.

```r
install.packages(c("dplyr", "tidyr"))
```

---

## Repository contents

| File | Description |
|------|-------------|
| `power_analysis_noreport.R` | The complete pipeline: analytic functions, Monte Carlo simulators, validation runs, plots, and summary tables. |
| `results_section.tex` *(optional)* | LaTeX results section reporting the tables and figures. |
| `power_analysis_summary.tex` *(optional)* | Standalone LaTeX summary with a worked study-planning example. |

Running the script writes the figure PDFs (`plot_validation.pdf`,
`plot_between_subject.pdf`, `plot_violations.pdf`, `plot_continuous_mclass.pdf`,
`plot_pvalue_qq.pdf`, `plot_pvalue_hist.pdf`, `plot_stress_qq.pdf`,
`plot_stress_hist.pdf`) to the working directory and prints the summary tables to
the console.

---

## Quick start

> **Before running:** the top of the script contains a hard-coded
> `setwd("/Users/...")`. **Edit or delete that line** so plots are written where
> you want them.

```r
source("power_analysis_noreport.R")   # runs the full study: validation + plots + tables
```

To use the functions on their own without running the whole pipeline, source the
file (or copy the function definitions) and call them directly:

```r
# Analytic power for a trial-based design
power_analytic(n = 20, m = 100, acc = 0.75, delta = 0.5, sigma_within = 1)

# Smallest n reaching 0.8 power given two effective counts (continuous viewing)
required_n_2comp(m_class = 200, m_noise = 294, acc = 0.70,
                 delta = 0.5, sigma_within = 1, sigma_between = 0.2)

# Plan a continuous-viewing study from physical inputs
plan_continuous(session_min = 20, mean_dom_s = 6, sampling_hz = 1,
                tau_noise_s = 2, acc = 0.70, d = 0.5, tau = 0.2)
```

---

## The two effective counts (continuous viewing)

In continuous viewing the samples are autocorrelated, so a single "effective
sample size" splits into **two** counts that shrink independently:

- **`m_class`** — the number of *independent classification decisions*. Bounded
  above by the number of dominance periods (`session_length / mean_dominance`)
  and set within that by how fast the decoder's *errors* decorrelate. Governs the
  classification-attenuation term.
- **`m_noise`** — the number of *independent noise samples* in the secondary
  measure, governed by that measure's own autocorrelation:
  `m_noise = m·(1 − ρ)/(1 + ρ)`. Governs the measurement-variance term.

Because the number of independent noise samples is set by the noise *correlation
time*, not by how densely you poll, specify the noise by **`tau_noise_s`**
(seconds) rather than a fixed lag-1 `ρ`. The helper `rho_from_tau()` converts:
`ρ = exp(−1 / (sampling_hz · tau_noise_s))`. With a correlation time, `m_noise`
is invariant to sampling rate above the resolving rate (→ `T / (2·tau_noise)`);
with a fixed `ρ` it scales spuriously with the sampling rate.

---

## Function reference

**Analytic power and planning**

| Function | Purpose |
|----------|---------|
| `power_analytic(n, m, acc, delta, sigma_within, sigma_between, alpha)` | Power for a trial-based design (single effective count `m`). |
| `power_analytic_2comp(n, m_class, m_noise, acc, delta, sigma_within, sigma_between, alpha)` | Power with the two continuous-viewing counts; reduces to `power_analytic` when `m_class = m_noise = m`. |
| `required_n_2comp(m_class, m_noise, acc, delta, sigma_within, sigma_between, target_power, alpha, n_max)` | Smallest `n` reaching `target_power`; `NA` if unreachable within `n_max`. |
| `rho_from_tau(tau_noise_s, sampling_hz)` | Convert a noise correlation time (s) to the AR(1) lag-1 `ρ`. |
| `plan_continuous(session_min, mean_dom_s, sampling_hz, tau_noise_s / rho, acc, d, tau, target_power, ...)` | End-to-end study planner: physical inputs → effective counts → required `n`, with conservative/optimistic `m_noise` brackets. Vectorised over any argument. |

**Monte Carlo simulators** (each returns power, or the full p-value vector with
`return_pvalues = TRUE`)

| Function | Scenario |
|----------|----------|
| `power_simulate(...)` | Baseline: independent classification errors; optional class imbalance (`balanced = FALSE`). |
| `power_simulate_correlated(...)` | Errors concentrated on **ambiguous** trials near the decision boundary (benign). |
| `power_simulate_correlated_far(...)` | Errors concentrated on **confident** trials far from the boundary (maximally damaging). |
| `power_simulate_autocorr(...)` | Continuous viewing, **one decision per dominance period** (the reference the formula assumes). |
| `power_simulate_persample(...)` | Continuous viewing, realistic **per-sample** decoder on a distinct feature channel; `m_class` emerges from the feature timescale. |
| `measure_persample(...)` | Diagnostic: reports emergent accuracy, number of periods, and the effective `m_class`/`m_noise` for a per-sample decoder. |

The continuous-viewing simulators accept non-Gaussian secondary-measure noise
(`noise_dist = "gaussian" | "t" | "lognorm"`), generated via a Gaussian-copula
(NORTA) construction so that the marginal shape is exact while the autocorrelation
is retained. They also accept `tau_noise_s` (and `feat_tau_s`) as sample-rate-
invariant alternatives to `rho`/`feat_rho`.

---

## What the pipeline demonstrates

1. **Validation** — the analytic approximation matches simulation across `n`,
   `m`, and accuracy (largest discrepancy ≈ 3×10⁻³).
2. **Between-subject variability** — heterogeneity `τ` imposes a power ceiling
   that more trials cannot overcome; only more subjects help.
3. **Assumption violations** — the approximation is robust to independent errors
   and class imbalance, but *structured* errors matter: near-boundary errors are
   benign, far-boundary (confident-but-wrong) errors sharply reduce power.
4. **Continuous viewing** — the effective classification count varies widely with
   the feature timescale, but power barely does, so the simple
   "`m_class` = number of periods" rule is adequate for sizing studies.
5. **Null calibration** — the group test yields Uniform(0,1) p-values under
   Gaussian noise; under heavy-tailed/skewed noise it stays calibrated *except*
   at low effective count, where a robust or permutation-based test is advisable.

---

## Planning a study (recipe)

1. **Fix the effect:** standardized per-trial `d = Δ/σ` and between-subject SD
   `τ`. Both are best bracketed rather than trusted to a small pilot.
2. **Anticipate accuracy** `acc` from the decoder (deflate optimistic pilot
   values; report speeds reversals).
3. **Convert to effective counts** with `plan_continuous()`
   (`m_class ≈ T/mean_dom`, `m_noise` from `tau_noise_s`).
4. **Solve for `n`** and read the conservative/optimistic bracket; sweep `d` and
   `τ` to check the required `n` is stable.

---

## Notes and caveats

- **`setwd()` at the top is machine-specific** — remove it before running.
- The estimator is a group one-sample *t* test on per-subject mean differences;
  at very low effective count with non-Gaussian noise, prefer a robust or
  permutation-based group test.
- The single-`τ` (AR(1)) noise model is cleanest for eye/M-EEG and only
  approximate for BOLD (the HRF makes the residual spectrum non-AR(1)); for BOLD,
  estimate `m_noise` from the residual integrated autocorrelation time directly.
- Simulations are stochastic; set a seed for reproducibility (the script does so
  per section).

---

## Citation

If you use this code, please cite the accompanying paper. *(Add citation / DOI
here.)*

## License

*(Add a license, e.g. MIT, here.)*
