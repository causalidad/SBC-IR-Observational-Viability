# SBC IR Observational Viability
Repository Version: v1.0-observational-benchmark
This repository accompanies the manuscript:

**Observational Viability of an IR-Safe Perturbative Cosmological Sector Inspired by the Stochastic Bidirectional Causality Framework**

## Purpose

This repository documents the numerical material used to characterize an observationally viable infrared-safe perturbative cosmological sector inspired by the SBC framework.

The analysis should be interpreted as a phenomenological observational-viability study, not as a full cosmological parameter-estimation pipeline and not as a first-principles derivation of SBC.

## Main observational benchmark

Proxy-posterior median values:

- epsilon = 0.068985
- k0 = 0.023872
- n = 0.658078

16th--84th percentile ranges:

- epsilon: 0.061998 -- 0.083889
- k0: 0.017374 -- 0.032180
- n: 0.238041 -- 0.897020

Background deformation:

- alpha_bg ≈ 0

## MCMC diagnostics

Integrated autocorrelation times:

- epsilon: 25.887822
- k0: 21.709094
- n: 23.151373

Effective samples:

- epsilon: 370.830734
- k0: 442.210995
- n: 414.662240

## Repository structure

- `paper/`: manuscript source files.
- `scripts/`: scripts used to generate posterior, figures, and diagnostics.
- `chains/`: proxy-posterior MCMC chain.
- `figures/`: final manuscript figures.
- `reports/`: posterior summaries and diagnostic reports.
- `results/`: numerical CSV outputs.
- `docs/`: notes and future research directions.

## Scientific status

The results establish observational viability of the adopted effective infrared parametrization. They do not demonstrate observational preference over LambdaCDM, Bayesian evidence, or a first-principles derivation from SBC.

## Open question

Can projection-based causal-memory constructions, such as those motivated by Mori--Zwanzig approaches, naturally reproduce the observationally viable infrared regime identified here?

Citation and Archive

The reproducibility repository associated with this work is publicly available on GitHub and permanently archived on Zenodo.

Zenodo DOI: 10.5281/zenodo.20720308

This DOI corresponds to the archived release of the observational benchmark, including scripts, MCMC chains, figures, reports, and reproducibility materials used in the analysis.

