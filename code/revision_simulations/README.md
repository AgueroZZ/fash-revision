# Revision Simulations

This directory contains the code for the four reviewer-facing revision
simulations. Run all commands from the workflowr repository root.

## Current Analyses

| Analysis | Workflowr page | Simulation entry point |
|---|---|---|
| R1: random B-spline effects | `analysis/revision_genotype_level_simulation.rmd` | `r1_random_bspline/run_random_bspline_mc_replication.R` |
| R2: spiky transient effects | `analysis/revision_combined_spiky_genotype_simulation.rmd` | `r2_spiky_transient/run_spiky_transient_mc_replication.R` |
| R3: functional testing | `analysis/revision_functional_testing_simulation.rmd` | `r3_functional_testing/run_matched_functional_testing_mc_replication.R` |
| R4: correlated expression errors | `analysis/revision_correlated_error_simulation.rmd` | `r4_correlated_errors/estimate_real_data_error_correlations.R`; `r4_correlated_errors/run_full_empirical_correlation_mc.R`; `r4_correlated_errors/run_lag1_correlation_sweep.R` |

The R2 folder also contains the single-replicate driver used while developing
the five-seed analysis.

R4 first estimates two complete 16 by 16 correlation patterns from the 500
most null-like one-pair-per-gene real-data trajectories, then inserts both
nearest-positive-definite matrices into paired R1 simulations. Its separate
lag-1 sweep reports IWP1 power, empirical FDR, and estimated pi0 before and
after BF adjustment. The earlier fixed-rho runner remains in the R4 directory
for auditability but is no longer the reviewer-page entry point.

## Shared Code

`shared/simulation_functions.R` contains the reusable data generation,
FASH fitting, direct-interaction testing, functional testing, summarization,
and plotting functions used across R1-R4.

Each analysis directory contains a `reporting.R` file. These scripts load and
validate the versioned cache and prepare report-only tables and plotting
objects. Keeping that code outside the R Markdown setup chunk lets the pages
show the scientifically relevant simulation and inference calls without
exposing routine report plumbing.

## Cached Builds

The workflowr pages read versioned results under
`output/revision_simulations/`. Building a page does not rerun a Monte Carlo
simulation. Each page contains an unevaluated command showing how to
regenerate its cache and evaluated chunks showing how every displayed table
and figure is produced from that cache.

## Internal Code

The `internal/` directory retains exploratory work for auditability:

- `pilots/`: parameter screens and truth previews;
- `diagnostics/`: calibration investigations, ablations, and summaries;
- `r1_real_genotype/`: the internal
  `analysis/revision_internal_real_genotype_simulation.rmd` experiment, with
  its real-dosage locus-block sampler, resumable driver, reporting helper, and
  cache at
  `output/revision_simulations/internal/r1_real_genotype_locus_blocks_pilot5/`;
- `archived_experiments/`: earlier analysis entry points superseded by R1-R3.

These scripts are not reviewer-facing entry points and are not required to
build the four retained workflowr pages.
