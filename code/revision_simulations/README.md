# Revision Simulations

This directory contains the code for current and historical revision
simulations. Run all commands from the workflowr repository root.

## Current Analyses

| Analysis | Workflowr page | Simulation entry point |
|---|---|---|
| R3: two-stage dynamic discovery and functional testing | `analysis/revision_functional_testing_simulation.rmd` | `internal/r3_stage1_method_comparison/run_r3_stage1_method_comparison.R` (paired Stage 1); `internal/r3_two_stage_functional_testing_pilot10/run_r3_two_stage_pilot10.R` (baseline and Stage 2) |
| R4: unit-specific z-null and full-model residual correlations | `analysis/revision_correlated_error_simulation.rmd` | `r4_correlated_errors/run_unit_specific_residual_permutation.R`; `r4_correlated_errors/compute_full_model_residual_correlation.R`; `r4_correlated_errors/unit_specific_residual_permutation_helpers.R` |

R1/R2 are retired from the formal index. Their source pages are preserved in
`archived_pages/2026-09-17/`, and their original URLs link to R3 and historical
HTML snapshots. Their original experiment code, figures, and caches are retained.

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

Historical R1/R2 and current R3 use the versioned YRI `DS` dosage cache built by
`shared/build_real_genotype_one_per_gene_cache.R`: for each seed, one tested
cis variant is sampled uniformly from each of 6,362 tested genes. The same
ordered genotype matrix for the original five seeds is shared across R1, R2,
and both R3 truth mechanisms and audited by
`shared/validate_r1_r2_real_genotype_pairing.R`. R3 extends that frozen cache to
ten seeds using `internal/r3_two_stage_functional_testing_pilot10/` and retains
the first five samples exactly. Its original two-stage result cache is
`output/revision_simulations/mc/r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_qdyn005_qfunc005_fashr0143_pilot10/`.
The paired Stage-1 comparison cache is
`output/revision_simulations/mc/r3_stage1_iwp1_linear_direct_paired_raw_bf_fashr0143_pilot10/`.
It contains 20 paired replicates with both Raw and BF all-unit scores and
curves for IWP1 and FASH-linear, plus both direct LRTs. Main figures show the
two BF-adjusted FASH methods and the direct linear/quadratic tests. Direct-test
eFDR uses 100 donor permutations and the true null proportion `5090/6362`.
Stage 1 evaluates dynamic power/FDR across alpha; Stage 2 retains the original
IWP1 candidate sets and reports FSR, conditional classification power, and
calls at `q_dyn = q_func = 0.05`. The reporting gate reproduces the saved curves
from unit-level scores and verifies unchanged IWP1 BF and Stage-2 results.

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

The `internal/` directory includes the current formal R3 entry point listed
above and retains exploratory work for auditability:

- `pilots/`: parameter screens and truth previews;
- `diagnostics/`: calibration investigations, ablations, and summaries;
- `r1_real_genotype/`: the internal
  `analysis/revision_internal_real_genotype_simulation.rmd` experiment, with
  its real-dosage locus-block sampler, resumable driver, reporting helper, and
  cache at
  `output/revision_simulations/internal/r1_real_genotype_locus_blocks_pilot5/`;
- `archived_experiments/`: earlier analysis entry points superseded by R1-R3.

The exploratory scripts listed above are not current formal entry points. The
`internal/r3_two_stage_functional_testing_pilot10/` reporting module is required
to build the formal R3 page; the earlier R3 runner remains for auditability.
