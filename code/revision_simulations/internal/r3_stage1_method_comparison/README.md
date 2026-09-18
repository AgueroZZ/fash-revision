# Paired R3 Stage-1 method comparison

This experiment adds matched FASH-linear and individual-level direct linear
and quadratic interaction tests to the existing ten-seed two-stage R3
simulation. The completed R3 truth, selected gene-variant pairs, observations,
IWP1 BF decisions, and Stage-2 classification are retained.

## Scientific contract

- Twenty paired replicates: ten seeds under each of the two R3 truth mechanisms.
- All 6,362 units are tested for dynamic status before any Stage-2 screening.
- The main figures compare BF-adjusted IWP1, BF-adjusted FASH-linear, direct
  linear LRT, and the joint linear-plus-quadratic LRT.
- Both FASH models retain Raw and BF all-unit lfdr/cumulative-FDR scores,
  prior weights, and curves over all 40 alpha levels. Adding Raw-versus-BF
  comparisons later requires only reporting changes.
- The direct tests share 100 donor permutations per replicate, use the original
  `seed + 10000L` rule, permute expression and covariate rows together, and use
  the known dynamic-null proportion `5090 / 6362` for permutation eFDR. These
  q-values are not ordinary BH adjustment.
- Power and empirical FDR are means of per-replicate rates. Ribbons are
  pointwise replicate ranges, not confidence intervals.
- Stage 2 remains the original BF-IWP1 candidate-restricted analysis at
  `q_dyn = q_func = 0.05`.

The comparator source is the frozen R1/R2 fashr 0.1.43 snapshot, sourced in a
separate environment. No shared scientific helper is edited. R3 inputs are
read from their retained server checkpoints when available. Otherwise, only
the inputs are reconstructed using the frozen R3 generator and original
seeds. In either case the runner reproduces all four retained input digests
and checks the expression and covariates against their original generation.
Original IWP1 Raw/BF fits must remain available on the server: the runner
extracts their scores and checks BF equality without refitting them.

## Local checks

From the workflowr repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla code/revision_simulations/internal/r3_stage1_method_comparison/test_comparison.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla code/revision_simulations/internal/r3_stage1_method_comparison/smoke_comparison.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla code/revision_simulations/internal/r3_stage1_method_comparison/run_r3_stage1_method_comparison.R \
  --mode preflight --num-cores 2
```

Smoke results contain only 60 units and three permutations per truth, run
sequentially in one R process with at most two workers. They are explicitly
excluded from production reporting. Full production refuses to run outside
a Linux Slurm job or with an R version other than the original R 4.4.1.

## Manual Midway3 submission

The mounted workspace rules reserve job submission and monitoring to the user:

```bash
cd /project/mstephens/ziangzhang/fash/workspace
sbatch slurm/run_r3_stage1_method_comparison.sbatch
```

The script requests account `pi-mstephens`, partition `caslake`, 16 CPUs,
128 GB, and 24 hours. Full fits and permutation null matrices remain below
`/project/mstephens/ziangzhang/fash/full_results/revision_simulations/r3_stage1_iwp1_linear_direct_paired_raw_bf_fashr0143_pilot10/`.
Compact results are saved below the mounted workspace at
`results/revision_simulations/r3_stage1_iwp1_linear_direct_paired_raw_bf_fashr0143_pilot10/`.

The runner checkpoints the retained IWP1 scores, new linear fits, and direct
permutations independently and reuses completed paired replicates. An
interrupted permutation batch restarts only that batch. A `.partial` result
is never rendered. Resubmission requires the same configuration and source
revision. An empty `.running` directory prevents concurrent writers; after
an abnormal scheduler termination, verify the job has ended before removing
only that empty directory with `rmdir`. Completed results are never overwritten.

## Result schema and promotion

Each compact replicate contains six named all-unit score tables, pair keys,
truth, input digests, Raw/BF prior weights, direct LRT and eFDR threshold
diagnostics, donor permutation indices, and full-checkpoint paths and hashes.
The cache contains 4,800 per-seed curve rows and 480 Monte Carlo summary rows;
the four-method display selects 3,200 and 320 of those rows, respectively.

After the manual production run, copy only the completed compact result to
`output/revision_simulations/mc/` with the same result ID and byte-verify every
file. `r3c_load_report()` validates source/configuration/record hashes,
reconstructs all curves from saved all-unit scores, verifies IWP1 BF equality,
and requires exact identity of the Stage-2 tables to the original cache.

Production job `59259939` completed on 2026-09-17 at 23:41 America/Chicago.
All 20 records passed the reporting gate, and 32 compact files were imported
byte-for-byte. The candidate was promoted to
`analysis/revision_functional_testing_simulation.rmd` after validation.
R1/R2 were removed from the formal index; their source pages are archived in
`code/revision_simulations/archived_pages/2026-09-17/`, and their original URLs
link to R3 and byte-identical historical HTML snapshots. All original figures
and result caches are retained. The rebuilt pages passed desktop and mobile
checks for figures, captions, code folding, and rendered mathematics.

See the paired workspace plan/log dated 2026-09-17 for validation evidence,
actual method comparisons, and the browser-console limitation.

## Six representative trajectories

Figure 1 shows one broad example for each temporal category and one transient
example for each generating peak count. Points and plus/minus two reported
standard-error bars accompany the orange dense true curve. Examples are ordered
by unit index without selecting on discoveries.

The original seed-12345 illustration pool lacks two-peak units. Run
`Rscript --vanilla code/revision_simulations/internal/r3_stage1_method_comparison/prepare_two_peak_example.R`
once with one CPU thread to recover the first two-peak unit using the frozen
input generator. No FASH fitting, direct testing, or permutations are repeated.
All unit-level true functionals and all 12 stored compact-example curves,
estimates, and standard errors must agree within 1e-12. The cache preserves both
original and reconstructed input digests because numerical agreement across
Mac/Linux does not imply byte-identical floating-point matrices. The reporting
gate validates this separate cache before plotting.

## Peak-stratified compact-truth power

Figure 3 reuses the global adjusted scores and separates true dynamic units
by their generating raised-cosine component count. Recover the original labels
once, without expression generation or fitting, using:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla code/revision_simulations/internal/r3_stage1_method_comparison/prepare_peak_stratification.R
```

The immutable label cache is stored in
`output/revision_simulations/internal/r3_peak_stratified_power_20260918/`.
The producer checks frozen source hashes, seeds, genotype order, all per-unit
true functionals (absolute tolerance 1e-12), and all retained original example
curves and peak counts. It records both original and reconstructed beta hashes:
Mac/Linux floating-point tails prevent bitwise beta equality. See the paired
2026-09-18 plan/log for measured discrepancies and validation limits.

`r3c_peak_power()` validates labels against each current replicate, keeps all
six retained methods, and requires weighted subgroup power to reproduce every
original compact-truth curve. Main figures select the four existing display
methods. Power is averaged across seeds after using the stratum's dynamic-unit
count as the denominator in each seed. No within-stratum FDR adjustment is made.

## Temporal-stratified power under both truth mechanisms

Figure 4 uses original baseline true-functional values to partition dynamic
units into Early/Middle/Late (369/534/369 per seed and mechanism), including
both Switch statuses. Rows show broad B-spline and compact raised-cosine truth.
`r3c_location_power()` uses the original global scores and verifies that
count-weighted strata reproduce all original powers separately by mechanism.
No reconstructed labels are required for this grouping. All six methods
remain in expanded derived tables under
`output/revision_simulations/internal/r3_location_stratified_power_both_truths_20260918/`;
the earlier compact-only tables are preserved in their original directory.
The page displays four methods. This is Stage-1 discovery power by truth category.
