# Internal R3 two-stage functional testing

This experiment retains the formal R3 truth and observation settings and changes
only the inference sequence: fit all 6,362 units, select BF-adjusted dynamic
discoveries at cumulative FDR 0.05, then select each functional at cumulative
FSR 0.05 within that candidate set. The formal R3 analysis remains unchanged.

The primary tables report five-seed means of empirical FSR, end-to-end power,
and call counts. Power uses all true positive functionals before screening.
Conditional classification power is a secondary diagnostic. A zero conditional
denominator gives `NA`, with the number of contributing seeds reported.

## Run on Midway3

The mounted workspace rules reserve submission and monitoring to the user.
After staging, submit from the server:

```bash
cd /project/mstephens/ziangzhang/fash/workspace
sbatch slurm/run_r3_two_stage_functional_testing.sbatch
```

The job uses the existing R3 account/partition, R 4.4.1, 16 CPUs, 128 GB RAM,
and a 24-hour limit. It verifies fashr 0.1.43 and its exact Git SHA before
computing. The formal posterior sampler uses 3,000 draws on the 0.1-day grid.

Result ID:

```text
r3_two_stage_functional_testing_iwp1_mixture_qdyn005_qfunc005_fashr0143_pilot5
```

Small results are saved below `workspace/results/revision_simulations/internal/`.
Large input/fit checkpoints are saved below
`full_results/revision_simulations/internal/`, outside the SSHFS workspace.
Both locations use the new result ID. An incomplete result uses a `.partial`
suffix; resubmitting the same script reuses compatible input, fit and completed
replicate checkpoints. Completed final results are never overwritten.

## Provenance and reuse

- Reuse the original sampled YRI genotype cache, frozen R3 configuration,
  simulation functions and all component seeds.
- The original ten formal replicate files do not contain full fits or all-unit
  posterior summaries. Reconstruct missing inputs rather than treating the
  original selected-call table as a full-universe posterior cache.
- Before fitting, match genotype, true-beta, regression-beta-hat, and adjusted-SE
  digests against the original-regression records retained by the ideal-Gaussian
  experiment. This reuses its validation records, not its ideal observations.
- Fit only IWP1 and apply its BF update; omit unused linear comparators and
  additional alpha evaluations. Require the BF prior null weight and dynamic
  discovery count to reproduce formal R3.
- Reuse the frozen posterior function and component seed, and recalculate
  cumulative FSR within the Stage-1 subset. Since only that subset is sampled,
  posterior Monte Carlo realizations need not be identical to the previous
  full-universe calls. The sampling distribution and draw count are unchanged.
- Retain new input matrices and raw/BF fits in server-only checkpoints. Small
  replicate records retain all-unit FDR tables and truth, every candidate lfsr,
  final selections, metrics, source/input checks and reuse information.

## Small result artifacts

- `configuration.rds`, `manifest.rds`, `sessionInfo.txt`, `complete.flag`.
- `replicates/<mechanism>_seed_<seed>.rds`: auditable per-seed inference records.
- `summary_tables.rds`: all CSV tables in one object.
- `summary/stage1_by_seed.csv` and `stage1_mc_summary.csv`.
- `summary/functional_by_seed.csv` and `functional_mc_summary.csv`.
- `summary/primary_random_bspline.csv` and `primary_raised_cosine.csv`: the two
  three-row primary tables, with numeric rates and mean call counts.
- `summary/all_stage1_units.csv`: all-unit lfdr/FDR, selected status and truth.
- `summary/all_stage2_candidates.csv`: candidate lfsr, recalculated cumulative
  FSR, final selections and actual functional truth.
- `summary/input_validation.csv`: four digest checks and full-data fit agreement.

## Local validation and reporting

Run from the local `coderepo-local` directory:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla code/revision_simulations/internal/r3_two_stage_functional_testing/test_two_stage.R
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript --vanilla code/revision_simulations/internal/r3_two_stage_functional_testing/run_r3_two_stage_functional_testing.R --mode preflight --num-cores 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript --vanilla code/revision_simulations/internal/r3_two_stage_functional_testing/smoke_two_stage.R
```

The smoke script is one reduced process with 60 units, at most two fitting
workers and one posterior worker. Its results are explicitly excluded from
production reporting. Full production execution is refused on macOS.

After the server cache completes, copy only the small result directory into
local `output/revision_simulations/internal/`, preserving its result ID and
manifest. Do not copy server-only fit checkpoints. Then render:

```r
workflowr::wflow_build(
  "analysis/revision_internal_r3_two_stage_functional_testing.rmd",
  view = FALSE
)
```

The page is absent from the workflowr index. Before production completion it
shows an explicit pending status and no numerical results. A present but invalid
final cache causes rendering to fail. Inspect the six tables, displayed values,
code folding, math, and table overflow after rendering the completed results.
