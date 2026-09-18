# R3 two-stage functional testing with ten seeds

This experiment is the production implementation for the revised formal R3
page. It retains the validated R3 data-generating process and evaluates:

1. BF-adjusted IWP1 dynamic discovery power and empirical FDR across nominal
   alpha from 0.005 through 0.20.
2. Functional classification at `q_dyn = q_func = 0.05`, summarized by empirical
   FSR, conditional classification power, and call counts.

The ten seeds are `12345, 22345, 32345, 42345, 52345, 62345, 72345, 82345,
92345, 102345`. The completed pilot-5 experiment already contains both truth
mechanisms for the first five seeds. Those ten records are validated and
migrated without recomputing them. Only the two mechanisms for the five added
seeds are generated and fitted.

## Exact genotype extension

The ordered pair keys were extracted once from the pair-summary RDS whose size
and MD5 match the frozen pilot-5 genotype-cache provenance. The resulting
4.9 MB artifact has SHA-256:

```text
a4b83e9828b0caeebffa226cc14100116599853a5c94f97ff4fcd2b8f8310145
```

Before extracting any new dosage, the builder resamples all five original seeds
from these keys and requires their complete selection tables to equal the
frozen cache exactly. It then streams the server VCF only for variants selected
by the five added seeds. The output is a separately named ten-seed genotype
cache; the five original sample objects must remain exactly identical to their
parent objects.

To recreate the compact pair-key artifact from the frozen pair-summary RDS:

```bash
Rscript --vanilla \
  code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10/build_pilot10_genotype_cache.R \
  --mode pair-source \
  --project-root . \
  --pair-summary <path-to-frozen-eqtl_summary.rds> \
  --pair-source output/revision_simulations/internal/r3_two_stage_pilot10_development_20260904/pair_key_source.rds
```

## Midway3 execution

The mounted workspace rules reserve Slurm submission and monitoring to the
user. After the files are staged, submit from Midway3:

```bash
cd /project/mstephens/ziangzhang/fash/workspace
sbatch slurm/run_r3_two_stage_pilot10.sbatch
```

The job requests the established `pi-mstephens` account and `caslake`
partition, 16 CPUs, 128 GB RAM, and 24 hours. It loads R 4.4.1, verifies the
exact `fashr` 0.1.43 Git SHA, builds or validates the new genotype cache, checks
the formal and pilot-5 manifests, migrates ten existing records, and computes
ten new records. Resubmission reuses compatible input, fit, and result
checkpoints. A completed result is never overwritten.

Large checkpoints are written below:

```text
/project/mstephens/ziangzhang/fash/full_results/revision_simulations/r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_qdyn005_qfunc005_fashr0143_pilot10/
```

The completed small result is written below:

```text
/project/mstephens/ziangzhang/fash/workspace/results/revision_simulations/r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_qdyn005_qfunc005_fashr0143_pilot10/
```

## Result contract

The cache contains 20 mechanism-by-seed records. Each record retains the full
6,362-unit BF FDR table, Stage-1 calls at 0.05, candidate-level functional lfsr,
Stage-2 calls, truth, and replicate-specific metrics. Summaries include:

- `stage1_alpha_by_seed.csv` and `stage1_mc_summary.csv` for the 40-point alpha
  grid;
- `functional_by_seed.csv` and `functional_mc_summary.csv` at the two 0.05
  thresholds;
- `stage2_random_bspline.csv` and `stage2_raised_cosine.csv`, each restricted to
  empirical FSR, conditional classification power, and mean calls;
- all-unit and all-candidate audit tables, genotype/input validation, exact
  source hashes, the session record, a manifest, and `complete.flag`.

The replicate records retain end-to-end power only as a non-displayed audit
quantity. It is not part of the formal Stage-2 tables or interpretation.

## Local checks and formal-page promotion

From `coderepo-local`, run:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla \
  code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10/test_pilot10_helpers.R

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  Rscript --vanilla \
  code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10/test_genotype_cache_helpers.R
```

`revision_functional_testing_simulation_candidate.rmd` is held beside the code
until the production cache is complete. After copying the completed result and
genotype cache locally, `reporting.R` reconstructs every Stage-1 alpha curve and
Stage-2 summary from the 20 records, validates all artifact hashes, and checks
that the fitted genotype cache is identical. Only then should the candidate
replace `analysis/revision_functional_testing_simulation.rmd` and be rendered as
the index-linked formal R3 page.
