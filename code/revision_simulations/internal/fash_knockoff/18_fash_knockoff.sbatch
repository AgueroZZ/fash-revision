#!/bin/bash

#SBATCH --job-name=fash_knockoff
#SBATCH --account=pi-mstephens
#SBATCH --partition=caslake
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/project/mstephens/ziangzhang/fash/workspace/logs/out/%x_%j.out
#SBATCH --error=/project/mstephens/ziangzhang/fash/workspace/logs/err/%x_%j.err

set -euo pipefail

# Full-data FASH target-decoy (fash-knockoff) analysis.
#
# One global donor permutation is applied to the entire genotype matrix. For
# each of the 1,009,173 tested gene-variant pairs a permuted decoy unit is
# built, and target plus decoy units are fitted jointly so that a single
# empirical-Bayes prior scores both sides of every pair. The competition
# statistic is W = log BF - log BF_permuted, aggregated to the gene level as a
# maximum over the same variant set on both sides, and selection uses
# knockoff+.
#
# The target likelihood rows are NOT recomputed: they are taken verbatim from
# the source fit. BF_update() does not modify L_matrix, so those rows are
# bit-identical whether the raw or the BF-updated fit is supplied. Only the
# decoy units require new likelihood evaluations.
#
# Input resolution is delegated to the R preflight, which searches a documented
# candidate list and prints exactly what it resolved. Pass a project root as
# the first argument if the inputs live outside the search roots:
#
#   sbatch slurm/18_fash_knockoff.sbatch /path/to/project/root
#
# Runtime basis: the decoy likelihood measured 33.4 ms per unit on 8 local
# cores under fashr 0.1.43, i.e. about 0.27 core-seconds per unit, so roughly
# 75 core-hours for all 1,009,173 decoy units, or about 3-5 hours on 32 cores
# allowing for imperfect scaling. The 24-hour request is deliberate headroom,
# and the likelihood stage is resumable by chunk if it is killed anyway.
#
# Memory basis: the merged likelihood matrix is 2,018,346 x 52 (840 MB), and
# collapse_L_fash() materializes two further copies of that size, so the merge
# stage peaks near 6 GB on top of the loaded fit.
#
# The three stages run in three separate R processes. This is required, not
# cosmetic: the full FASH object is several gigabytes and forking likelihood
# workers from a session that has loaded it is what produced the anomalous
# runtimes recorded in earlier logs.

readonly PROJECT_ROOT="/project/mstephens/ziangzhang/fash"
readonly WORKSPACE="${PROJECT_ROOT}/workspace"
readonly CODE_DIR="${WORKSPACE}/code/revision_simulations/internal/fash_knockoff"

readonly R_SCRIPT="${CODE_DIR}/run_fash_knockoff.R"
readonly HELPERS="${CODE_DIR}/fash_knockoff_helpers.R"
readonly TESTS="${CODE_DIR}/test_fash_knockoff_helpers.R"
readonly DONOR_HELPERS="${WORKSPACE}/code/revision_simulations/internal/covariance_estimation/donor_null_permutation_helpers.R"
readonly SELECTED_HELPERS="${WORKSPACE}/code/revision_simulations/internal/selected_signal_genotype_permutation/selected_signal_genotype_permutation_helpers.R"
readonly SIMULATION_FUNCTIONS="${WORKSPACE}/code/revision_simulations/shared/simulation_functions.R"

readonly SEED="${FASH_KNOCKOFF_SEED:-20260823}"
readonly GENE_LIMIT="${FASH_KNOCKOFF_GENE_LIMIT:-0}"
if [[ "${GENE_LIMIT}" -gt 0 ]]; then
  readonly OUTPUT_ID="fash_knockoff_genes${GENE_LIMIT}_seed${SEED}"
else
  readonly OUTPUT_ID="fash_knockoff_full_seed${SEED}"
fi

# Bulk output (the target likelihood cache, the prepared decoy beta/SE, and the
# per-chunk decoy likelihood blocks) is about 1.1 GB and belongs in
# full_results, not in the SSHFS-mounted workspace.
readonly RESULT_PARENT="${PROJECT_ROOT}/full_results/application"
readonly FINAL_DIR="${RESULT_PARENT}/${OUTPUT_ID}"
readonly COMPLETE_FLAG="${FINAL_DIR}/complete.flag"

# Small diagnostic tables are copied here so they are visible on the mount.
readonly SUMMARY_DIR="${WORKSPACE}/results/application/fash_knockoff/${OUTPUT_ID}"

export FASH_KNOCKOFF_OUTPUT_PARENT="${RESULT_PARENT}"
export FASH_KNOCKOFF_SUMMARY_DIR="${SUMMARY_DIR}"

# Optional: an extra search root for the analysis inputs, highest priority.
if [[ "$#" -ge 1 && -n "${1:-}" ]]; then
  export FASH_KNOCKOFF_PROJECT_ROOT="$1"
fi

# Forked likelihood workers must not each spawn a BLAS thread pool.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p \
  "${WORKSPACE}/logs/out" \
  "${WORKSPACE}/logs/err" \
  "${RESULT_PARENT}" \
  "${SUMMARY_DIR}"

# Only the code is checked here. The analysis inputs are resolved and reported
# by the R preflight below, which fails with the full candidate list if any of
# them cannot be found.
missing_count=0
for required_file in \
  "${R_SCRIPT}" \
  "${HELPERS}" \
  "${TESTS}" \
  "${DONOR_HELPERS}" \
  "${SELECTED_HELPERS}" \
  "${SIMULATION_FUNCTIONS}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "MISSING required code file: ${required_file}" >&2
    missing_count=$((missing_count + 1))
  fi
done
if [[ "${missing_count}" -gt 0 ]]; then
  echo "${missing_count} required code file(s) are missing; see above." >&2
  exit 1
fi

# Resuming an interrupted run is intended; rerunning a finished one is not.
if [[ -f "${COMPLETE_FLAG}" ]]; then
  echo "Refusing to rerun a completed analysis: ${COMPLETE_FLAG}" >&2
  echo "Delete the flag, or submit with a different FASH_KNOCKOFF_SEED." >&2
  exit 1
fi

module load R/4.4.1

# find_workflowr_root() resolves against the working directory.
cd "${WORKSPACE}"

echo "=== fash-knockoff full-data target-decoy analysis ==="
echo "Seed:        ${SEED}"
echo "Gene limit:  ${GENE_LIMIT} (0 = all tested genes)"
echo "Extra root:  ${FASH_KNOCKOFF_PROJECT_ROOT:-<none>}"
echo "Slurm CPUs:  ${SLURM_CPUS_PER_TASK}"
echo "Bulk output: ${FINAL_DIR}"
echo "Summaries:   ${SUMMARY_DIR}"
echo

echo "--- Focused helper tests ---"
Rscript --vanilla "${TESTS}"
echo

echo "--- Package and input preflight ---"
Rscript --vanilla "${R_SCRIPT}" --stage prepare --seed "${SEED}" \
  --gene-limit "${GENE_LIMIT}" --preflight-only true
echo

if [[ -f "${FINAL_DIR}/prepared_decoy_input.rds" && \
      -f "${FINAL_DIR}/target_likelihood.rds" ]]; then
  echo "--- Stage 1/3 prepare: cached, skipping ---"
else
  echo "--- Stage 1/3 prepare ---"
  Rscript --vanilla "${R_SCRIPT}" --stage prepare --seed "${SEED}" \
    --gene-limit "${GENE_LIMIT}"
fi
echo

echo "--- Stage 2/3 decoy likelihood (resumable by chunk) ---"
Rscript --vanilla "${R_SCRIPT}" --stage likelihood --seed "${SEED}" \
  --gene-limit "${GENE_LIMIT}" --num-cores "${SLURM_CPUS_PER_TASK}" \
  --chunk-units 25000
echo

echo "--- Stage 3/3 merged empirical-Bayes refit and knockoff+ ---"
Rscript --vanilla "${R_SCRIPT}" --stage merge --seed "${SEED}" \
  --gene-limit "${GENE_LIMIT}"
echo

if [[ ! -f "${COMPLETE_FLAG}" ]]; then
  echo "The merge stage returned without creating ${COMPLETE_FLAG}." >&2
  exit 1
fi

echo "=== Completed: ${FINAL_DIR} ==="
cat "${COMPLETE_FLAG}"
