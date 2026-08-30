#!/bin/bash

#SBATCH --job-name=fash_r3_ideal_gaussian
#SBATCH --account=pi-mstephens
#SBATCH --partition=caslake
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --output=/project/mstephens/ziangzhang/fash/workspace/logs/out/%x_%j.out
#SBATCH --error=/project/mstephens/ziangzhang/fash/workspace/logs/err/%x_%j.err

set -euo pipefail

readonly PROJECT_ROOT="/project/mstephens/ziangzhang/fash"
readonly WORKSPACE="${PROJECT_ROOT}/workspace"
readonly R_SCRIPT="${WORKSPACE}/scripts/run_r3_ideal_gaussian_measurement.R"
readonly EXPERIMENT_CODE="${WORKSPACE}/code/revision_simulations/r3_ideal_gaussian_measurement"
readonly HELPER="${EXPERIMENT_CODE}/ideal_gaussian_measurement.R"
readonly FROZEN_R3_CODE="${WORKSPACE}/code/revision_simulations/r3_full_universe_functional_fashr0143"
readonly SIMULATION_FUNCTIONS="${FROZEN_R3_CODE}/r3_full_universe_functional_simulation_functions.R"
readonly REAL_GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
readonly GENOTYPE_CACHE="${WORKSPACE}/inputs/r1_r2_fashr0143/genotype_samples.rds"
readonly RESULT_PARENT="${WORKSPACE}/results/revision_simulations"
readonly RESULT_ID="r3_ideal_gaussian_known_t_adjusted_se_matched_truth_open_middle_3_12_full_universe_fashr0143_pilot5"
readonly FINAL_DIR="${RESULT_PARENT}/${RESULT_ID}"
readonly COMPLETE_FLAG="${FINAL_DIR}/complete.flag"

mkdir -p \
  "${WORKSPACE}/logs/out" \
  "${WORKSPACE}/logs/err" \
  "${RESULT_PARENT}"

for required_file in \
  "${R_SCRIPT}" \
  "${HELPER}" \
  "${SIMULATION_FUNCTIONS}" \
  "${REAL_GENOTYPE_HELPER}" \
  "${GENOTYPE_CACHE}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Required file is missing: ${required_file}" >&2
    exit 1
  fi
done

if [[ -e "${FINAL_DIR}" ]]; then
  echo "Refusing to overwrite completed output: ${FINAL_DIR}" >&2
  exit 1
fi

module load R/4.4.1

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "R3 ideal-Gaussian entry point: ${R_SCRIPT}"
echo "Observation model: beta_hat = true_beta + adjusted_se * iid N(0,1)"
echo "SE source: final t-adjusted R3 SE, regenerated and then treated as fixed"
echo "Truth mechanisms: random_bspline, raised_cosine"
echo "Middle definition: 3 < t < 12"
echo "Functional candidates: all 6362 units"
echo "Slurm CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Final output: ${FINAL_DIR}"

echo "Running package, source, seed, and input preflight."
FASH_R3_IDEAL_HELPER="${HELPER}" \
FASH_R3_SIMULATION_FUNCTIONS="${SIMULATION_FUNCTIONS}" \
FASH_R3_REAL_GENOTYPE_HELPER="${REAL_GENOTYPE_HELPER}" \
FASH_R3_GENOTYPE_CACHE="${GENOTYPE_CACHE}" \
FASH_R3_IDEAL_RESULT_PARENT="${RESULT_PARENT}" \
  Rscript --vanilla "${R_SCRIPT}" --preflight-only true

echo "Starting or resuming the ten ideal-Gaussian R3 replicates."
FASH_R3_IDEAL_HELPER="${HELPER}" \
FASH_R3_SIMULATION_FUNCTIONS="${SIMULATION_FUNCTIONS}" \
FASH_R3_REAL_GENOTYPE_HELPER="${REAL_GENOTYPE_HELPER}" \
FASH_R3_GENOTYPE_CACHE="${GENOTYPE_CACHE}" \
FASH_R3_IDEAL_RESULT_PARENT="${RESULT_PARENT}" \
  Rscript --vanilla "${R_SCRIPT}"

if [[ ! -f "${COMPLETE_FLAG}" ]]; then
  echo "The runner returned without creating ${COMPLETE_FLAG}." >&2
  exit 1
fi

echo "R3 ideal-Gaussian measurement experiment completed: ${FINAL_DIR}"
