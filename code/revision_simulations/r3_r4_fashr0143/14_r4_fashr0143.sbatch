#!/bin/bash

#SBATCH --job-name=fash_r4_0143
#SBATCH --account=pi-mstephens
#SBATCH --partition=caslake
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/project/mstephens/ziangzhang/fash/workspace/logs/out/%x_%j.out
#SBATCH --error=/project/mstephens/ziangzhang/fash/workspace/logs/err/%x_%j.err

set -euo pipefail

readonly PROJECT_ROOT="/project/mstephens/ziangzhang/fash"
readonly WORKSPACE="${PROJECT_ROOT}/workspace"
readonly R_SCRIPT="${WORKSPACE}/scripts/14_r4_fashr0143.R"
readonly SIMULATION_FUNCTIONS="${WORKSPACE}/code/revision_simulations/shared/simulation_functions.R"
readonly REAL_GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
readonly R4_R1_HELPER="${WORKSPACE}/code/revision_simulations/r4_correlated_errors/real_genotype_r1_helpers.R"
readonly FULL_DRIVER="${WORKSPACE}/code/revision_simulations/r4_correlated_errors/run_full_empirical_correlation_mc.R"
readonly SWEEP_DRIVER="${WORKSPACE}/code/revision_simulations/r4_correlated_errors/run_lag1_correlation_sweep.R"
readonly PLOTTING="${WORKSPACE}/code/revision_simulations/r4_correlated_errors/plotting.R"
readonly GENOTYPE_CACHE="${WORKSPACE}/inputs/r1_r2_fashr0143/genotype_samples.rds"
readonly MATRIX_CACHE="${WORKSPACE}/inputs/r3_r4_fashr0143/simulation_correlation_matrices.rds"
readonly R1_CACHE="${WORKSPACE}/results/revision_simulations/r1_r2_fashr0143"
readonly RESULT_PARENT="${WORKSPACE}/results/revision_simulations"
readonly FINAL_DIR="${RESULT_PARENT}/r4_fashr0143"
readonly COMPLETE_FLAG="${FINAL_DIR}/complete.flag"

mkdir -p \
  "${WORKSPACE}/logs/out" \
  "${WORKSPACE}/logs/err" \
  "${RESULT_PARENT}"

for required_file in \
  "${R_SCRIPT}" \
  "${SIMULATION_FUNCTIONS}" \
  "${REAL_GENOTYPE_HELPER}" \
  "${R4_R1_HELPER}" \
  "${FULL_DRIVER}" \
  "${SWEEP_DRIVER}" \
  "${PLOTTING}" \
  "${GENOTYPE_CACHE}" \
  "${MATRIX_CACHE}" \
  "${R1_CACHE}/manifest.rds" \
  "${R1_CACHE}/complete.flag"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Required file is missing: ${required_file}" >&2
    exit 1
  fi
done

if [[ -e "${FINAL_DIR}" ]]; then
  echo "Refusing to overwrite completed R4 output: ${FINAL_DIR}" >&2
  exit 1
fi

module load R/4.4.1

echo "R4 wrapper: ${R_SCRIPT}"
echo "Formal genotype cache: ${GENOTYPE_CACHE}"
echo "R4 matrix cache: ${MATRIX_CACHE}"
echo "R1 reference cache: ${R1_CACHE}"
echo "Slurm CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Final output: ${FINAL_DIR}"

echo "Running the exact R4 package, source, and input preflight."
FASH_WORKSPACE="${WORKSPACE}" \
  FASH_R4_RESULT_PARENT="${RESULT_PARENT}" \
  Rscript --vanilla "${R_SCRIPT}" --preflight-only true

echo "Starting the resumable formal R4 refresh."
FASH_WORKSPACE="${WORKSPACE}" \
  FASH_R4_RESULT_PARENT="${RESULT_PARENT}" \
  Rscript --vanilla "${R_SCRIPT}"

if [[ ! -f "${COMPLETE_FLAG}" ]]; then
  echo "The R4 wrapper returned without creating ${COMPLETE_FLAG}." >&2
  exit 1
fi

echo "R4 fashr 0.1.43 refresh completed: ${FINAL_DIR}"
