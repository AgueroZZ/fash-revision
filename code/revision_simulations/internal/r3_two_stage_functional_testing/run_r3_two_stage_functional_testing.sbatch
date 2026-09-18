#!/bin/bash

#SBATCH --job-name=fash_r3_two_stage
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
readonly R_SCRIPT="${WORKSPACE}/scripts/run_r3_two_stage_functional_testing.R"
readonly RESULT_ID="r3_two_stage_functional_testing_iwp1_mixture_qdyn005_qfunc005_fashr0143_pilot5"

export FASH_R3_TWO_STAGE_CODE="${WORKSPACE}/code/revision_simulations/internal/r3_two_stage_functional_testing"
export FASH_R3_SIMULATION_FUNCTIONS="${WORKSPACE}/code/revision_simulations/r3_full_universe_functional_fashr0143/r3_full_universe_functional_simulation_functions.R"
export FASH_R3_REAL_GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
export FASH_R3_DIGEST_HELPER="${WORKSPACE}/code/revision_simulations/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R"
export FASH_R3_GENOTYPE_CACHE="${WORKSPACE}/inputs/r1_r2_fashr0143/genotype_samples.rds"
export FASH_R3_FORMAL_DIR="${WORKSPACE}/results/revision_simulations/r3_real_genotype_one_per_gene_J6362_matched_functional_open_middle_3_12_center_aligned_iwp1_geometry_mixture_relative_location_clearance_full_universe_paired_posterior_fashr0143_pilot5"
export FASH_R3_IDEAL_DIR="${WORKSPACE}/results/revision_simulations/r3_ideal_gaussian_known_t_adjusted_se_matched_truth_open_middle_3_12_center_aligned_iwp1_geometry_mixture_full_universe_fashr0143_pilot5"
export FASH_R3_TWO_STAGE_RESULT_PARENT="${WORKSPACE}/results/revision_simulations/internal"
export FASH_R3_TWO_STAGE_CHECKPOINT_ROOT="${PROJECT_ROOT}/full_results/revision_simulations/internal/${RESULT_ID}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p "${WORKSPACE}/logs/out" "${WORKSPACE}/logs/err"
for required in \
  "${R_SCRIPT}" \
  "${FASH_R3_TWO_STAGE_CODE}/two_stage_helpers.R" \
  "${FASH_R3_TWO_STAGE_CODE}/reconstruct_r3.R" \
  "${FASH_R3_SIMULATION_FUNCTIONS}" \
  "${FASH_R3_REAL_GENOTYPE_HELPER}" \
  "${FASH_R3_DIGEST_HELPER}" \
  "${FASH_R3_GENOTYPE_CACHE}" \
  "${FASH_R3_FORMAL_DIR}/complete.flag" \
  "${FASH_R3_IDEAL_DIR}/summary/replicate_input_digests.csv"; do
  if [[ ! -f "${required}" ]]; then
    echo "Required input is missing: ${required}" >&2
    exit 1
  fi
done
if [[ -e "${FASH_R3_TWO_STAGE_RESULT_PARENT}/${RESULT_ID}" ]]; then
  echo "Refusing to overwrite completed two-stage results." >&2
  exit 1
fi

module load R/4.4.1
cd "${WORKSPACE}"
Rscript --vanilla "${R_SCRIPT}" --mode preflight --project-root "${WORKSPACE}"
Rscript --vanilla "${R_SCRIPT}" --mode run --project-root "${WORKSPACE}"

if [[ ! -f "${FASH_R3_TWO_STAGE_RESULT_PARENT}/${RESULT_ID}/complete.flag" ]]; then
  echo "The runner returned without a completed ten-replicate cache." >&2
  exit 1
fi
echo "Completed internal R3 two-stage results: ${FASH_R3_TWO_STAGE_RESULT_PARENT}/${RESULT_ID}"
