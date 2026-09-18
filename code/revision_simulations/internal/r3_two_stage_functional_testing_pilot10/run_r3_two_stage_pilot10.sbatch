#!/bin/bash

#SBATCH --job-name=fash_r3_pilot10
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
readonly EXPERIMENT="${WORKSPACE}/code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10"
readonly BASE_EXPERIMENT="${WORKSPACE}/code/revision_simulations/internal/r3_two_stage_functional_testing"
readonly GENOTYPE_BUILDER="${WORKSPACE}/scripts/build_r3_two_stage_pilot10_genotypes.R"
readonly RUNNER="${WORKSPACE}/scripts/run_r3_two_stage_pilot10.R"
readonly RESULT_ID="r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_qdyn005_qfunc005_fashr0143_pilot10"
readonly GENOTYPE_ID="real_genotype_one_per_gene_J6362_pilot10"

export FASH_R3_PILOT10_CODE="${EXPERIMENT}"
export FASH_R3_TWO_STAGE_BASE_CODE="${BASE_EXPERIMENT}"
export FASH_R3_SIMULATION_FUNCTIONS="${WORKSPACE}/code/revision_simulations/r3_full_universe_functional_fashr0143/r3_full_universe_functional_simulation_functions.R"
export FASH_R3_REAL_GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
export FASH_R3_DIGEST_HELPER="${WORKSPACE}/code/revision_simulations/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R"
export FASH_R3_BASE_GENOTYPE_CACHE="${WORKSPACE}/inputs/r1_r2_fashr0143/genotype_samples.rds"
export FASH_R3_PILOT10_PAIR_SOURCE="${WORKSPACE}/inputs/r3_two_stage_pilot10/pair_key_source.rds"
export FASH_R3_PILOT10_GENOTYPE_DIR="${WORKSPACE}/inputs/r3_two_stage_pilot10/${GENOTYPE_ID}"
export FASH_R3_VCF="${PROJECT_ROOT}/data/raw/YRI_genotype.vcf.gz"
export FASH_R3_FORMAL_DIR="${WORKSPACE}/results/revision_simulations/r3_real_genotype_one_per_gene_J6362_matched_functional_open_middle_3_12_center_aligned_iwp1_geometry_mixture_relative_location_clearance_full_universe_paired_posterior_fashr0143_pilot5"
export FASH_R3_TWO_STAGE_BASE_RESULT="${WORKSPACE}/results/revision_simulations/internal/r3_two_stage_functional_testing_iwp1_mixture_qdyn005_qfunc005_fashr0143_pilot5"
export FASH_R3_PILOT10_RESULT_PARENT="${WORKSPACE}/results/revision_simulations"
export FASH_R3_PILOT10_CHECKPOINT_ROOT="${PROJECT_ROOT}/full_results/revision_simulations/${RESULT_ID}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p \
  "${WORKSPACE}/logs/out" \
  "${WORKSPACE}/logs/err" \
  "${WORKSPACE}/inputs/r3_two_stage_pilot10"

for required in \
  "${GENOTYPE_BUILDER}" \
  "${RUNNER}" \
  "${EXPERIMENT}/pilot10_helpers.R" \
  "${EXPERIMENT}/genotype_cache_helpers.R" \
  "${BASE_EXPERIMENT}/two_stage_helpers.R" \
  "${BASE_EXPERIMENT}/reconstruct_r3.R" \
  "${FASH_R3_SIMULATION_FUNCTIONS}" \
  "${FASH_R3_REAL_GENOTYPE_HELPER}" \
  "${FASH_R3_DIGEST_HELPER}" \
  "${FASH_R3_BASE_GENOTYPE_CACHE}" \
  "${FASH_R3_PILOT10_PAIR_SOURCE}" \
  "${FASH_R3_VCF}" \
  "${FASH_R3_FORMAL_DIR}/complete.flag" \
  "${FASH_R3_TWO_STAGE_BASE_RESULT}/complete.flag"; do
  if [[ ! -f "${required}" ]]; then
    echo "Required input is missing: ${required}" >&2
    exit 1
  fi
done

if [[ -e "${FASH_R3_PILOT10_RESULT_PARENT}/${RESULT_ID}" ]]; then
  echo "Refusing to overwrite completed pilot-10 results." >&2
  exit 1
fi

module load R/4.4.1
cd "${WORKSPACE}"

Rscript --vanilla "${GENOTYPE_BUILDER}" \
  --mode preflight \
  --project-root "${WORKSPACE}" \
  --base-cache "${FASH_R3_BASE_GENOTYPE_CACHE}" \
  --pair-source "${FASH_R3_PILOT10_PAIR_SOURCE}" \
  --vcf-path "${FASH_R3_VCF}" \
  --output-dir "${FASH_R3_PILOT10_GENOTYPE_DIR}"

Rscript --vanilla "${GENOTYPE_BUILDER}" \
  --mode build \
  --project-root "${WORKSPACE}" \
  --base-cache "${FASH_R3_BASE_GENOTYPE_CACHE}" \
  --pair-source "${FASH_R3_PILOT10_PAIR_SOURCE}" \
  --vcf-path "${FASH_R3_VCF}" \
  --output-dir "${FASH_R3_PILOT10_GENOTYPE_DIR}"

if [[ ! -f "${FASH_R3_PILOT10_GENOTYPE_DIR}/complete.flag" ]]; then
  echo "The ten-seed genotype cache did not complete." >&2
  exit 1
fi

Rscript --vanilla "${RUNNER}" \
  --mode preflight \
  --project-root "${WORKSPACE}"
Rscript --vanilla "${RUNNER}" \
  --mode run \
  --project-root "${WORKSPACE}"

if [[ ! -f "${FASH_R3_PILOT10_RESULT_PARENT}/${RESULT_ID}/complete.flag" ]]; then
  echo "The runner returned without a completed 20-replicate cache." >&2
  exit 1
fi
echo "Completed formal two-stage R3 cache: ${FASH_R3_PILOT10_RESULT_PARENT}/${RESULT_ID}"
