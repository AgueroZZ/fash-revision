#!/bin/bash
#SBATCH --job-name=fash_r3_compare
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
readonly RESULT_ID="r3_stage1_iwp1_linear_direct_paired_raw_bf_fashr0143_pilot10"
readonly EXPERIMENT="${WORKSPACE}/code/revision_simulations/internal/r3_stage1_method_comparison"
readonly RUNNER="${WORKSPACE}/scripts/run_r3_stage1_method_comparison.R"

export FASH_R3_COMPARE_BASE="${WORKSPACE}/results/revision_simulations/r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_qdyn005_qfunc005_fashr0143_pilot10"
export FASH_R3_COMPARE_GENOTYPES="${WORKSPACE}/inputs/r3_two_stage_pilot10/real_genotype_one_per_gene_J6362_pilot10/genotype_samples.rds"
export FASH_R3_COMPARE_SIMULATION="${WORKSPACE}/code/revision_simulations/r3_full_universe_functional_fashr0143/r3_full_universe_functional_simulation_functions.R"
export FASH_R3_COMPARE_GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
export FASH_R3_COMPARE_DIGEST_HELPER="${WORKSPACE}/code/revision_simulations/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R"
export FASH_R3_COMPARE_RESULT_PARENT="${WORKSPACE}/results/revision_simulations"
export FASH_R3_COMPARE_CHECKPOINT_ROOT="${PROJECT_ROOT}/full_results/revision_simulations/${RESULT_ID}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p "${WORKSPACE}/logs/out" "${WORKSPACE}/logs/err"
for required in "${RUNNER}" "${EXPERIMENT}/comparison_helpers.R" \
  "${FASH_R3_COMPARE_BASE}/complete.flag" "${FASH_R3_COMPARE_GENOTYPES}" \
  "${FASH_R3_COMPARE_SIMULATION}" "${FASH_R3_COMPARE_GENOTYPE_HELPER}" \
  "${FASH_R3_COMPARE_DIGEST_HELPER}" \
  "${WORKSPACE}/code/revision_simulations/r1_r2_fashr0143/source_snapshots/r1_r2_fashr0143_simulation_functions.R"; do
  if [[ ! -f "${required}" ]]; then
    echo "Required source or completed input is missing: ${required}" >&2
    exit 1
  fi
done
if [[ -e "${FASH_R3_COMPARE_RESULT_PARENT}/${RESULT_ID}" ]]; then
  echo "Refusing to overwrite completed R3 comparison results." >&2
  exit 1
fi
module load R/4.4.1
cd "${WORKSPACE}"
Rscript --vanilla "${RUNNER}" --mode preflight --project-root "${WORKSPACE}"
Rscript --vanilla "${RUNNER}" --mode run --project-root "${WORKSPACE}"
test -f "${FASH_R3_COMPARE_RESULT_PARENT}/${RESULT_ID}/complete.flag"
echo "Completed paired R3 method comparison: ${RESULT_ID}"
