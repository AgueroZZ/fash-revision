#!/bin/bash

#SBATCH --job-name=fash_r3_geo_0143
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
readonly R_SCRIPT="${WORKSPACE}/scripts/18_r3_prior_geometry_fashr0143.R"
readonly WRAPPER_CORE="${WORKSPACE}/scripts/18_r3_prior_geometry_core.R"
readonly PRIOR_GEOMETRY_CODE="${WORKSPACE}/code/revision_simulations/r3_prior_geometry_fashr0143"
readonly R3_DRIVER="${PRIOR_GEOMETRY_CODE}/r3_prior_geometry_driver.R"
readonly SIMULATION_FUNCTIONS="${PRIOR_GEOMETRY_CODE}/r3_prior_geometry_simulation_functions.R"
readonly REAL_GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
readonly GENOTYPE_CACHE="${WORKSPACE}/inputs/r1_r2_fashr0143/genotype_samples.rds"
readonly RESULT_PARENT="${WORKSPACE}/results/revision_simulations"
readonly RESULT_ID="r3_real_genotype_one_per_gene_J6362_matched_functional_open_middle_3_12_full_support_iwp1_geometry_mixture_relative_clearance_main_effect_fashr0143_pilot5"
readonly FINAL_DIR="${RESULT_PARENT}/${RESULT_ID}"
readonly COMPLETE_FLAG="${FINAL_DIR}/complete.flag"

mkdir -p \
  "${WORKSPACE}/logs/out" \
  "${WORKSPACE}/logs/err" \
  "${RESULT_PARENT}"

for required_file in \
  "${R_SCRIPT}" \
  "${WRAPPER_CORE}" \
  "${R3_DRIVER}" \
  "${SIMULATION_FUNCTIONS}" \
  "${REAL_GENOTYPE_HELPER}" \
  "${GENOTYPE_CACHE}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Required file is missing: ${required_file}" >&2
    exit 1
  fi
done

if [[ -e "${FINAL_DIR}" ]]; then
  echo "Refusing to overwrite completed prior-geometry R3 output: ${FINAL_DIR}" >&2
  exit 1
fi

module load R/4.4.1

echo "R3 prior-geometry entry point: ${R_SCRIPT}"
echo "R3 wrapper core: ${WRAPPER_CORE}"
echo "Middle functional: 3 < t < 12"
echo "Raised-cosine primary centers: Early 1.5; Middle 4.5-10.5; Late 13.5"
echo "Every primary raised-cosine support is inside the observed domain and target region."
echo "Temporal truth-category probabilities: Early 0.29; Middle 0.42; Late 0.29"
echo "Balanced-category historical results remain the prior-shift stress test."
echo "Formal genotype cache: ${GENOTYPE_CACHE}"
echo "Slurm CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Final output: ${FINAL_DIR}"

echo "Running the prior-geometry R3 package, source, and input preflight."
FASH_WORKSPACE="${WORKSPACE}" \
  FASH_R3_RESULT_PARENT="${RESULT_PARENT}" \
  FASH_R3_DRIVER="${R3_DRIVER}" \
  FASH_R3_SIMULATION_FUNCTIONS="${SIMULATION_FUNCTIONS}" \
  Rscript --vanilla "${R_SCRIPT}" --preflight-only true

echo "Starting the resumable formal prior-geometry R3 refresh."
FASH_WORKSPACE="${WORKSPACE}" \
  FASH_R3_RESULT_PARENT="${RESULT_PARENT}" \
  FASH_R3_DRIVER="${R3_DRIVER}" \
  FASH_R3_SIMULATION_FUNCTIONS="${SIMULATION_FUNCTIONS}" \
  Rscript --vanilla "${R_SCRIPT}"

if [[ ! -f "${COMPLETE_FLAG}" ]]; then
  echo "The prior-geometry R3 wrapper returned without creating ${COMPLETE_FLAG}." >&2
  exit 1
fi

echo "R3 prior-geometry fashr 0.1.43 refresh completed: ${FINAL_DIR}"
