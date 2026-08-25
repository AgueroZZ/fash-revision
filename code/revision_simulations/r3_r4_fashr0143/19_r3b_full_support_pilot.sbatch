#!/bin/bash

#SBATCH --job-name=fash_r3b_fs_pilot
#SBATCH --account=pi-mstephens
#SBATCH --partition=caslake
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --output=/project/mstephens/ziangzhang/fash/workspace/logs/out/%x_%j.out
#SBATCH --error=/project/mstephens/ziangzhang/fash/workspace/logs/err/%x_%j.err

set -euo pipefail

readonly WORKSPACE="/project/mstephens/ziangzhang/fash/workspace"
readonly CODE_DIR="${WORKSPACE}/code/revision_simulations/r3_prior_geometry_fashr0143"
readonly DRIVER="${CODE_DIR}/r3_prior_geometry_driver.R"
readonly SIMULATION_FUNCTIONS="${CODE_DIR}/r3_prior_geometry_simulation_functions.R"
readonly VALIDATOR="${WORKSPACE}/scripts/19_r3b_full_support_pilot_validate.R"
readonly GENOTYPE_HELPER="${WORKSPACE}/code/revision_simulations/shared/real_genotype_one_per_gene.R"
readonly GENOTYPE_CACHE="${WORKSPACE}/inputs/r1_r2_fashr0143/genotype_samples.rds"
readonly RESULT_ID="r3b_real_genotype_one_per_gene_J6362_open_middle_3_12_full_support_iwp1_geometry_mixture_fashr0143_seed12345_pilot"
readonly RESULT_DIR="${WORKSPACE}/results/revision_simulations/${RESULT_ID}"

mkdir -p "${WORKSPACE}/logs/out" "${WORKSPACE}/logs/err"

for required_file in \
  "${DRIVER}" \
  "${SIMULATION_FUNCTIONS}" \
  "${VALIDATOR}" \
  "${GENOTYPE_HELPER}" \
  "${GENOTYPE_CACHE}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "Required pilot file is missing: ${required_file}" >&2
    exit 1
  fi
done

echo "32a0c6e3fb99e355d6c0ab0d68acbf7818b01d12f4516236addbfda2b33d48f9  ${DRIVER}" | sha256sum --check --status
echo "3c5a97612591817a7edbbd3c331b380fae1d3a2fb468b4a78c0aed13ee3b8db9  ${SIMULATION_FUNCTIONS}" | sha256sum --check --status

module load R/4.4.1

echo "Starting the single-seed R3B full-support design pilot."
echo "This pilot is not a formal calibration result."
echo "Result directory: ${RESULT_DIR}"

FASH_R3_SIMULATION_FUNCTIONS="${SIMULATION_FUNCTIONS}" \
FASH_R3_REAL_GENOTYPE_HELPER="${GENOTYPE_HELPER}" \
Rscript --vanilla "${DRIVER}" \
  --J 6362 \
  --n-donors 19 \
  --n-covariates 5 \
  --noise-sd 1 \
  --dynamic-main-effect-sd 1 \
  --truth-mechanisms raised_cosine \
  --seed-list 12345 \
  --genotype-cache "${GENOTYPE_CACHE}" \
  --output-id "${RESULT_ID}" \
  --output-dir "${RESULT_DIR}" \
  --num-cores "${SLURM_CPUS_PER_TASK}" \
  --expected-fashr-version 0.1.43 \
  --expected-fashr-remote-sha bf223df75da6e41ae48607a56b4cd12d7c3b24e7 \
  --expected-genotype-content-md5 526a7318aa2af901e09252f5a6ca3c46 \
  --middle-window 3,12 \
  --middle-boundary open \
  --temporal-category-probs 0.29,0.42,0.29

Rscript --vanilla "${VALIDATOR}" --result-dir "${RESULT_DIR}"

if [[ ! -f "${RESULT_DIR}/pilot_complete.flag" ]]; then
  echo "Pilot validation did not create pilot_complete.flag." >&2
  exit 1
fi

echo "R3B full-support design pilot completed."
