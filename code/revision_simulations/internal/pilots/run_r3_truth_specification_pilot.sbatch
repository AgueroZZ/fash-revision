#!/bin/bash

#SBATCH --job-name=fash_r3_truth_pilot
#SBATCH --account=pi-mstephens
#SBATCH --partition=caslake
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --array=0-4
#SBATCH --output=/project/mstephens/ziangzhang/fash/workspace/logs/out/%x_%A_%a.out
#SBATCH --error=/project/mstephens/ziangzhang/fash/workspace/logs/err/%x_%A_%a.err

# Submit only after reviewing the J=300 exploratory summary. This job produces
# an exploratory five-seed real-genotype pilot, not a replacement formal cache.

set -euo pipefail

readonly PROJECT_ROOT="/project/mstephens/ziangzhang/fash"
readonly WORKFLOWR_ROOT="${PROJECT_ROOT}/fash-revision"
readonly RUNNER="${WORKFLOWR_ROOT}/code/revision_simulations/internal/pilots/run_r3_truth_specification_pilot.R"
readonly SUMMARY_RUNNER="${WORKFLOWR_ROOT}/code/revision_simulations/internal/pilots/summarize_r3_truth_specification_pilot.R"
readonly R_BINARY="/software/R-4.4.1-el8-x86_64+gcc-13.2.0/bin/Rscript"
readonly R_LIBRARY="/home/ziangzhang/R/x86_64-pc-linux-gnu-library/4.4"
readonly GCC_LIBRARY="/software/gcc-13.2.0-el8-x86_64/lib64"
readonly GENOTYPE_CACHE="${PROJECT_ROOT}/workspace/inputs/r1_r2_fashr0143/genotype_samples.rds"
readonly TRUTH_MECHANISM="${FASH_R3_TRUTH_MECHANISM:-r3a_constrained_global_bspline}"
readonly RESULT_PREFIX="${FASH_R3_PILOT_PREFIX:-r3a_constrained_global_bspline_J6362}"
readonly SEEDS=(12345 22345 32345 42345 52345)
readonly SEED="${SEEDS[${SLURM_ARRAY_TASK_ID}]}"

if [[ ! -x "${R_BINARY}" || ! -f "${RUNNER}" || ! -f "${GENOTYPE_CACHE}" ]]; then
  echo "Missing R binary, pilot runner, or genotype cache." >&2
  exit 1
fi
if [[ "${TRUTH_MECHANISM}" != "r3a_targeted_broad_bspline" &&
      "${TRUTH_MECHANISM}" != "r3a_constrained_global_bspline" &&
      "${TRUTH_MECHANISM}" != "r3b_interior_single_lobe" ]]; then
  echo "Unsupported FASH_R3_TRUTH_MECHANISM: ${TRUTH_MECHANISM}" >&2
  exit 1
fi

export LD_LIBRARY_PATH="${GCC_LIBRARY}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export R_LIBS_USER="${R_LIBRARY}"

cd "${WORKFLOWR_ROOT}"
"${R_BINARY}" --vanilla "${RUNNER}" \
  --truth-mechanism "${TRUTH_MECHANISM}" \
  --seed "${SEED}" \
  --J 6362 \
  --num-cores "${SLURM_CPUS_PER_TASK}" \
  --genotype-cache "${GENOTYPE_CACHE}" \
  --output-id "${RESULT_PREFIX}_seed${SEED}"

if [[ "${SLURM_ARRAY_TASK_ID}" == "0" ]]; then
  cat <<EOF
Each array task saves one seed. After all five succeed, aggregate them with:
${R_BINARY} --vanilla ${SUMMARY_RUNNER} \
  --output-prefix ${RESULT_PREFIX} \
  --summary-file output/revision_simulations/internal/${RESULT_PREFIX}_alpha005_summary.csv
EOF
fi
