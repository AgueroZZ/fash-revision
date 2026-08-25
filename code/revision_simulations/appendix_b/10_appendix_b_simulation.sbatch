#!/bin/bash

#SBATCH --job-name=fash_appendix_b_0143
#SBATCH --account=pi-mstephens
#SBATCH --partition=caslake
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/project/mstephens/ziangzhang/fash/workspace/logs/out/%x_%j.out
#SBATCH --error=/project/mstephens/ziangzhang/fash/workspace/logs/err/%x_%j.err

set -euo pipefail

readonly PROJECT_ROOT="/project/mstephens/ziangzhang/fash"
readonly WORKSPACE="${PROJECT_ROOT}/workspace"
readonly R_SCRIPT="${WORKSPACE}/scripts/10_appendix_b_simulation.R"
readonly HELPER_SCRIPT="${WORKSPACE}/code/appendix_b_helpers.R"
readonly RESULT_PARENT="${WORKSPACE}/results/revision_simulations"
readonly FINAL_DIR="${RESULT_PARENT}/appendix_b_fashr0143"
readonly COMPLETE_FLAG="${FINAL_DIR}/complete.flag"

mkdir -p \
  "${WORKSPACE}/logs/out" \
  "${WORKSPACE}/logs/err" \
  "${RESULT_PARENT}"

for required_file in "${R_SCRIPT}" "${HELPER_SCRIPT}"; do
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

echo "R script: ${R_SCRIPT}"
echo "Helper script: ${HELPER_SCRIPT}"
echo "Slurm CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Final output: ${FINAL_DIR}"

Rscript --vanilla -e '
expected_version <- "0.1.43"
expected_sha <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is not available in this R environment.")
}
description <- packageDescription("fashr")
actual_version <- as.character(packageVersion("fashr"))
actual_sha <- if (is.null(description$RemoteSha)) NA_character_ else description$RemoteSha
cat("fashr version: ", actual_version, "\n", sep = "")
cat("fashr RemoteSha: ", actual_sha, "\n", sep = "")
cat("fashr path: ", find.package("fashr"), "\n", sep = "")
cat("R library paths:\n")
cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")
if (!identical(actual_version, expected_version)) {
  stop("Expected fashr ", expected_version, "; found ", actual_version, ".")
}
if (!identical(actual_sha, expected_sha)) {
  stop("Expected fashr RemoteSha ", expected_sha, "; found ", actual_sha, ".")
}
'

echo "Starting the resumable Appendix B simulation."
FASH_WORKSPACE="${WORKSPACE}" Rscript --vanilla "${R_SCRIPT}"

if [[ ! -f "${COMPLETE_FLAG}" ]]; then
  echo "The runner returned without creating ${COMPLETE_FLAG}." >&2
  exit 1
fi

echo "Appendix B completed: ${FINAL_DIR}"
