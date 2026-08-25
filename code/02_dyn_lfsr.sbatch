#!/bin/bash

#SBATCH --job-name=fash-middle-lfsr
#SBATCH --account=pi-mstephens
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=02_dyn_lfsr_%j.out

set -euo pipefail

# Submit from the project root, or pass the project root as the first argument:
#   sbatch code/02_dyn_lfsr.sbatch /path/to/fash-project
#
# The optional second argument is a new output directory. Existing paths are
# rejected to prevent accidental replacement of retained manuscript results.
PROJECT_ROOT="${1:-${SLURM_SUBMIT_DIR}}"
OUTPUT_DIRECTORY="${2:-results/lfsr_open_middle_3_12_seed20260820}"

module load R/4.4.1
cd "${PROJECT_ROOT}"

echo "Project root: ${PROJECT_ROOT}"
echo "Output directory: ${OUTPUT_DIRECTORY}"
echo "Slurm workers: ${SLURM_CPUS_PER_TASK}"

Rscript code/02_dyn_lfsr.R \
  --output-dir "${OUTPUT_DIRECTORY}" \
  --categories middle \
  --grid-step 0.10 \
  --posterior-draws 3000 \
  --num-cores "${SLURM_CPUS_PER_TASK}" \
  --seed 20260820 \
  --alpha 0.05 \
  --expected-pairs 9205
