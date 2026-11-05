#!/bin/bash -l

#SBATCH --account=evomics_projects
#SBATCH --job-name=fastoma_large_dataset
#SBATCH --output="/home/saschabrown/evomics_projects/saschabrown/HOGprop_validation_project/logs/slurm/fastoma_large_%j.out"
#SBATCH --error="/home/saschabrown/evomics_projects/saschabrown/HOGprop_validation_project/logs/slurm/fastoma_large_%j.err"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=72G
# Launcher resources are separate from the 256 GB maximum for pipeline tasks.
#SBATCH --time=72:00:00
#SBATCH --export=NONE

source /home/saschabrown/miniforge3/etc/profile.d/conda.sh || exit 1
conda activate fastoma_omamer212 || exit 1
export SLURM_EXPORT_ENV=ALL

set -e
set -o pipefail

RUN_DIR="/home/saschabrown/evomics_projects/saschabrown/HOGprop_validation_project"
RUN_TAG="$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID}"

cd "$RUN_DIR"

# Create the report directory if needed.
mkdir -p "${RUN_DIR}/reports/nextflow"

# Confirm the configuration and required inputs exist.
# The Slurm log directory must already exist before sbatch is called.
test -r "${RUN_DIR}/config/fastoma_genomedk_slurm.config"
test -d "${RUN_DIR}/fastoma_input/proteome"
test -s "${RUN_DIR}/fastoma_input/species_tree.nwk"
test -s "${RUN_DIR}/fastoma_input/LUCA.h5"
test -s "${RUN_DIR}/FastOMA/FastOMA.nf"

nextflow run "${RUN_DIR}/FastOMA/FastOMA.nf" \
    -profile slurm,large \
    -c "${RUN_DIR}/config/fastoma_genomedk_slurm.config" \
    --input "${RUN_DIR}/fastoma_input" \
    --output_folder "${RUN_DIR}/fastoma_output/${RUN_TAG}" \
    --omamer_db "${RUN_DIR}/fastoma_input/LUCA.h5" \
    -work-dir "${RUN_DIR}/fastoma_work/${RUN_TAG}" \
    -ansi-log false \
    -with-trace "${RUN_DIR}/reports/nextflow/fastoma_trace_${SLURM_JOB_ID}.tsv" \
    -with-report "${RUN_DIR}/reports/nextflow/fastoma_report_${SLURM_JOB_ID}.html" \
    -with-timeline "${RUN_DIR}/reports/nextflow/fastoma_timeline_${SLURM_JOB_ID}.html"
