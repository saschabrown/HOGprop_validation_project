#!/usr/bin/env bash
set -euo pipefail

: "${WP1_ROOT:?Set WP1_ROOT to the absolute workpackage1 path}"
if [[ "$WP1_ROOT" != /* || $# != 2 ]]; then
    echo "Usage: run_fastoma.sh PREPARED_RUN_ID NEW_RUN_ID" >&2
    exit 2
fi
for run_id in "$@"; do
    if [[ ! "$run_id" =~ ^[a-zA-Z0-9][a-zA-Z0-9._-]*$ ]]; then
        echo "Invalid run ID: $run_id" >&2
        exit 2
    fi
done

repo="${WP1_ROOT}/HOGprop_validation_project"
stage="${repo}/analysis/01_fastoma"
input="${stage}/data/$1/fastoma_input"
output="${stage}/data/$2"
reports="${stage}/results/$2"
scratch="${WP1_ROOT}/working_directory/01_fastoma/$2"
workflow="${repo}/software/fastoma/FastOMA.nf"
config="${repo}/config/fastoma_genomedk_slurm.config"

for file in "$workflow" "$config" "$input/species_tree.nwk" "$input/LUCA.h5"; do
    if [[ ! -s "$file" ]]; then
        echo "Required file is missing or empty: $file" >&2
        exit 1
    fi
done
if [[ ! -d "$input/proteome" ]]; then
    echo "Required proteome directory is missing: $input/proteome" >&2
    exit 1
fi
command -v nextflow >/dev/null
for path in "$output" "$reports" "$scratch"; do
    if [[ -e "$path" || -L "$path" ]]; then
        echo "Run path already exists: $path" >&2
        exit 1
    fi
done
mkdir -p "${stage}/data" "${stage}/results" "${WP1_ROOT}/working_directory/01_fastoma"
mkdir "$output" "$reports" "$scratch"
cd "$scratch"
nextflow run "$workflow" \
    -profile slurm,large \
    -c "$config" \
    --input "$input" \
    --output_folder "$output" \
    --omamer_db "$input/LUCA.h5" \
    -work-dir "$scratch/work" \
    -ansi-log false \
    -with-trace "$reports/nextflow_trace.tsv" \
    -with-report "$reports/nextflow_report.html" \
    -with-timeline "$reports/nextflow_timeline.html"
