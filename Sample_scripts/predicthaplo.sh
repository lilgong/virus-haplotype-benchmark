#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

dummy='/home/project/data/dummy.fasta'

zone="fmdv_mixed_simulate"
base="/home/project"
input_file="${base}/data/${zone}/aligned_simulated_reads.sam"
ref_file="${base}/data/${zone}/reference.fasta"
output_dir="${base}/results/predicthaplo/${zone}"
log_dir="${base}/logs/predicthaplo/${zone}"
log_file="${log_dir}/run.log"
time_log="${log_dir}/time.log"


echo "==== PredictHaplo test run started at $(date) ===="

# check if input files exist
if [[ ! -f "${input_file}" ]]; then
    echo "ERROR: Input SAM file not found: ${input_file}"
    exit 1
fi

if [[ ! -f "${ref_file}" ]]; then
    cho "ERROR: Reference FASTA file not found: ${ref_file}"
    exit 1
fi

# create output and log directory for data analyzed with predicthaplo
mkdir -p "${output_dir}"
mkdir -p "${log_dir}"

# create dummy.fasta if it does not exist
if [[ ! -f "${dummy}" ]]; then
    echo -e ">dummy\nACGT" > "${dummy}"
    echo "Created dummy fasta file: ${dummy}"
fi

echo "Running PredictHaplo with parameters:"
echo "  SAM: ${input_file}"
echo "  Reference: ${ref_file}"
echo "  Output prefix: ${output_dir}"
echo "  Dummy haplotype: ${dummy}"
echo

# run PredictHaplo
/usr/bin/time -v predicthaplo \
    --sam "${input_file}" \
    --reference "${ref_file}" \
    --prefix "${output_dir}/run1" \
    --visualization_level 0 \
    --have_true_haplotypes 0 \
    --true_haplotypes "${dummy}" \
    --do_local_Analysis 1 \
    --entropy_threshold 0.3 \
    --reconstruction_start 1 \
    --local_window_size_factor 0.2 \
    --max_reads_in_window 10000 \
    > "${log_file}" 2> "${time_log}"

echo
echo "==== PredictHaplo finished at $(date) ===="
