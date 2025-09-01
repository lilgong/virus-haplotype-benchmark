#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

zone="fmdv_mixed_simulate"
base="/home/project"
input_bam="${base}/data/${zone}/sorted_simulated_reads.bam"
ref="${base}/data/${zone}/reference.fasta"
output_dir="${base}/results/cliquesnv/${zone}"
log_dir="${base}/logs/cliquesnv/${zone}"
log_file="${log_dir}/run.log"
time_log="${log_dir}/time.log"

mkdir -p "${output_dir}"
mkdir -p "${log_dir}"

echo "Running CliqueSNV..."

/usr/bin/time -v cliquesnv \
  -m snv-illumina \
  -in "${input_bam}" \
  -ref "${ref}" \
  -outDir "${output_dir}" \
  -t 5 \
  -tf 0.001 \
  -fdf extended \
  -threads 4 \
  > "${log_file}" 2> "${time_log}"

echo "==== CliqueSNV finished at $(date) ===="