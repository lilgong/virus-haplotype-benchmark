#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# === Activate regressHaplo conda environment ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate regresshaplo_env

zone="fmdv_mixed_simulate"
base="/home/project"
input_bam="${base}/data/${zone}/sorted_simulated_reads.bam"
output_dir="${base}/results/regresshaplo/${zone}"
log_dir="${base}/logs/regresshaplo/${zone}"

log="${log_dir}/run.log"
time_log="${log_dir}/time.log"

mkdir -p "${output_dir}" "${log_dir}"

echo "Running RegressHaplo..."

/usr/bin/time -v Rscript /home/project/scripts/regresshaplo.R "${input_bam}" "${output_dir}" > "${log}" 2> "${time_log}"

echo "=== RegressHaplo Finished at $(date) ==="

conda deactivate
