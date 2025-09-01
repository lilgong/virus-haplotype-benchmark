#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# The read file and the reference genome file must be in FASTA format (nucleotides).
zone="fmdv_mixed_simulate"
base="/home/project"
input_file="${base}/data/${zone}/reads_combined.fasta"
ref_file="${base}/data/${zone}/reference.fasta"
output_dir="${base}/results/qure/${zone}"
log_dir="${base}/logs/qure/${zone}"
log_file="${log_dir}/run.log"
time_log="${log_dir}/time.log"

mkdir -p "${output_dir}" "${log_dir}"

# ==== parameters ====
# homo_err=0.01
# nonhomo_err=0.005
# iterations=3000


# ==== Run QuRe ====
reads="${output_dir}/reads.fasta"
ref="${output_dir}/ref.fasta"
ln -sf "${input_file}" "${reads}"
ln -sf "${ref_file}" "${ref}"

echo "Running QuRe..."
/usr/bin/time -v java -Xmx3G -cp /home/2899321l/project/tools/QuRe QuRe \
    "${reads}" "${ref}" \
    > "${log_file}" 2> "${time_log}"
    # "${homo_err}" "${nonhomo_err}" "${iterations}" \
    # > "${log_file}" 2> "${time_log}"

# ==== convert to FASTA ====
cp "${output_dir}"/*_reconstructedVariants.txt "${output_dir}/reconstructedVariants.fasta"
# ==== Done ====
echo "==== QuRe finished at $(date) ===="