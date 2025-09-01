#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# zone="seperate_mr0.01"
zone="fmdv_mixed_simulate"
base="/home/2899321l/project"
input_file="${base}/data/${zone}/aligned_simulated_reads.sam"
# ref_file="${base}/data/${zone}/parent.fasta"
ref_file="${base}/data/${zone}/reference.fasta"
output_dir="${base}/results/tensqr/${zone}"
log_dir="${base}/logs/tensqr/${zone}"
log_file="${log_dir}/run.log"
time_log="${log_dir}/time.log"
config_file="${base}/data/${zone}/config"

# ===== parameters =====
snv_thres="0.01"
start_pos=1
# stop_pos=8239
stop_pos=8200
min_map_qual=40
min_read_len=100
max_insert=500
seq_err=0.2
mec_thresh=0.0312
init_size=5

# ===== create directories =====
mkdir -p "${output_dir}"
mkdir -p "${log_dir}"

# ===== generate config file =====
cat > "$config_file" <<EOF
filename of reference sequence (FASTA) : ${ref_file}
filname of the aligned reads (sam format) : ${input_file}
SNV_thres : ${snv_thres}
reconstruction_start : ${start_pos}
reconstruction_stop: ${stop_pos}
min_mapping_qual : ${min_map_qual}
min_read_length : ${min_read_len}
max_insert_length : ${max_insert}
characteristic zone name : ${zone}
seq_err (assumed sequencing error rate(%)) : ${seq_err}
MEC improvement threshold : ${mec_thresh}
initial population size : ${init_size}
EOF

echo "${config_file} generated"

# ===== run TenSQR=====
cd "${output_dir}"

echo "Start Running TenSQR..." | tee "$log_file"

/usr/bin/time -v tensqr "${config_file}" > "${log_file}" 2> "${time_log}"

# ===== convert to fasta format ======
output_seq_txt="${output_dir}/${zone}_ViralSeq.txt"
output_seq_fasta="${output_dir}/${zone}_ViralSeq.fasta"
sed 's/Viral Quasispecies/>ViralQuasispecies/g' "${output_seq_txt}" > "${output_seq_fasta}"
# ==== replace * to N ====
sed -i 's/\*/N/g' "${output_seq_fasta}"
# ==== delete space in the header, replace :- to _ ====
sed -i '/^>/ {
    s/^>\(.*\)/\1/;
    s/ //g;
    s/[^A-Za-z0-9_]/_/g;
    s/^/>/
}' "${output_seq_fasta}"
echo "==== Convert Result to Fasta ===="
# ==== Done ====
echo "==== TenSQR finished at $(date) ===="