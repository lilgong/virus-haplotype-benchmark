#!/bin/bash
# ----input files------
zone="fmdv_mixed_simulate"
base="/home/project"
input_fq_r1="${base}/data/${zone}/simulated_R1.fastq"
input_fq_r2="${base}/data/${zone}/simulated_R2.fastq"
ref="${base}/data/${zone}/reference.fasta"

cd "${base}/data/${zone}"
echo "Data pre-processing..."
# --- use bwa-mem2 align the simulated reads with the reference reads
bwa-mem2 index "${ref}"
# align and sort to BAM file
bwa-mem2 mem "${ref}" "${input_fq_r1}" "${input_fq_r2}" > aligned_simulated_reads.sam
samtools sort -o sorted_simulated_reads.bam aligned_simulated_reads.sam
# generate index for the BAM file
samtools index sorted_simulated_reads.bam 

cat "${input_fq_r1}" "${input_fq_r2}" | seqtk seq -a - > reads_combined.fasta

echo "Done"
