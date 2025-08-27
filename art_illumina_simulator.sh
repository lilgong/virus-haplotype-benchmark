#!/bin/bash
# ------parameters----------
input_fasta_dir=$1
total_coverage=$2
f_parameters=$3
read_length=250
output_dir="./simulate_reads"

# ------ check if the input is valid-----------
## ----- check the input directory storing fasta files exists----------
if [ ! -d "${input_fasta_dir}" ]; then
    echo "Error: input directory '${input_fasta_dir}' does not exist."
    exit 1
fi
## Convert to absolute path
input_fasta_dir="$(cd "$1" && pwd)"

## ------get the input f parameters and check if it is valid------
IFS=' ' read -r -a f_parameters <<< ${f_parameters}
sum=0
for frac in ${f_parameters[@]}
do
	#check if frac > 0
	if ! awk "BEGIN {exit !(${frac} > 0 && ${frac} < 1)}"
	then
		echo "Invalid input, f parameters should be between 0 and 1"
		exit 1
	fi
	sum=$(echo "${sum} + ${frac}" | bc)
done
if ! awk "BEGIN {exit !(${sum} >= 0.99 && ${sum} <= 1.01)}"
then
    echo "The sum of f parameters should be 1"
    exit 1
fi		

## ------check if the number of f parameters matches the file number------
fasta_files=(${input_fasta_dir}/*.fasta)
variant_count=${#fasta_files[@]}
fraction_count=${#f_parameters[@]}
if [ "$variant_count" -ne "$fraction_count" ]
then
    echo "the number of f parameters do not match the file number"
    exit 1
fi

# ------ create output directory -------
mkdir -p ${output_dir}
cd ${output_dir}
# ------start art illumina------
i=0
for fasta_file in ${fasta_files[@]}
do
    frac=${f_parameters[$i]}
    cov=$(echo "$total_coverage * $frac" | bc)
    cov_int=$(printf "%.0f" "$cov")

    base_name=$(basename ${fasta_file} .fasta)
    base_name=${base_name%.fa}
    output_name="paired_${base_name}_R"

    echo "start simulate ${fasta_file}, coverage=${cov_int}"

    art_illumina\
	-sam \
        -i ${fasta_file} \
	-p \
        -l ${read_length} \
	-ss MSv3 \
        -f ${cov_int} \
	-m 500 \
	-s 20 \
        -o ${output_name}

    ((i++))
done

cat paired_*_R1.fq > simulated_R1.fastq
cat paired_*_R2.fq > simulated_R2.fastq
# remove the original paired reads (keep them if want)
read -p "Only keep combined FASTQ files? (y/n): " confirm && [[ $confirm == [Yy] ]] && rm paired*

# --- use bwa-mem2 align the simulated reads with the reference (parent) reads
bwa-mem2 index ${input_fasta_dir}/parent.fasta
# align and sort to BAM file
bwa-mem2 mem ${input_fasta_dir}/parent.fasta simulated_R1.fastq simulated_R2.fastq > aligned_simulated_reads.sam
samtools sort -o sorted_simulated_reads.bam aligned_simulated_reads.sam
# generate index for the BAM file
samtools index sorted_simulated_reads.bam

# convert to fasta file for some tools
cat simulated_R1.fastq simulated_R2.fastq | seqtk seq -a - > reads_combined.fasta

echo "All done"
