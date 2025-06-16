#!/bin/bash
# ------parameters----------
input_fasta_dir=$1
total_coverage=$2
f_parameters=$3
art_illumina_path="/Users/nyl/Downloads/art_bin_MountRainier/art_illumina" ###change to your art_illumina path
read_length=250
output_dir="./simulate_reads"

mkdir -p ${output_dir}
cd ${output_dir}

# ------get the input f parameters and check if it is valid------
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
# ------check if the number of f parameters matches the file number------
fasta_files=(${input_fasta_dir}/*.fasta)
variant_count=${#fasta_files[@]}
fraction_count=${#f_parameters[@]}
if [ "$variant_count" -ne "$fraction_count" ]
then
    echo "the number of f parameters do not match the file number"
    exit 1
fi

# ------start art illumina------
i=0
for fasta_file in ${fasta_files[@]}
do
    frac=${f_parameters[$i]}
    cov=$(echo "$total_coverage * $frac" | bc)
    cov_int=$(printf "%.0f" "$cov")

    base_name=$(basename ${fasta_file} .fasta)
    base_name=${base_name%.fa}
    output_name="paired_${base_name}_"

    echo "start simulate ${fasta_file}，coverage=${cov_int}"

    ${art_illumina_path} \
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

echo "All done"
		