# Benchmarking of Virus Haplotype Reconstruction Tools
## Step 1 Generate varients sequence
Download `simulator_gui` to generate varients from your uploaded fasta file.

**Input:**

1. **Parent FASTA file** (only the first sequence in the file will be used)

    sample sequence: KC503937.fasta
      
2. **Output directory** (*existing files will be deleted)

      example of output: parent.fasta, variant1.fasta, variant2.fasta...

3. **Mutation rate (0~1)**

4. **Transition rate (0~1)**

5. **The number of variants**

6. **Choose the Mode**

    Seperate: generate X number of mutated sequences independently from the parent.
    
    Linked:Build up mutations
    
        Variant 1 – mutate from parent
        Variant 2 – mutate from varient 1
        Variant 3 – mutate from varient 2
        Variant 4 – mutate from varient 3

## Step 2 ART Illumina Paired-end sequencing simulation
1. Download ART-illumina from https://www.niehs.nih.gov/research/resources/software/biostatistics/art

    * add art_illumina to your path **(soft link)**

2. Download `art_illumina_simulator.sh`

3. Give the permission to execute

    `chmod +x art_illumina_simulator.sh`

3. Run in the terminal

    usage: 
    ```
    your_path/art_illumina_simulator.sh \
      <input_fasta_dir> \          # Directory from Step 1 output
      <overall_coverage> \         # Total sequencing depth
      "<fraction_of_each_genome>"  # Space-delimited fractions (sum≈1)
    ```
    
    example: `./art_illumina_simulator.sh ./variants_output 10000 "0.72 0.22 0.05 0.01 0.001"`

4. Output

    - simulated_R1.fq
    - simulated_R2.fq
    - aligned_simulated_reads.sam
    - sorted_simulated_reads.bam
    - sorted_simulated_reads.bam.bai

## Step 3 Test Haplotype Reconstruction Tools
### 1. PredictHaplo

```
gtime -v predicthaplo \
  --sam aligned_simulated_reads.sam \      
  --reference parent.fasta \
  --prefix output \              
  --visualization_level 0 \
  --have_true_haplotypes 0 \
  --true_haplotypes dummy.fasta \
  --do_local_Analysis 0 # run 1 first
```

output:
 - haplotypes ?
 - User time (seconds): 2177.60
    
    System time (seconds): 2.48

    Elapsed (wall clock) time (h:mm:ss or m:ss): 36:22.81
 - 	Percent of CPU this job got: 99%

	Maximum resident set size (kbytes): 467848
### 2. cliqueSNV 
input: bam file
```
gtime -v cliqueSNV \
  -m snv-illumina \
  -in simulate_reads/sorted_simulated_reads.bam \
  -ref variants_output/parent.fasta \
  -outdir clique_output \
  -threads 4
```
output:
 - 6 haplotypes (0.6392, 0.0908, 0.0728, 0.0681, 0.0655, 0.0635)
 - User time (seconds): 75.74
    
    System time (seconds): 3.01

    Elapsed (wall clock) time (h:mm:ss or m:ss): 0:29.39
- 
	Percent of CPU this job got: 267%

	Maximum resident set size (kbytes): 1108976

### 3. QuRe
input: only support FASTA file, cannot do paired-end

1.  use `seqtk` convert paired-end FASTQ to FASTA file
`cat simulated_R1.fq simulated_R2.fq | seqtk seq -a - > reads_combined.fasta`
2. create a QuRE output directory
3. Run QuRe
`gtime -v java -cp path/QuRe QuRe path/reads_combined.fasta path/parent.fasta`
### 4. SAVAGE
python 2.6+ but not python 3.x
### 5. TenSQR
1. `git clone https://github.com/SoYeonA/TenSQR.git`
2. Enter TenSQR directory run `make`
3. Check config file format (configure to your setting)