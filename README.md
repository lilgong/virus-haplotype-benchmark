# Benchmarking of Virus Haplotype Reconstruction Tools
## Step 1 Generate varients sequence
Download `simulator_gui` to generate varients from your uploaded fasta file.

**Input:**

1. **Parent FASTA file** (only the first sequence in the file will be used)

    sample sequence: KC503937.fasta
      
2. **Output directory** (existing files will not be deleted. recommend to make it empty to avoid overwriting)

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

    [Option] add art_illumina to the system path

2. Download `art_illumina_simulator.sh`, edit the `art_illumina_path` to your own path.

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



