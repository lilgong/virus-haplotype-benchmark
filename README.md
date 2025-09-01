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

    * add art_illumina to your path

2. Packages: bwa-mem2, samtools, seqtk
3. Download `art_illumina_simulator.sh` and Give the permission to execute `chmod +x art_illumina_simulator.sh`
4. Run in the terminal

    usage: 
    ```
    your_path/art_illumina_simulator.sh \
      <input_fasta_dir> \          # Directory from Step 1 output
      <overall_coverage> \         # Total sequencing depth
      "<fraction_of_each_genome>"  # Space-delimited fractions (sum≈1)
    ```
    
    example: `./art_illumina_simulator.sh ./variants_output 10000 "0.72 0.22 0.05 0.01 0.001"`

5. Output

    - simulated_R1.fastq
    - simulated_R2.fastq
    - aligned_simulated_reads.sam
    - sorted_simulated_reads.bam
    - sorted_simulated_reads.bam.bai
    - reads_combined.fasta

## Step 3 Test Haplotype Reconstruction Tools

### Tool List (detailed list is in Summary_of_Haplotype_Tools.xlsx)

| Tool | Source Link |
|------|-------------|
| SAVAGE | https://github.com/HaploConduct/HaploConduct |
| aBayesQR | https://github.com/SoYeonA/aBayesQR |
| CliqueSNV | https://github.com/vtsyvina/CliqueSNV |
| HaploClique | https://github.com/armintoepfer/haploclique |
| PredictHaplo | https://github.com/cbg-ethz/PredictHaplo |
| QuasiRecomb | https://github.com/cbg-ethz/QuasiRecomb |
| QuRe | https://sourceforge.net/projects/qure/ |
| RegressHaplo | https://github.com/SLeviyang/RegressHaplo |
| ShoRAH | https://github.com/cbg-ethz/shorah |
| TenSQR | https://github.com/SoYeonA/TenSQR |
| ViQuaS | https://sourceforge.net/projects/viquas/ |
| BIOA (with VirA) | https://alan.cs.gsu.edu/vira/#welcome-to-vira-and-bioa-s-documentation |
| HaploDMF | https://github.com/dhcai21/HaploDMF |
| HaROLD | https://github.com/RichardAGoldstein/HaROLD |
| NeurHap | https://github.com/xuehansheng/NeurHap |
| paphpipe | https://github.com/gwcbi/haphpipe |
| Qcolors | not found |
| QSdpR | https://sourceforge.net/projects/qsdpr/ |
| ShotMCF | https://alan.cs.gsu.edu/NGS/?q=content/shotmcf |
| VILOCA | https://github.com/cbg-ethz/VILOCA |
| ViSpA | https://alan.cs.gsu.edu/NGS/?q=content/vispa |
| V-Phaser2 | https://github.com/broadinstitute/v-phaser2 |
| 2SNV | https://alan.cs.gsu.edu/NGS/?q=content/2snv |
| V-Phaser + V-Profiler1 | https://www.broadinstitute.org/viral-genomics/viral-genomics-analysis-software-registration |
| VGA | http://genetics.cs.ucla.edu/vga/ |
| StrainLine | https://github.com/HaploKit/Strainline |
| ViQUF | https://github.com/borjaf696/ViQUF |
| VStrains | https://github.com/metagentools/VStrains |
| Mutant-Bin | inaccessible |
| PEHaplo | https://github.com/chjiao/PEHaplo |
| VG-Flow | https://bitbucket.org/jbaaijens/vg-flow/src/master/ |
| viaDBG | https://github.com/borjaf696/viaDBG |
| ViPRA-Haplo (MLEHaplo) | https://github.com/raunaq-m/MLEHaplo |
| Virus-VG | https://bitbucket.org/jbaaijens/virus-vg/src/master/ |

---

### Project Structure
```
project/
├── data/                  # Input files
│   ├── zone1/
│   │   ├── simulated_R1.fastq
│   │   ├── simulated_R2.fastq
│   │   ├── aligned_simulated_reads.sam
│   │   ├── sorted_simulated_reads.bam
│   │   ├── sorted_simulated_reads.bam.bai
│   │   ├── reads_combined.fasta
│   │   └── reference.fasta
│   ├── zone2/
│   │   └── ...
│   └── ...
├── logs/                  # Log files for each tool run
│   ├── zone1/
│   │   ├── time.log
│   │   └── run.log
│   ├── zone2/
│   │   ├── time.log
│   │   └── run.log
│   └── ...
├── results/               # outputs
│   ├── tool1/
│   ├── tool2/
│   └── ...
├── scripts/               # Run scripts (bash, python, etc.)
└── tools/                 # Installed tools (added to PATH/bin)
```
---
### Execution
- Each tool is run in its own conda environment. 
- Results and logs are stored in tool-specific folders.  

Below is a template bash script to standardize execution:

```

#!/bin/bash

# === Activate conda environment ===
source ~/miniconda3/etc/profile.d/conda.sh
conda activate haploclique_env   # change to your env name

# Define paths
zone="zone_name"        # change to your zone name
tool=“tool_name”        # change to your tool name
base="/home/project"
input_bam="${base}/data/${zone}/sorted_simulated_reads.bam"
output_dir="${base}/results/${tool}/${zone}"
log_dir="${base}/logs/${tool}/${zone}"

log="${log_dir}/run.log"
time_log="${log_dir}/time.log"

# Create folders
mkdir -p "${output_dir}" "${log_dir}"

# Run tool
echo "Running Tool_Name at $(date)..."
/usr/bin/time -v command        # give the command based on your settings
    > "${log}" 2> "${time_log}"

echo "=== Finished at $(date) ==="

conda deactivate

```

---
Note that you may need to format the output file. (e.g., `sed`)

## Step 4 Evaluation
- Merge ground truth and predicted haplotypes into a single FASTA file.

  `cat true.fasta predict.fasta > all.fasta`
- Perform multiple sequence alignment using MAFFT.

  `mafft --auto all.fasta > all_aln.fasta`
- Construct a phylogenetic tree with IQ-TREE.

  `iqtree -s all_aln.fasta -m GTR+G -nt AUTO -pre all_sequences`
- Compute the distance (e.g., using Clustal Omega):

  `clustalo -i all.fasta -o all_aln.fasta --outfmt=clu --force`
- Alternatively, sequence distance, weighted UniFrac distance, as well as precision and recall can be computed using the provided Python scripts (distance.py & precision.py)