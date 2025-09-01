# load RegressHaplo
library("RegressHaplo")

# BAM file containing the dataset
args <- commandArgs(trailingOnly = TRUE)
bam_file <- args[1]
out_dir <- args[2]

# run the RegressHaplo pipeline
full_pipeline(bam_file, out_dir, start_pos=NULL, end_pos=NULL, num_trials=700)
