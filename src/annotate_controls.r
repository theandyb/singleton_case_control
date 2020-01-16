# Script to annotate our sampled controls with genomic feature information

# pass in chromosome number as input plz
args = commandArgs(trailingOnly=TRUE)
chr = args[1]

print(paste0("Now annotating controls from chromosome ", chr))

infile <- paste0("/net/snowwhite/home/beckandy/research/singleton_case_control/data/chr", chr, ".csv")
outfile <- paste0("/net/snowwhite/home/beckandy/research/singleton_case_control/data/chr", chr, "_annotated.csv")

library(tidyverse)

source("/net/snowwhite/home/beckandy/research/singleton_case_control/src/smaug_funcs.r")

col_names_controls <- c("mut_type", "chr", "seq" , "pos", "ref")
df <- read_csv(infile, col_names = col_names_controls)

# paths to reference files
ref_dir <- "/net/snowwhite/home/beckandy/research/smaug-redux/reference_data/"
exon_file <- paste0(ref_dir, "GRCh37_RefSeq_sorted.bed")
cpg_island_file <- paste0(ref_dir, "cpg_islands_sorted.bed")
rcr_file <- paste0(ref_dir, "recomb_rate.bed")
lamin_file <- paste0(ref_dir, "lamin_B1_LADS2.bed")
dhs_file <- paste0(ref_dir, "DHS.bed")
time_file <- paste0(ref_dir, "lymph_rep_time.txt")
gc_file <- paste0(ref_dir, "gc10kb.bed")

# Annotate the data frame
df$exon <- binaryCol(df, exon_file)
df$CpGI <- binaryCol(df, cpg_island_file)
df$RR <- rcrCol(df, rcr_file, chr)
df$LAMIN <- binaryCol(df, lamin_file)
df$DHS <- binaryCol(df, dhs_file)
df$TIME <- repCol(df, time_file)
df$GC <- gcCol(df, gc_file)

write_csv(df, outfile)