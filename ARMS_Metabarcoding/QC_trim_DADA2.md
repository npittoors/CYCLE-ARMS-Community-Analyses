#!/bin/bash
# =============================================================================
# COI Metabarcoding Pipeline — Primer Trimming, Quality Control, and DADA2
# =============================================================================
# Author: Nicole Pittoors, Lehigh Oceans Research Center, Lehigh University
# Manuscript: "Environmental filtering shapes patch dynamics across isolated
#              mesophotic reefs"
#
# Amplicon target: ~313 bp fragment of mitochondrial COI (cytochrome c oxidase
#                  subunit I); ~390 bp including linkers, adapters, and primers
#
# Primers (Geller/Leray):
#   Forward 5'->3': GGWACWGGWTGAACWGTWTAYCCYCC
#   Reverse 5'->3': TAIACYTCIGGRTGICCRAARAAYCA  (contains inosine)
#
# References:
#   Geller et al. 2013, Mol Ecol Res 13(5):851-861
#   Leray et al. 2013, Front Zool 10(34):1-14
#
# Pipeline steps:
#   1. Quality control    — FastQC + MultiQC
#   2. Primer trimming    — Cutadapt
#   3. Denoising + ASVs   — DADA2 (R script embedded at end of this file)
#
# Usage:
#   Run interactively on the Herrera Lab Thelio server or any HPC environment.
#   Adjust paths in the CONFIGURATION section below before running.
#   DADA2 steps are in R — copy that section into an R session or .Rmd file.
#
# Note on DADA2 truncation lengths:
#   truncLen = c(220, 190) confirmed for this dataset. If re-running on new
#   data, inspect quality profiles first and adjust accordingly.
# =============================================================================

# =============================================================================
# CONFIGURATION — edit these paths before running
# =============================================================================

SEQS_DIR="/media/meta/ARMS/ncp/CYCLE_ARMS_COI_run1/miseq/seqs"
TRIMMED_DIR="${SEQS_DIR}/fastq-primers-trimmed"
WORKING_DIR="${SEQS_DIR}"

# =============================================================================
# STEP 1: QUALITY CONTROL
# =============================================================================

cd "${SEQS_DIR}"

mamba activate readsqc

# Run FastQC on all raw fastq files
fastqc ./*.fastq.gz >> fastqc.out

# Summarize with MultiQC
multiqc ./ -o multiqc_report --interactive

mamba deactivate

# =============================================================================
# STEP 2: PRIMER TRIMMING WITH CUTADAPT
# =============================================================================

mamba activate edna

mkdir -p "${TRIMMED_DIR}"
cd "${SEQS_DIR}"

# --- Test on a single sample first ---
# Uncomment and edit the sample name to verify primer trimming before looping.
#
# cutadapt \
#   -g ^GGWACWGGWTGAACWGTWTAYCCYCC \
#   -G ^TAIACYTCIGGRTGICCRAARAAYCA \
#   -n 1 \
#   -O 5 \
#   --action trim --discard-untrimmed \
#   --minimum-length 5 \
#   --pair-filter=any \
#   -o "${TRIMMED_DIR}/TEST_R1_trimmed.fastq.gz" \
#   -p "${TRIMMED_DIR}/TEST_R2_trimmed.fastq.gz" \
#   TEST_R1_001.fastq.gz TEST_R2_001.fastq.gz

# --- Run primer trimming on all samples ---
for i in *_R1_001.fastq.gz; do
    SAMPLE=$(echo "${i}" | sed "s/_R1_001\.fastq\.gz//")
    echo "Trimming: ${SAMPLE}"
    cutadapt \
        -g ^GGWACWGGWTGAACWGTWTAYCCYCC \
        -G ^TAIACYTCIGGRTGICCRAARAAYCA \
        -n 1 \
        -O 5 \
        --action trim --discard-untrimmed \
        --minimum-length 5 \
        --pair-filter=any \
        -o "${TRIMMED_DIR}/primers-trimmed-${SAMPLE}_R1_001.fastq.gz" \
        -p "${TRIMMED_DIR}/primers-trimmed-${SAMPLE}_R2_001.fastq.gz" \
        "${SAMPLE}_R1_001.fastq.gz" \
        "${SAMPLE}_R2_001.fastq.gz"
done

mamba deactivate

# --- Unzip trimmed files and remove zero-size files ---
cd "${TRIMMED_DIR}"
gunzip *.fastq.gz
find . -name '*.fastq' -size 0 -print0 | xargs -0 rm
cd ..

# =============================================================================
# STEP 3: DADA2 — QUALITY FILTERING, DENOISING, AND ASV TABLE GENERATION
# =============================================================================
# Switch to R for the following steps.
# Copy this block into an R session or .Rmd file.
# Adjust the working_directory path before running.
# =============================================================================

: '
-------------------------------------------------------------------------------
### R CODE BEGINS HERE ###
-------------------------------------------------------------------------------

library(dada2)
library(tidyverse)
library(phyloseq)
library(decontam)
library(ggplot2)

# --- Paths ---
working_directory <- "/path/to/fastq-primers-trimmed"   # <-- edit this
setwd(working_directory)

fastq_path <- working_directory
fnFs <- sort(list.files(fastq_path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastq_path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# --- Visualize read quality (inspect before filtering) ---
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

# --- Quality filtering ---
# truncLen = c(220, 190) confirmed for this dataset.
# maxEE = c(2, 3) applied to forward and reverse reads respectively.
filtFs <- file.path(working_directory, "filtered", basename(fnFs))
filtRs <- file.path(working_directory, "filtered", basename(fnRs))

filtering_output <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen = c(220, 190),
    maxEE = c(2, 3),
    compress = TRUE,
    multithread = TRUE,
    verbose = TRUE
)
print(filtering_output)

# Remove any samples that did not survive filtering
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

# --- Learn error rates ---
# Uses loessErrfun (quality-score-ignoring error estimation)
errF <- learnErrors(filtFs, multithread = TRUE, randomize = TRUE,
                    verbose = TRUE, errorEstimationFunction = loessErrfun)
errR <- learnErrors(filtRs, multithread = TRUE, randomize = TRUE,
                    verbose = TRUE, errorEstimationFunction = loessErrfun)

# Inspect error models
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# --- Denoise ---
# Pseudopooling enables detection of low-abundance singletons across samples
ddF <- dada(filtFs, err = errF, pool = "pseudo", multithread = TRUE)
ddR <- dada(filtRs, err = errR, pool = "pseudo", multithread = TRUE)

# --- Merge paired reads ---
# justConcatenate = FALSE: reads are merged (require overlap), not concatenated
merged <- mergePairs(ddF, filtFs, ddR, filtRs,
                     verbose = TRUE, justConcatenate = FALSE)

# --- Construct sequence table and remove chimeras ---
table.chimeras <- makeSequenceTable(merged)
dim(table.chimeras)

table.no.chimeras <- removeBimeraDenovo(table.chimeras,
                                        multithread = TRUE, verbose = TRUE)
cat("Fraction of reads retained after chimera removal:",
    round(sum(table.no.chimeras) / sum(table.chimeras), 3), "\n")

# --- Export ASV table and representative sequences ---
write.table(t(table.no.chimeras), file = "table.tsv",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

uniquesToFasta(table.no.chimeras, fout = "rep-seqs.fasta",
               ids = colnames(table.no.chimeras))

# --- Rename ASVs (ASV_1, ASV_2, ...) and export clean tables ---
asv_seqs <- colnames(table.no.chimeras)
asv_headers <- paste0(">ASV_", seq_along(asv_seqs))

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(table.no.chimeras)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# --- Pipeline tracking summary ---
getN <- function(x) sum(getUniques(x))
track <- cbind(
    filtering_output[exists, ],
    sapply(ddF, getN),
    sapply(ddR, getN),
    sapply(merged, getN),
    rowSums(table.no.chimeras)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- basename(fnFs)
print(head(track))
write.csv(track, "pipeline_stats.csv")

### R CODE ENDS HERE ###
'

#!/bin/bash
# =============================================================================
# 18S rRNA V4 Metabarcoding Pipeline — Primer Trimming, Quality Control,
# and DADA2
# =============================================================================
# Author: Nicole Pittoors, Lehigh Oceans Research Center, Lehigh University
# Manuscript: "Environmental filtering shapes patch dynamics across isolated
#              mesophotic reefs"
#
# Amplicon target: ~536 bp fragment of the 18S V4 region
#
# Primers (V4_18SNext; Piredda et al. 2017; Tragin et al. 2018):
#   Full fusion primers used in PCR (Illumina Nextera adapter + biological):
#     V4_18SNext.For: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG[CCAGCASCYGCGGTAATTCC]
#     V4_18SNext.Rev: GTCTCGTGGGCTCGGAGATCTGTATAAGAGACAG[ACTTTCGTTCTTGATYRATGA]
#
#   Biological primer sequences trimmed by Cutadapt:
#     Forward: CCAGCASCYGCGGTAATTCC  (required, anchored)
#     Reverse: ACTTTCGTTCTTGATYRATGA (required, anchored)
#     Internal reverse complements also trimmed (optional, linked adapter syntax)
#
# References:
#   Piredda R et al. 2017, FEMS Microbiol Ecol 93(2)
#   Tragin M et al. 2018, PeerJ 6:e4781
#
# Pipeline steps:
#   1. Quality control    — FastQC + MultiQC
#   2. Primer trimming    — Cutadapt
#   3. Denoising + ASVs   — DADA2 (R script embedded at end of this file)
#
# Usage:
#   Run interactively on the Herrera Lab Thelio server or any HPC environment.
#   Adjust paths in the CONFIGURATION section below before running.
#   DADA2 steps are in R — copy that section into an R session or .Rmd file.
#
# Note on DADA2 truncation lengths:
#   Truncation lengths for 18S were not recoverable from retained code.
#   Inspect quality profiles for your data and set truncLen accordingly.
#   See manuscript methods for reported values used in this study.
# =============================================================================

# =============================================================================
# CONFIGURATION — edit these paths before running
# =============================================================================

SEQS_DIR="/path/to/seqs"                          # <-- edit this
TRIMMED_DIR="${SEQS_DIR}/fastq-primers-trimmed"

# =============================================================================
# STEP 1: QUALITY CONTROL
# =============================================================================

cd "${SEQS_DIR}"

mamba activate readsqc

# Run FastQC on all raw fastq files
fastqc ./*.fastq.gz >> fastqc.out

# Summarize with MultiQC — all reads, forward only, reverse only
multiqc . -o multiqc_report_all --interactive
multiqc . --ignore "*R2_001_fastqc.html" -o multiqc_report_forward --interactive
multiqc . --ignore "*R1_001_fastqc.html" -o multiqc_report_reverse --interactive

mamba deactivate

# =============================================================================
# STEP 2: PRIMER TRIMMING WITH CUTADAPT
# =============================================================================

mamba activate cutadaptenv

mkdir -p "${TRIMMED_DIR}"
cd "${SEQS_DIR}"

# --- Test on a single sample first ---
# Uncomment and edit the sample name to verify primer trimming before looping.
#
# cutadapt \
#   -a "^CCAGCASCYGCGGTAATTCC;required;...TCATYRATCAAGAACGAAAGT$;optional" \
#   -A "^ACTTTCGTTCTTGATYRATGA;required;...GGAATTACCGCRGSTGCTGG$;optional" \
#   -n 1 \
#   -O 5 \
#   --action trim --discard-untrimmed \
#   --minimum-length 5 \
#   --pair-filter=any \
#   -o "${TRIMMED_DIR}/TEST_R1_trimmed.fastq.gz" \
#   -p "${TRIMMED_DIR}/TEST_R2_trimmed.fastq.gz" \
#   TEST_R1_001.fastq.gz TEST_R2_001.fastq.gz
#
# Verify trimming was successful:
#   head -n 2 TEST_R1_001.fastq.gz
#   head -n 2 TEST_R1_trimmed.fastq.gz

# --- Run primer trimming on all samples ---
for i in *_R1_001.fastq.gz; do
    SAMPLE=$(echo "${i}" | sed "s/_R1_001\.fastq\.gz//")
    echo "Trimming: ${SAMPLE}"
    cutadapt \
        -a "^CCAGCASCYGCGGTAATTCC;required;...TCATYRATCAAGAACGAAAGT$;optional" \
        -A "^ACTTTCGTTCTTGATYRATGA;required;...GGAATTACCGCRGSTGCTGG$;optional" \
        -n 1 \
        -O 5 \
        --action trim --discard-untrimmed \
        --minimum-length 5 \
        --pair-filter=any \
        -o "${TRIMMED_DIR}/primers-trimmed-${SAMPLE}_R1_001.fastq.gz" \
        -p "${TRIMMED_DIR}/primers-trimmed-${SAMPLE}_R2_001.fastq.gz" \
        "${SAMPLE}_R1_001.fastq.gz" \
        "${SAMPLE}_R2_001.fastq.gz"
done

conda deactivate

# =============================================================================
# STEP 3: DADA2 — QUALITY FILTERING, DENOISING, AND ASV TABLE GENERATION
# =============================================================================
# Switch to R for the following steps.
# Copy this block into an R session or .Rmd file.
# Adjust the working_directory path before running.
#
# NOTE: truncLen values are not confirmed for 18S — inspect quality profiles
# and set appropriate truncation lengths before running filterAndTrim().
# =============================================================================

: '
-------------------------------------------------------------------------------
### R CODE BEGINS HERE ###
-------------------------------------------------------------------------------

library(dada2)
library(tidyverse)
library(phyloseq)
library(decontam)
library(ggplot2)

# --- Paths ---
working_directory <- "/path/to/fastq-primers-trimmed"   # <-- edit this
setwd(working_directory)

fastq_path <- working_directory
fnFs <- sort(list.files(fastq_path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(fastq_path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# --- Visualize read quality BEFORE setting truncLen ---
# Inspect these plots carefully. The 18S V4 amplicon (~536 bp) is longer than
# COI (~313 bp), so truncation lengths will differ. Set truncLen to where
# quality scores drop appreciably (typically Q < 30).
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

# --- Quality filtering ---
# !! SET truncLen BASED ON QUALITY PROFILES ABOVE !!
# Values used in the original study are reported in the manuscript methods.
filtFs <- file.path(working_directory, "filtered", basename(fnFs))
filtRs <- file.path(working_directory, "filtered", basename(fnRs))

filtering_output <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    truncLen = c(220,180),
    maxEE = c(2, 3),
    compress = TRUE,
    multithread = TRUE,
    verbose = TRUE
)
print(filtering_output)

# Remove any samples that did not survive filtering
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

# --- Learn error rates ---
errF <- learnErrors(filtFs, multithread = TRUE, randomize = TRUE, verbose = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE, randomize = TRUE, verbose = TRUE)

# Inspect error models — points should follow the black line closely
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# --- Denoise ---
# Pseudopooling enables detection of low-abundance singletons across samples
ddF <- dada(filtFs, err = errF, pool = "pseudo", multithread = TRUE)
ddR <- dada(filtRs, err = errR, pool = "pseudo", multithread = TRUE)

# --- Merge paired reads ---
# justConcatenate = FALSE: reads are merged (require overlap), not concatenated
# Note: 18S V4 (~536 bp) with standard MiSeq 2x300 should have sufficient
# overlap for merging. If merge rates are poor, check truncLen settings.
merged <- mergePairs(ddF, filtFs, ddR, filtRs,
                     verbose = TRUE, justConcatenate = FALSE)

# --- Construct sequence table and remove chimeras ---
table.chimeras <- makeSequenceTable(merged)
dim(table.chimeras)

# Inspect amplicon length distribution — expect peak around 536 bp for 18S V4
table(nchar(getSequences(table.chimeras)))

table.no.chimeras <- removeBimeraDenovo(table.chimeras,
                                        multithread = TRUE, verbose = TRUE)
cat("Fraction of reads retained after chimera removal:",
    round(sum(table.no.chimeras) / sum(table.chimeras), 3), "\n")

# --- Export ASV table and representative sequences ---
write.table(t(table.no.chimeras), file = "table.tsv",
            sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

uniquesToFasta(table.no.chimeras, fout = "rep-seqs.fasta",
               ids = colnames(table.no.chimeras))

# --- Rename ASVs (ASV_1, ASV_2, ...) and export clean tables ---
asv_seqs <- colnames(table.no.chimeras)
asv_headers <- paste0(">ASV_", seq_along(asv_seqs))

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

asv_tab <- t(table.no.chimeras)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# --- Pipeline tracking summary ---
getN <- function(x) sum(getUniques(x))
track <- cbind(
    filtering_output[exists, ],
    sapply(ddF, getN),
    sapply(ddR, getN),
    sapply(merged, getN),
    rowSums(table.no.chimeras)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                     "merged", "nonchim")
rownames(track) <- basename(fnFs)
print(head(track))
write.csv(track, "pipeline_stats.csv")

### R CODE ENDS HERE ###
'
