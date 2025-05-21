#!/usr/bin/env Rscript
# ==============================================================================
# EWAStools-Based DNAm Preprocessing Script
# Script Name: preprocessingEwastools.R
# Description: Processes IDAT files using ewastools and performs quality
#              control, SNP handling, genotype calling, and leukocyte estimation.
# ==============================================================================
# Usage Example (Full version):
# ==============================================================================
# Rscript preprocessingEwastools.R \
#   --phenoFile data/preprocessingMinfi/pheno.csv \
#   --idatFolder data/idats/ \
#   --outputLogs logs/preprocessingEwastools/ \
#   --sampleIdCols Basename \
#   --matchColumn Basename \
#   --detectionPcutoff 0.01 \
#   --lcReference salivaEPIC \
#   --snpProbeType rs \
#   --outlierThreshold -4 \
#   --outPhenoQC data/preprocessingEwastools/pheno_ewasQC.csv
# ==============================================================================
# Usage Example (Default values with key parameters):
# ==============================================================================
# Rscript preprocessingEwastools.R \
#   --phenoFile data/preprocessingMinfi/pheno.csv \
#   --idatFolder data/idats/

# ==============================================================================
# DNAm Preprocessing Script — Input Arguments
# ==============================================================================
#   --phenoFile         [FILE]   Path to phenotype/sample sheet (CSV or TSV)  
#   --idatFolder        [DIR]    Path to raw .idat files  
#   --outputLogs        [DIR]    Folder to write logs and messages  
#   --sampleIdCols      [STR]    Columns to combine into SampleID (comma-separated)  
#   --matchColumn       [STR]    Column used to match IDAT files and metadata  
#   --detectionPcutoff  [NUM]    Detection P-value cutoff (e.g., 0.01)  
#   --lcReference       [STR]    Reference panel for leukocyte estimation (e.g., "salivaEPIC")  
#   --snpProbeType      [STR]    Probe prefix for SNP filtering (e.g., "rs")  
#   --outlierThreshold  [NUM]    Threshold to flag sample SNP outliers (e.g., -4)  
#   --outPhenoQC        [FILE]   Path to save the phenotype file with QC annotations  
# ==============================================================================

# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(optparse)
        library(data.table)
        library(ewastools)
        library(magrittr)
})

# ----------- Define Options -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--phenoFile", default = "data/preprocessingMinfi/pheno.csv", help = "Phenotype file"),
        make_option("--idatFolder", default = "data/idats/", help = "IDAT folder"),
        make_option("--outputLogs", default = "logs/", help = "Log output folder"),
        make_option("--scriptLabel", default = "preprocessingEwastools", help = "Script label used in log filename"),
        make_option("--sampleIdCols", default = "Basename", help = "Columns to construct SampleID (comma-separated)"),
        make_option("--matchColumn", default = "Basename", help = "Column in pheno used to match meth$meta$sample_id"),
        make_option("--detectionPcutoff", type = "double", default = 0.01, help = "Detection P-value cutoff"),
        make_option("--lcReference", default = "salivaEPIC", help = "Reference used in estimateLC()"),
        make_option("--snpProbeType", default = "rs", help = "SNP probe type for filtering"),
        make_option("--outlierThreshold", type = "double", default = -4, help = "Outlier threshold for SNP"),
        make_option("--outPhenoQC", default = "data/preprocessingEwastools/pheno_ewasQC.csv", help = "Path to save phenotype file with QC annotations [default: %default]", metavar = "FILE")
        
)))

# ----------- Logging Setup -----------
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)

logFilePath <- file.path(opt$outputLogs, paste0("log_", opt$scriptLabel, ".txt"))
logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
#===============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting Preprocessing ====\n")
cat("Start time: ", format(Sys.time()), "\n\n")
cat("Phenotype file:       ", opt$phenoFile, "\n")
cat("IDAT folder:          ", opt$idatFolder, "\n")
cat("Logging output path:  ", opt$outputLogs, "\n")
cat("Sample ID columns:    ", opt$sampleIdCols, "\n")
cat("Match column:         ", opt$matchColumn, "\n")
cat("Detection P cutoff:   ", opt$detectionPcutoff, "\n")
cat("Leukocyte reference:  ", opt$lcReference, "\n")
cat("SNP probe type:       ", opt$snpProbeType, "\n")
cat("Outlier threshold:    ", opt$outlierThreshold, "\n\n")
cat("Output QC phenotype file: ", opt$outPhenoQC, "\n")
# =============================================================================

dir.create(dirname(opt$outPhenoQC), recursive = TRUE, showWarnings = FALSE)
cat("=======================================================================\n")

# ----------- Load Sample Sheet -----------
cat("Loading phenotype file:", opt$phenoFile, "\n")
pheno <- read.csv(opt$phenoFile)
cat("=======================================================================\n")

# ----------- Create SampleID -----------
sampleIdCols <- strsplit(opt$sampleIdCols, ",")[[1]]
pheno$SampleID <- do.call(paste, c(pheno[sampleIdCols], sep = ""))
pheno$SampleID <- file.path(opt$idatFolder, pheno$SampleID)
cat("=======================================================================\n")

# ----------- Read IDATs -----------
cat("Reading IDAT files from:", opt$idatFolder, "\n")
meth <- read_idats(pheno$SampleID, quiet = TRUE)
cat("=======================================================================\n")

# ----------- Match Sample Metadata -----------
pheno <- pheno[order(match(pheno[[opt$matchColumn]], meth$meta$sample_id)), ]
stopifnot(all(pheno[[opt$matchColumn]] == meth$meta$sample_id))
cat("=======================================================================\n")

# ----------- Control Metrics & Failure Flags -----------
cat("Running control metrics QC...\n")
ctrls <- control_metrics(meth)
pheno$failed <- sample_failure(ctrls)
fails <- subset(pheno, failed == "TRUE")
cat("Failed samples:", nrow(fails), "\n")
cat("=======================================================================\n")

# ----------- Estimate LC Composition -----------
cat("Estimating leukocyte composition...\n")
beta <- meth %>% 
        detectionP %>%  
        mask(opt$detectionPcutoff) %>% 
        correct_dye_bias %>% 
        dont_normalize
LC <- estimateLC(beta, ref = opt$lcReference, constrained = FALSE)
pheno <- cbind(pheno, LC)
print(round(head(LC), 3))
cat("=======================================================================\n")

# ----------- Genotype Calling and Outlier Detection -----------
snps <- meth$manifest[probe_type == opt$snpProbeType, index]
snps <- beta[snps, ]
stopifnot(all(grepl(opt$snpProbeType, rownames(snps))))

cat("Calling genotypes and detecting SNP outliers...\n")
genotypes <- call_genotypes(snps, learn = FALSE)
pheno$outlier <- snp_outliers(genotypes)

if (!is.data.table(pheno)) pheno <- data.table(pheno)

pheno[outlier > opt$outlierThreshold, .(SampleID,outlier)]

pheno$outlierYN <- pheno$outlier > opt$outlierThreshold

pheno$donor_id <- enumerate_sample_donors(genotypes)
pheno[, n := .N, by = donor_id]
pheno[n > 1, .(SampleID,donor_id)] 
cat("=======================================================================\n")

# ----------- Save Output -----------
cat("Saving output phenotype with QC annotations to:", opt$outPhenoQC, "\n")
write.csv(pheno, file = opt$outPhenoQC, row.names = FALSE)

cat("Done. Final phenotype rows:", nrow(pheno), "\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# =============================================================================
# ----------- Close Logging -----------
sink(type = "message")  
sink()                  
close(logCon)           

