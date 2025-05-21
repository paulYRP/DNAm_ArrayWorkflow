#!/usr/bin/env Rscript
# ==============================================================================
# Epigenetic Clock Integration Script (Timepoint 1 and 3)
# Script Name: epigeneticAge_T1T2.R
# Description: This script integrates epigenetic clock estimates (Methylclock, 
# Clock Foundation) with psychological phenotype data at two timepoints 
# (T1 and T2).
# ==============================================================================
# Usage Example (Full version):

# Rscript epigeneticAge_T1T2.R \
#   --outputLogs logs/epigeneticAge_T1T2 \
#   --betaT1 rData/preprocessingPheno/metrics/betaT1.RData \
#   --betaT2 rData/preprocessingPheno/metrics/betaT2.RData \
#   --phenoT1 data/preprocessingPheno/phenoT1.csv \
#   --phenoT2 data/preprocessingPheno/phenoT2.csv \
#   --clockDictT1 data/clockFundation/T1_DNAmAgeCalcProject_17761_DataDict.csv \
#   --clockResT1 data/clockFundation/T1_DNAmAgeCalcProject_17761_Results.csv \
#   --clockDictT2 data/clockFundation/T2_DNAmAgeCalcProject_17685_DataDict.csv \
#   --clockResT2 data/clockFundation/T2_DNAmAgeCalcProject_17685_Results.csv \
#   --minfiPheno data/preprocessingMinfi/pheno.csv \
#   --phenoSep "\t" \
#   --idRenameFrom SID \
#   --idRenameTo id \
#   --sexZeroValue 0 \
#   --sex0Label Female \
#   --sex1Label Male \
#   --phenotypes DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore \
#   --covariates Sex,Age,Ethnicity,TraumaDefinition,Leukocytes.EWAS,Epithelial.cells.EWAS \
#   --factorVars Sex,Ethnicity,TraumaDefinition \
#   --clockPattern "\\.Methylclock$|\\.DNAm_Age\\.ClockF$" \
#   --prsMap PCL_SUM:PTSD_PRS,PCL5_B:PTSD_PRS,PCL5_C:PTSD_PRS,PCL5_D:PTSD_PRS,PCL5_E:PTSD_PRS,PTGIX_SUM:PTSD_PRS,DASS_D:MDD_PRS,DASS_S:PTSD_PRS,DASS_A:GAD_PRS,SSS8_SUM:MDD_PRS \
#   --modelOutDir rData/epigeneticAge_T1T2/models \
#   --dnAmAgeLM_T1Out rData/epigeneticAge_T1T2/models/phenoDNAmAgeLM_T1.RData \
#   --dnAmAgeLM_T2Out rData/epigeneticAge_T1T2/models/phenoDNAmAgeLM_T2.RData \
#   --outputPlots figures/epigeneticAge_T1T2 \
#   --plotWidth 2000 --plotHeight 1000 --plotDPI 150
# ==============================================================================

# Usage Example (Default values with key parameters):

# Rscript epigeneticAge_T1T2.R \
#   --betaT1 rData/preprocessingPheno/metrics/betaT1.RData \
#   --betaT2 rData/preprocessingPheno/metrics/betaT2.RData
# ==============================================================================

suppressPackageStartupMessages({
        library(optparse)
        library(rstudioapi)
        library(tidyverse)
        library(readr)
        library(dplyr)
        library(stringr)
        library(purrr)
        library(methylclock)

})
# ==============================================================================

# DNAm Epigenetic Age Analysis — Input Arguments

#   --outputLogs             [DIR]    Directory to save log file
#   --betaT1                 [RDATA]  RData file with DNAm beta values at Timepoint 1
#   --betaT2                 [RDATA]  RData file with DNAm beta values at Timepoint 3
#   --phenoT1                [CSV]    CSV file with phenotype metadata for T1
#   --phenoT2                [CSV]    CSV file with phenotype metadata for T2
#   --clockDictT1            [CSV]    Data dictionary file for Clock Foundation results (T1)
#   --clockResT1             [CSV]    Clock Foundation results file (T1)
#   --clockDictT2            [CSV]    Data dictionary file for Clock Foundation results (T2)
#   --clockResT2             [CSV]    Clock Foundation results file (T2)
#   --minfiPheno             [CSV]    Original phenotype metadata from Minfi preprocessing
#   --phenoSep               [CHAR]   Field separator used in phenotype CSV (default: tab)
#   --idRenameFrom           [STR]    Column name to be renamed as ID (e.g., "SID")
#   --idRenameTo             [STR]    New column name used by Clock calculators (e.g., "id")
#   --sexZeroValue           [STR]    Numeric value representing Female (default: "0")
#   --sex0Label              [STR]    Label for Female (default: "Female")
#   --sex1Label              [STR]    Label for Male (default: "Male")
#   --probeIdCol             [STR]    Column name for CpG identifiers [default: ProbeID]
#   --colPatternT1           [REGEX]  Pattern to clean column names from betaT1
#   --colPatternT2           [REGEX]  Pattern to clean column names from betaT2
#   --rowPattern             [REGEX]  Pattern to clean rownames (sample IDs)
#   --zipFlags               [STR]    Flags to be passed to `zip()` function
#   --epigeneticDir          [DIR]    Directory to save intermediate beta CSVs and ZIPs
#   --csvT1                  [CSV]    CSV output for betaT1
#   --csvT2                  [CSV]    CSV output for betaT2
#   --zipT1                  [ZIP]    ZIP output for betaT1
#   --zipT2                  [ZIP]    ZIP output for betaT2
#   --phenoClockF            [CSV]    Clock-compatible phenotype file after ID/Sex adjustments
#   --dropClockCols          [INT]    Index range to remove non-clock columns (e.g., "2:19")
#   --clockJoinId            [STR]    Column in clock results used for join
#   --clockPhenoJoinCol      [STR]    Column in pheno used for join (usually "SID")
#   --clockCleanPattern      [REGEX]  Pattern to clean special characters from clock column names
#   --methylclockSuffix      [STR]    Suffix to append to methylclock column names
#   --methylclockJoinCol     [STR]    Column in methylclock data to join with pheno
#   --phenoDNAmAgeT1Out      [CSV]    Final output for T1 pheno + DNAm age
#   --phenoDNAmAgeT2Out      [CSV]    Final output for T2 pheno + DNAm age
#   --phenotypes             [STR]    Comma-separated list of mental health phenotypes
#   --covariates             [STR]    Comma-separated list of covariates to include in LM
#   --factorVars             [STR]    Comma-separated list of variables to convert to factor
#   --clockPattern           [REGEX]  Regex pattern to detect clock columns
#   --prsMap                 [STR]    Optional: Comma-separated mapping of phenotype to PRS covariates
#   --timeLabelT1            [STR]    Label for Timepoint 1
#   --timeLabelT2            [STR]    Label for Timepoint 3
#   --modelOutDir            [DIR]    Output directory for GLM model objects
#   --dnAmAgeLM_T1Out        [RDATA]  Output file path for linear model results (T1)
#   --dnAmAgeLM_T2Out        [RDATA]  Output file path for linear model results (T2)
#   --modelSummaryDir        [DIR]    Output directory for summary text files
#   --outputPlots            [DIR]    Directory for saving diagnostic plots
#   --plotWidth              [INT]    Width (in pixels) for diagnostic plots
#   --plotHeight             [INT]    Height (in pixels) for diagnostic plots
#   --plotDPI                [INT]    Resolution (DPI) for diagnostic plots
# ==============================================================================

opt <- parse_args(OptionParser(option_list = list(
        make_option("--outputLogs", default = "logs/epigeneticAge_T1T2", help = "Directory to save log file [default: %default]"),
        make_option("--betaT1", default = "rData/preprocessingPheno/metrics/betaT1.RData", help = "RData file with DNAm beta values for Timepoint 1 [default: %default]"),
        make_option("--betaT2", default = "rData/preprocessingPheno/metrics/betaT2.RData", help = "RData file with DNAm beta values for Timepoint 2 [default: %default]"),
        make_option("--phenoT1", default = "data/preprocessingPheno/phenoT1.csv", help = "CSV with phenotype metadata for T1 [default: %default]"),
        make_option("--phenoT2", default = "data/preprocessingPheno/phenoT2.csv", help = "CSV with phenotype metadata for T2 [default: %default]"),
        make_option("--clockDictT1", default = "data/clockFundation/T1_DNAmAgeCalcProject_17761_DataDict.csv", help = "CSV dictionary for T1 epigenetic age results [default: %default]"),
        make_option("--clockResT1", default = "data/clockFundation/T1_DNAmAgeCalcProject_17761_Results.csv", help = "CSV results for T1 epigenetic age [default: %default]"),
        make_option("--clockDictT2", default = "data/clockFundation/T2_DNAmAgeCalcProject_17685_DataDict.csv", help = "CSV dictionary for T2 epigenetic age results [default: %default]"),
        make_option("--clockResT2", default = "data/clockFundation/T2_DNAmAgeCalcProject_17685_Results.csv", help = "CSV results for T2 epigenetic age [default: %default]"),
        make_option("--probeIdCol", default = "ProbeID", help = "Column name to use for CpG ID rownames [default: %default]"),
        make_option("--colPatternT1", default = "\\.1$", help = "Pattern to clean suffix from T1 column names [default: %default]"),
        make_option("--colPatternT2", default = "\\.2$", help = "Pattern to clean suffix from T2 column names [default: %default]"),
        make_option("--rowPattern", default = "_.*$", help = "Pattern to remove suffix from rownames [default: %default]"),
        make_option("--zipFlags", default = "-j", help = "Flags passed to `zip()` function [default: %default]"),
        make_option("--epigeneticDir", default = "data/epigeneticAge_T1T2", help = "Directory to save csv and zip files for epigenetic age [default: %default]"),
        make_option("--csvT1", default = "data/epigeneticAge_T1T2/betaT1.csv", help = "CSV file output for betaT1 [default: %default]"),
        make_option("--csvT2", default = "data/epigeneticAge_T1T2/betaT2.csv", help = "CSV file output for betaT2 [default: %default]"),
        make_option("--zipT1", default = "data/epigeneticAge_T1T2/betaT1.zip", help = "ZIP file output for betaT1 [default: %default]"),
        make_option("--zipT2", default = "data/epigeneticAge_T1T2/betaT2.zip", help = "ZIP file output for betaT2 [default: %default]"),
        make_option("--minfiPheno", default = "data/preprocessingMinfi/pheno.csv", help = "Original phenotype file from Minfi preprocessing [default: %default]"),
        make_option("--phenoSep", default = NULL, help = "Field separator for phenotype file: ',' for comma, ';' for semicolon, '\\t' for tab [default: %default]"),
        make_option("--idRenameFrom", default = "SID", help = "Column to rename as ID for Horvath [default: %default]"),
        make_option("--idRenameTo", default = "id", help = "New column name to use for ID (Horvath) [default: %default]"),
        make_option("--sex0Label", default = "Female",  help = "Label to assign when Sex == 0 [default: %default]"),
        make_option("--sex1Label", default = "Male",  help = "Label to assign when Sex == 1 [default: %default]"),
        make_option("--sexZeroValue", default = "0",  help = "Numeric or string value indicating 'Female' in Sex column [default: %default]"),
        make_option("--phenoClockF", default = "data/epigeneticAge_T1T2/phenoClockF.csv", help = "Output CSV for Horvath-compatible phenotype file [default: %default]"),
        make_option("--dropClockCols", default = "2:19", help = "Index range to drop columns from clock output [default: %default]"),
        make_option("--clockJoinId", default = "SID.SampleAnnotation.ClockF", help = "Column name in clock results to join with SID in phenotype [default: %default]"),
        make_option("--clockCleanPattern", default = "[\\.\\:\\ ]", help = "Pattern to clean clock column names [default: %default]"),
        make_option("--clockPhenoJoinCol", default = "SID", help = "Column in pheno data to join to clock results [default: %default]"),
        make_option("--methylclockSuffix", default = ".Methylclock", help = "Suffix for methylclock columns [default: %default]"),
        make_option("--methylclockJoinCol", default = "id.Methylclock", help = "Column in methylclock result to join by [default: %default]"),
        make_option("--phenoDNAmAgeT1Out", default = "data/epigeneticAge_T1T2/phenoDNAmAgeT1.csv", help = "Path to save DNAm Age + pheno for T1 [default: %default]"),
        make_option("--phenoDNAmAgeT2Out", default = "data/epigeneticAge_T1T2/phenoDNAmAgeT2.csv", help = "Path to save DNAm Age + pheno for T2 [default: %default]"),
        make_option("--phenotypes", type = "character", default = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore", help = "Comma-separated list of mental health phenotypes [default: %default]"),
        make_option("--covariates", type = "character", default = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes.EWAS,Epithelial.cells.EWAS", help = "Comma-separated list of covariates to include in GLM [default: %default]"),
        make_option("--factorVars", default = "Sex,Ethnicity,TraumaDefinition", help = "Variables to convert to factor [default: %default]"),
        make_option("--clockPattern", default = "\\.Methylclock$|\\.DNAm_Age\\.ClockF$", help = "Regex pattern to identify clock columns [default: %default]"),
        make_option("--prsMap", default = NULL, help = "Optional: comma-separated mapping of phenotype to PRS covariates (e.g., phenotype:PRS)"),
        make_option("--timeLabelT1", default = "T1", help = "Time label for Timepoint 1 [default: %default]"),
        make_option("--timeLabelT2", default = "T2", help = "Time label for Timepoint 2 [default: %default]"),
        make_option("--modelOutDir", default = "rData/epigeneticAge_T1T2/models", help = "Output directory for GLM model objects [default: %default]"),
        make_option("--dnAmAgeLM_T1Out", default = "rData/epigeneticAge_T1T2/models/phenoDNAmAgeLM_T1.RData", help = "Output file path for T1 DNAm Age GLM results [default: %default]"),
        make_option("--dnAmAgeLM_T2Out", default = "rData/epigeneticAge_T1T2/models/phenoDNAmAgeLM_T2.RData", help = "Output file path for T2 DNAm Age GLM results [default: %default]"),
        make_option("--modelSummaryDir", default = "preliminaryResults/epigeneticAge_T1T2", help = "Base directory to save per-clock model summary text files [default: %default]"),
        make_option("--outputPlots", default = "figures/epigeneticAge_T1T2", help = "Directory to save diagnostic plots [default: %default]"),
        make_option("--plotWidth", default = 2000, type = "integer", help = "Plot width in pixels [default: %default]"),
        make_option("--plotHeight", default = 1000, type = "integer", help = "Plot height in pixels [default: %default]"),
        make_option("--plotDPI", default = 150, type = "integer", help = "Plot DPI (resolution) [default: %default]")

)))

opt$covariateList <- strsplit(opt$covariates, ",")[[1]]
opt$phenotypeList <- strsplit(opt$phenotypes, ",")[[1]]
opt$factorVarsList <- strsplit(opt$factorVars, ",")[[1]]

if (!is.null(opt$prsMap)) {
        opt$prsMapList <- setNames(
                sapply(strsplit(unlist(strsplit(opt$prsMap, ",")), ":"), `[`, 2),
                sapply(strsplit(unlist(strsplit(opt$prsMap, ",")), ":"), `[`, 1)
        )
} else {
        opt$prsMapList <- list()
}

dir.create(opt$epigeneticDir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$modelOutDir, recursive = TRUE, showWarnings = FALSE)

summaryDir <- file.path(opt$modelSummaryDir, "summary")

dir.create(summaryDir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)
# ==============================================================================

# ----------- Setup Logging -----------
logFilePath <- file.path(opt$outputLogs, "log_epigeneticAge_T1T2.txt")
logCon <- file(logFilePath, open = "wt")
sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("=======================================================\n")
cat("Epigenetic Age Script Log\n")
cat("Start Time:", as.character(Sys.time()), "\n")
cat("Output Directory:", opt$outputLogs, "\n")
cat("Log File Path:", logFilePath, "\n")
cat("Column Pattern T1:", opt$colPatternT1, "\n")
cat("Column Pattern T2:", opt$colPatternT2, "\n")
cat("Row Pattern:", opt$rowPattern, "\n")
cat("Probe ID Column:", opt$probeIdCol, "\n")
cat("Column to Rename as ID:", opt$idRenameFrom, "\n")
cat("Mental Health Variables:", opt$phenotypes, "\n")
cat("Covariate Variables:", opt$covariates, "\n")
cat("Factor variables: ", opt$factorVars, "\n")
cat("Clock Join ID:", opt$clockJoinId, "\n")
cat("Clock Phenotype Join Column:", opt$clockPhenoJoinCol, "\n")
cat("Clock Clean Pattern:", opt$clockCleanPattern, "\n")
cat("Methylclock Suffix:", opt$methylclockSuffix, "\n")
cat("Methylclock Join Column:", opt$methylclockJoinCol, "\n")
cat("Drop Clock Columns:", opt$dropClockCols, "\n")
cat("Clock Pattern:", opt$clockPattern, "\n")
cat("Time Label T1:", opt$timeLabelT1, "\n")
cat("Time Label T2:", opt$timeLabelT2, "\n")
cat("Model Output Directory:", opt$modelOutDir, "\n")
cat("DNAm Age LM T1 Output:", opt$dnAmAgeLM_T1Out, "\n")
cat("DNAm Age LM T2 Output:", opt$dnAmAgeLM_T2Out, "\n")
cat("Model Summary Directory:", opt$modelSummaryDir, "\n")
cat("Output Plots Directory:", opt$outputPlots, "\n")
cat("=======================================================\n")
cat("=======================================================\n\n")
# ----------- Load DNAm Beta Matrices -----------
load(opt$betaT1)
load(opt$betaT2)

# ----------- Load Phenotype Metadata -----------
phenoT1 <- read.csv(opt$phenoT1)
cat("T1 phenotype dimensions:\n"); print(dim(phenoT1))
cat("T1 phenotype preview:\n"); print(head(phenoT1))
cat("=========================\n")

phenoT2 <- read.csv(opt$phenoT2)
cat("T2 phenotype dimensions:\n"); print(dim(phenoT2))
cat("T2 phenotype preview:\n"); print(head(phenoT2))
cat("=========================\n")

# ----------- Load Epigenetic Clock Outputs -----------
DNAmAgeDicT1 <- read.csv(opt$clockDictT1)
DNAmAgeCalT1 <- read.csv(opt$clockResT1)

DNAmAgeDicT2 <- read.csv(opt$clockDictT2)
DNAmAgeCalT2 <- read.csv(opt$clockResT2)

# -- Clock Summary Display --
cat("Clock T1 Dictionary:\n"); print(dim(DNAmAgeDicT1)); print(colnames(DNAmAgeDicT1)); print(head(DNAmAgeDicT1))
cat("=========================\n")
cat("Clock T1 Results:\n"); print(dim(DNAmAgeCalT1)); print(colnames(DNAmAgeCalT1)); print(head(DNAmAgeCalT1))
cat("=========================\n")
cat("Clock T2 Dictionary:\n"); print(dim(DNAmAgeDicT2)); print(colnames(DNAmAgeDicT2)); print(head(DNAmAgeDicT2))
cat("=========================\n")
cat("Clock T2 Results:\n"); print(dim(DNAmAgeCalT2)); print(colnames(DNAmAgeCalT2)); print(head(DNAmAgeCalT2))
cat("=======================================================================\n")

# ----------- Preprocessing Betas for Horvath Calculator -----------

cat("Preprocessing beta matrices for Horvath clock format...\n")

# Clean column names
colnames(betaT1) <- str_replace(colnames(betaT1), opt$colPatternT1, "")
colnames(betaT2) <- str_replace(colnames(betaT2), opt$colPatternT2, "")

# Clean row names
rownames(betaT1) <- str_replace(rownames(betaT1), opt$rowPattern, "")
rownames(betaT2) <- str_replace(rownames(betaT2), opt$rowPattern, "")

# Convert to data frame and add Probe ID column
betaT1csv <- betaT1 |> as.data.frame() |>
        tibble::rownames_to_column(var = opt$probeIdCol)

betaT2csv <- betaT2 |> as.data.frame() |>
        tibble::rownames_to_column(var = opt$probeIdCol)

# Save CSV files
write_csv(betaT1csv, opt$csvT1)
write_csv(betaT2csv, opt$csvT2)

# Save ZIP files (no nested folder structure)
zip(zipfile = opt$zipT1, files = opt$csvT1, flags = opt$zipFlags)
zip(zipfile = opt$zipT2, files = opt$csvT2, flags = opt$zipFlags)

# Summary output
cat("T1 CSV + ZIP:\n"); print(dim(betaT1)); print(head(betaT1))
cat("=========================\n")
cat("T2 CSV + ZIP:\n"); print(dim(betaT2)); print(head(betaT2))
cat("=======================================================================\n")

# ----------- Preprocessing CSV for Horvath Calculator -----------

cat("Preprocessing phenotype file for Horvath calculator...\n")

# Validate and convert separator argument
if (!is.null(opt$phenoSep)) {
        validSeps <- c(",", ";", "\\t")
        if (!(opt$phenoSep %in% validSeps)) {
                stop(paste0("Invalid separator: '", opt$phenoSep, 
                            "'. Allowed values: ',', ';', '\\t'"))
        }
        sepChar <- if (opt$phenoSep == "\\t") "\t" else opt$phenoSep
} else {
        sepChar <- NULL
}

# Load the phenotype file
pheno <- if (!is.null(sepChar)) {
        read.csv(opt$minfiPheno, sep = sepChar)
} else {
        read.csv(opt$minfiPheno)
}

# Rename specified ID column
pheno <- pheno %>%
        rename(!!opt$idRenameTo := !!sym(opt$idRenameFrom))

# Convert sex column based on dynamic comparison value
pheno$Sex <- ifelse(as.character(pheno$Sex) == as.character(opt$sexZeroValue),
                    opt$sex0Label, opt$sex1Label)

# Save the formatted phenotype file
write.csv(pheno, file = opt$phenoClockF, row.names = FALSE)

cat("Saved cleaned phenotype for Horvath clock:", opt$phenoClockF, "\n")
cat("=======================================================================\n")

# ----------- Preprocess Clock Results and Join to Phenotypes -----------

cat("Preprocessing Horvath clock results...\n")

# Function to rename clock columns using dictionary metadata
renameClockCol <- function(clockData, clockDict, pattern) {
        cleanedNames <- gsub(pattern, "", colnames(clockData))
        newNames <- paste0(cleanedNames, ".", clockDict$FieldType, ".ClockF")
        colnames(clockData) <- newNames
        return(clockData)
}

# Apply to both T1 and T2 with configurable pattern
DNAmAgeCalT1 <- renameClockCol(DNAmAgeCalT1, DNAmAgeDicT1, opt$clockCleanPattern)
DNAmAgeCalT2 <- renameClockCol(DNAmAgeCalT2, DNAmAgeDicT2, opt$clockCleanPattern)

# Dynamically drop specified columns
dropCols <- eval(parse(text = paste0("c(", opt$dropClockCols, ")")))
DNAmAgeCalT1 <- DNAmAgeCalT1[, -dropCols]
DNAmAgeCalT2 <- DNAmAgeCalT2[, -dropCols]

# Join to phenotypes using configurable phenotype join column
phenoT1 <- phenoT1 %>%
                left_join(DNAmAgeCalT1, by = setNames(opt$clockJoinId, opt$clockPhenoJoinCol))

phenoT2 <- phenoT2 %>%
        left_join(DNAmAgeCalT2, by = setNames(opt$clockJoinId, opt$clockPhenoJoinCol))

# Preview structure
cat("Joined column names for phenoT1:\n")
print(colnames(phenoT1))
cat("========================================\n")
cat("Joined column names for phenoT2:\n")
print(colnames(phenoT2))
cat("=======================================================================\n")

# ----------- Methylclock Age Calculation and Merge -----------
cat("Running Methylclock age estimation...\n")

# Check CpGs required by clocks
betaT1missCpGs <- checkClocks(betaT1)
betaT2missCpGs <- checkClocks(betaT2)

# Run methylclock prediction
DNAmAgeT1 <- DNAmAge(betaT1)
DNAmAgeT2 <- DNAmAge(betaT2)


cat("Removing unnecessary columns from DNAmAgeT1 and DNAmAgeT2...\n")
# Remove unnecessary columns from DNAmAgeT1 and DNAmAgeT2
DNAmAgeT1 <- DNAmAgeT1 %>%
        dplyr::select(-c(BNN, TL))

DNAmAgeT2 <- DNAmAgeT2 %>%
        dplyr::select(-c(BNN, TL))

cat("Dim DNAmAgeT1:\n")
print(dim(DNAmAgeT1))
cat("=========================\n")
print(head(DNAmAgeT1))

cat("Dim DNAmAgeT2:\n")
print(dim(DNAmAgeT2))
cat("=========================\n")
print(head(DNAmAgeT2))

# Rename columns to include .Methylclock suffix
colnames(DNAmAgeT1) <- paste0(colnames(DNAmAgeT1), opt$methylclockSuffix)
colnames(DNAmAgeT2) <- paste0(colnames(DNAmAgeT2), opt$methylclockSuffix)

cat("Colnames DNAmAgeT1:\n")
print(colnames(DNAmAgeT1))
cat("=========================\n")
cat("Colnames DNAmAgeT2:\n")
print(colnames(DNAmAgeT2))
# ==============================================================================

# Merge DNAmAge with phenotype data
phenoDNAmAgeT1 <- inner_join(
        phenoT1,
        DNAmAgeT1,
        by = setNames(opt$methylclockJoinCol, opt$clockPhenoJoinCol)
)

phenoDNAmAgeT2 <- inner_join(
        phenoT2,
        DNAmAgeT2,
        by = setNames(opt$methylclockJoinCol, opt$clockPhenoJoinCol)
)

cat("Merged phenoDNAmAgeT1 dim:\n")
print(dim(phenoDNAmAgeT1))
cat("=========================\n")
cat("Merged phenoDNAmAgeT2 dim:\n")
print(dim(phenoDNAmAgeT2))

# Save outputs
write.csv(phenoDNAmAgeT1, file = opt$phenoDNAmAgeT1Out, row.names = FALSE)
write.csv(phenoDNAmAgeT2, file = opt$phenoDNAmAgeT2Out, row.names = FALSE)
cat("=======================================================================\n")

# ----------- Epigenetic Age GLM Model Function -----------

cat("Converting factor variables to factors...\n")
for (var in opt$factorVarsList) {
        if (var %in% colnames(phenoDNAmAgeT1)) {
                phenoDNAmAgeT1[[var]] <- as.factor(phenoDNAmAgeT1[[var]])
                cat("T1 Levels for", var, ":\n")
                print(levels(phenoDNAmAgeT1[[var]]))
                cat("=========================\n")
        }
        if (var %in% colnames(phenoDNAmAgeT2)) {
                phenoDNAmAgeT2[[var]] <- as.factor(phenoDNAmAgeT2[[var]])
                cat("T2 Levels for", var, ":\n")
                print(levels(phenoDNAmAgeT2[[var]]))
                cat("=========================\n")
        }
}

cat("Fitting epigenetic age models...\n")
epigeneticAgeModel <- function(
                dataFrame,
                timeLabel,
                clockList,
                phenotypes = opt$phenotypeList,
                covariates = opt$covariateList) 
        {
        modelResults <- list()
        
        for (clock in clockList) {
                for (phenotype in phenotypes) {
                        
                        modelFormula <- as.formula(
                                paste(clock, "~", phenotype, "+", 
                                      paste(covariates, collapse = "+"))
                        )
                        
                        model <- lm(modelFormula, data = dataFrame)
                        
                        phenotypeTerm <- broom::tidy(model) %>%
                                dplyr::filter(term == phenotype) %>%
                                dplyr::mutate(
                                        epigeneticClock = clock,
                                        timepoint = timeLabel,
                                        rSquared = summary(model)$r.squared
                                )
                        
                        resultKey <- paste(clock, phenotype, sep = "_")
                        
                        modelResults[[resultKey]] <- list(
                                model = model,
                                coef = phenotypeTerm,
                                residuals = residuals(model),
                                fittedValues = fitted.values(model),
                                summary = summary(model)
                        )
                }
        }
        return(modelResults)
}

# ----------- Prepare Clock List via Pattern Match -----------
clockList <- grep(opt$clockPattern,
                  colnames(phenoDNAmAgeT1),
                  value = TRUE)

# ----------- Fit Models for T1 and T2 -----------
cat("Fitting epigenetic age models for T1 and T2...\n")

for (pheno in opt$phenotypeList) {
        
        cat("Running LM for:", pheno, "\n")
        
        # Dynamic covariate inclusion
        prsVar <- if (pheno %in% names(opt$prsMapList)) opt$prsMapList[[pheno]] else NULL
        allCovariates <- if (!is.null(prsVar)) c(opt$covariateList, prsVar) else opt$covariateList
        
        # Log the formula
        modelFormula <- paste("~", pheno, "+", paste(allCovariates, collapse = " + "))
        cat("Formula:", modelFormula, "\n")
        
        # Run T1
        phenoDNAmAgeLM_T1[[pheno]] <- epigeneticAgeModel(
                dataFrame = phenoDNAmAgeT1,
                timeLabel = opt$timeLabelT1,
                clockList = clockList,
                phenotypes = pheno,
                covariates = allCovariates
        )
        
        # Run T2
        phenoDNAmAgeLM_T2[[pheno]] <- epigeneticAgeModel(
                dataFrame = phenoDNAmAgeT2,
                timeLabel = opt$timeLabelT2,
                clockList = clockList,
                phenotypes = pheno,
                covariates = allCovariates
        )
}


save(phenoDNAmAgeLM_T1, file = opt$dnAmAgeLM_T1Out)
save(phenoDNAmAgeLM_T2, file = opt$dnAmAgeLM_T2Out)
cat("Epigenetic age models saved: ", opt$modelOutDir, "\n")

cat("Clock list: ", length(names(phenoDNAmAgeLM_T1)) , "\n")
cat("=========================\n")
cat("First Clock in the list: ", names(phenoDNAmAgeLM_T1)[1], "\n")
cat("=========================\n")
cat("First model summary:\n")
print(phenoDNAmAgeLM_T1[[1]]$summary)
cat("=======================================================================\n")

# ----------- Save Model Summaries for Epigenetic Age Clocks -----------
saveClockSummaries <- function(
                modelList,
                label,
                outputBaseDir = opt$modelSummaryDir
) {
        outDir <- file.path(outputBaseDir, paste0("models", label))
        lapply(names(modelList), function(n) {
                outFile <- file.path(outDir, paste0(n, ".txt"))
                writeLines(capture.output(modelList[[n]]$summary), con = outFile)
        })
}

# ----------- Create Output Directories First -----------
labelT1 <- deparse(substitute(phenoDNAmAgeLM_T1))
labelT2 <- deparse(substitute(phenoDNAmAgeLM_T2))

outDirT1 <- file.path(opt$modelSummaryDir, paste0("models", labelT1))
outDirT2 <- file.path(opt$modelSummaryDir, paste0("models", labelT2))

dir.create(outDirT1, recursive = TRUE, showWarnings = FALSE)
dir.create(outDirT2, recursive = TRUE, showWarnings = FALSE)

# ----------- Save Model Summaries -----------
saveClockSummaries(phenoDNAmAgeLM_T1, label = labelT1)
saveClockSummaries(phenoDNAmAgeLM_T2, label = labelT2)

cat("Model summaries saved to: ", outDirT1, "\n")
cat("Model summaries saved to: ", outDirT2, "\n")
cat("=======================================================================\n")

# ----------- Extract and Save Epigenetic Clock Coefficients -----------
extractClockCoefficients <- function(modelList) {
        purrr::map_dfr(modelList, "coef")
}

coefT1 <- extractClockCoefficients(phenoDNAmAgeLM_T1)
coefT2 <- extractClockCoefficients(phenoDNAmAgeLM_T2)

head(coefT1)
head(coefT2)

# Save summary coefficient tables
write.table(
        coefT1,
        file = file.path(summaryDir, "summaryPhenoDNAmAgeLM_T1.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
        coefT2,
        file = file.path(summaryDir, "summaryPhenoDNAmAgeLM_T2.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
)
cat("Coefficient tables saved to: ", summaryDir, "\n")
cat("=======================================================================\n")

# ----------- Define Diagnostic Plot Function -----------
plotEpigeneticDiagnostics <- function(
                modelList,
                outputDir,
                timeLabel,
                plotWidth = opt$plotWidth,
                plotHeight = opt$plotHeight,
                plotDPI = opt$plotDPI
) {
        dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
        
        for (modelName in names(modelList)) {
                modelObj <- modelList[[modelName]]$model
                if (is.null(modelObj)) next
                
                clockName <- sub("_(\\w+)$", "", modelName)
                
                clockDir <- file.path(outputDir, clockName)
                dir.create(clockDir, recursive = TRUE, showWarnings = FALSE)
                
                outPath <- file.path(clockDir, paste0(modelName, "_",
                                                      timeLabel, ".tif"))
                
                tiff(filename = outPath,
                     width = plotWidth,
                     height = plotHeight,
                     res = plotDPI, type = "cairo")
                
                par(mfrow = c(1, 3))
                
                # Residuals vs Fitted
                plot(modelObj$fitted.values, resid(modelObj),
                     xlab = "Fitted Values", ylab = "Residuals",
                     main = "Residuals vs Fitted")
                abline(h = 0, col = "red", lty = 2)
                
                # Q-Q Plot
                qqnorm(resid(modelObj), main = "Q-Q Plot")
                qqline(resid(modelObj), col = "blue", lty = 2)
                
                # Histogram of residuals
                hist(residuals(modelObj), col = "lightblue",
                     main = "Histogram of Residuals",
                     xlab = "Residuals")
                
                dev.off()
                par(mfrow = c(1, 1))  
        }
}

# ----------- Create Output Directory for Plots -----------
# Apply to Timepoint 1
plotEpigeneticDiagnostics(
        modelList = phenoDNAmAgeLM_T1,
        outputDir = file.path(opt$outputPlots, labelT1),
        timeLabel = "T1"
)

# Apply to Timepoint 3
plotEpigeneticDiagnostics(
        modelList = phenoDNAmAgeLM_T2,
        outputDir = file.path(opt$outputPlots, labelT2),
        timeLabel = "T2"
)

cat("Diagnostic plots saved to: ", opt$outputPlots, "\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)




