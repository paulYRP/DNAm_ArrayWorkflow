# !/usr/bin/env Rscript
# ==============================================================================
# DNAm Phenotype Preprocessing Script
# Script Name: preprocessingPheno.R
# Description: Merges phenotype data with EWAS QC, subsets by timepoints,
#              and prepares beta/m matrices for downstream modeling.
# ==============================================================================
# Usage Example (Full version):
# ==============================================================================
# Rscript preprocessingPheno.R \
#   --phenoFile data/preprocessingMinfi/pheno.csv \
#   --phenoEWAS data/preprocessingEwastools/pheno_ewasQC.csv \
#   --betaPath rData/preprocessingMinfi/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
#   --mPath rData/preprocessingMinfi/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
#   --cnPath rData/preprocessingMinfi/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
#   --dropColumnsPhenoEWAS SampleID \
#   --colsToRenamePhenoEWAS failed,Leukocytes,Epithelial.cells,outlier,outlierYN,donor_id,n \
#   --mergeKey SID \
#   --factorVars Sex,Ethnicity,TraumaDefinition \
#   --factorPrefixes Sex,Ethn,TraD \
#   --timepoints 1,2,3 \
#   --combineTimepoints 1,3 \
#   --outputPheno data/preprocessingPheno/ \
#   --outputRData rData/preprocessingPheno \
#   --outputLogs logs/preprocessingPheno

# ==============================================================================
# Usage Example (Default values with key parameters):
# ==============================================================================
# Rscript preprocessingPheno.R \
#   --phenoFile data/preprocessingMinfi/pheno.csv \
#   --phenoEWAS data/preprocessingEwastools/pheno_ewasQC.csv

# ==============================================================================
# DNAm Preprocessing Script — Input Arguments
# ==============================================================================
#   --phenoFile              [FILE]   Path to primary phenotype CSV/TSV file
#   --phenoEWAS              [FILE]   Path to EWAS-based QC phenotype file
#   --betaPath               [FILE]   RData file containing beta matrix
#   --mPath                  [FILE]   RData file containing M-value matrix
#   --cnPath                 [FILE]   RData file containing CN matrix
#   --dropColumnsPhenoEWAS   [STR]    Columns to drop from phenoEWAS (comma-separated)
#   --colsToRenamePhenoEWAS  [STR]    Columns in phenoEWAS to append '.EWAS' suffix (comma-separated)
#   --mergeKey               [STR]    Column name to merge pheno and phenoEWAS (default: "SID")
#   --factorVars             [STR]    Variables to convert to factors (comma-separated)
#   --factorPrefixes         [STR]    Prefixes for factor levels (must match factorVars order)
#   --timepoints             [STR]    Comma-separated list of timepoints to split/subset (default: "1,2,3")
#   --combineTimepoints      [STR]    Timepoints to combine for longitudinal analysis (default: "1,3")
#   --outputPheno            [DIR]   Output CSV file for merged phenotype (default: "data/preprocessingPheno/")
#   --outputRData            [DIR]    Output folder for timepoint-subset and merged RData (default: "rData/preprocessingPheno")
#   --outputLogs             [DIR]    Folder to save logs (default: "logs/preprocessingPheno")
# ==============================================================================


# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(optparse)
        library(dplyr)
})

# ----------- Define Input Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--phenoFile", default = "data/preprocessingMinfi/pheno.csv", help = "Input phenotype CSV file [default: %default]", metavar = "FILE"),
        make_option("--phenoEWAS", default = "data/preprocessingEwastools/pheno_ewasQC.csv", help = "EWAS phenotype CSV with QC columns", metavar = "FILE"),
        make_option("--betaPath", default = "rData/preprocessingMinfi/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData", help = "Path to Beta matrix RData [default: %default]", metavar = "FILE"),
        make_option("--mPath", default = "rData/preprocessingMinfi/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData", help = "Path to M-values RData [default: %default]", metavar = "FILE"),
        make_option("--cnPath", default = "rData/preprocessingMinfi/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData", help = "Path to CN matrix RData [default: %default]", metavar = "FILE"),
        make_option("--dropColumnsPhenoEWAS", default = "SampleID", help = "Columns to drop from phenoEWAS [default: %default]", metavar = "COLS"),
        make_option("--colsToRenamePhenoEWAS", default = "failed,Leukocytes,Epithelial.cells,outlier,outlierYN,donor_id,n", help = "Comma-separated columns to rename with '.EWAS'", metavar = "COLS"),
        make_option("--mergeKey", default = "SID", help = "Column name used to merge pheno and phenoEWAS [default: %default]", metavar = "COL"),
        make_option("--factorVars", default = "Sex,Ethnicity,TraumaDefinition", help = "Comma-separated columns to convert to factor", metavar = "COLS"),
        make_option("--factorPrefixes", default = "Sex,Ethn,TraD", help = "Comma-separated prefixes to prepend to factor levels", metavar = "PREFS"),
        make_option("--timepoints", default = "1,2,3", help = "Timepoints to subset for beta and M values", metavar = "T1,T2,T3"),
        make_option("--combineTimepoints", default = "1,3", help = "Timepoints to combine for longitudinal analysis", metavar = "T1,Tn"),
        make_option("--outputPheno", default = "data/preprocessingPheno/", help = "Path to save final phenotype CSVs", metavar = "DIR"),
        make_option("--outputRData", default = "rData/preprocessingPheno/metrics", help = "Directory to save processed RData objects metrics", metavar = "DIR"),
        make_option("--outputRDataMerge", default = "rData/preprocessingPheno/mergeData", help = "Directory to save processed RData objects mergedata", metavar = "DIR"),
        make_option("--outputLogs", default = "logs/preprocessingPheno", help = "Directory for all log output [default: %default]", metavar = "DIR"),
        make_option("--scriptLabel", default = "preprocessingPheno", help = "Label for log file naming [default: %default]", metavar = "STR")
        
)))

dir.create(opt$outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputRDataMerge, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$outputPheno, recursive = TRUE, showWarnings = FALSE)
#===============================================================================

# ----------- Logging Setup -----------
logFilePath <- file.path(opt$outputLogs, paste0("log_", opt$scriptLabel, ".txt"))
logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
#===============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting Phenotype Preprocessing ====\n")
cat("Start Time: ", format(Sys.time()), "\n")
cat("Log file path: ", logFilePath, "\n\n")
cat("Phenotype file:           ", opt$phenoFile, "\n")
cat("EWAS file:                ", opt$phenoEWAS, "\n")
cat("Beta path:                ", opt$betaPath, "\n")
cat("M-values path:            ", opt$mPath, "\n")
cat("CN path:                  ", opt$cnPath, "\n")
cat("Merged phenotype output:  ", opt$outputPheno, "\n")
cat("RData output directory:   ", opt$outputRData, "\n\n")
cat("Drop EWAS columns:        ", opt$dropColumnsPhenoEWAS, "\n")
cat("Columns to rename:        ", opt$colsToRenamePhenoEWAS, "\n")
cat("Merge key:                ", opt$mergeKey, "\n")
cat("Factor columns:           ", opt$factorVars, "\n")
cat("Factor prefixes:          ", opt$factorPrefixes, "\n")
cat("Timepoints:               ", opt$timepoints, "\n")
cat("Combined timepoints:      ", opt$combineTimepoints, "\n")
cat("=======================================================================\n")

# ----------- Load Data -----------
load(opt$betaPath)
load(opt$mPath)
load(opt$cnPath)

cat("Beta dimensions: ", dim(beta), "\n")
cat("M dimensions: ", dim(m), "\n")
cat("CN dimensions: ", dim(cn), "\n")

pheno <- read.csv(opt$phenoFile)
phenoEWAS <- read.csv(opt$phenoEWAS)

cat("pheno colnames: ", paste(colnames(pheno), collapse = ", "), "\n")
cat("==============\n")
cat("phenoEWAS colnames: ", paste(colnames(phenoEWAS), collapse = ", "), "\n")

cat("=======================================================================\n")

# ----------- Preprocessing Merge (pheno + phenoEWAS) -----------
cat("Merging pheno with EWAS-based annotations...\n")

# Drop specified columns
colsToDrop <- strsplit(opt$dropColumnsPhenoEWAS, ",")[[1]]
phenoEWAS <- phenoEWAS[, !(colnames(phenoEWAS) %in% colsToDrop), drop = FALSE]

# Rename specific columns with suffix
colsToRename <- strsplit(opt$colsToRenamePhenoEWAS, ",")[[1]]
colnames(phenoEWAS)[colnames(phenoEWAS) %in% colsToRename] <- 
        paste0(colnames(phenoEWAS)[colnames(phenoEWAS) %in% colsToRename], 
               ".EWAS")

# Merge using specified mergeKey
mergeKey <- opt$mergeKey
ewasCols <- c(mergeKey, colnames(phenoEWAS)[grepl(".EWAS$", 
                                                  colnames(phenoEWAS))])
pheno <- merge(pheno, phenoEWAS[, ewasCols], by = mergeKey)

cat("Saving merged phenotype file to:", opt$outputPheno, "\n")
write.csv(pheno, file = file.path(opt$outputPheno, "phenoEWAS.csv"), row.names = FALSE)

cat("=======================================================================\n")

# ----------- Convert Variables to Factor with Prefixes -----------
varsToFactor <- strsplit(opt$factorVars, ",")[[1]]
prefixes     <- strsplit(opt$factorPrefixes, ",")[[1]]

for (i in seq_along(varsToFactor)) {
        var <- varsToFactor[i]
        prefix <- ifelse(i <= length(prefixes), prefixes[i], "")
        pheno[[var]] <- factor(paste0(prefix, pheno[[var]]))
}

# ----------- Subsetting Timepoints & Data Splitting -----------
timepoints <- as.numeric(strsplit(opt$timepoints, ",")[[1]])
cat("Subsetting to timepoints:", paste(timepoints, collapse = ", "), "\n")

for (tp in timepoints) {
        assign(paste0("phenoT", tp), subset(pheno, Timepoint == tp))
        assign(paste0("betaT", tp), beta[, grepl(paste0("\\.", tp, "$"), 
                                                 colnames(beta))])
        assign(paste0("mT", tp), m[, grepl(paste0("\\.", tp, "$"), 
                                           colnames(m))])
}

# Save each subset
for (tp in timepoints) {
        write.csv(get(paste0("phenoT", tp)), file = file.path(opt$outputPheno, 
                                                              paste0("phenoT", 
                                                                     tp, ".csv")), 
                  row.names = FALSE)
        save(list = paste0("betaT", tp), file = file.path(opt$outputRData, 
                                                          paste0("betaT", 
                                                                 tp, ".RData")))
        save(list = paste0("mT", tp), file = file.path(opt$outputRData, 
                                                       paste0("mT", 
                                                              tp, ".RData")))
}

# ----------- Merge Combined Timepoints for Longitudinal Analysis -----------
cat("Combining timepoints:", opt$combineTimepoints, "\n")
combineTPs <- as.numeric(strsplit(opt$combineTimepoints, ",")[[1]])

combinedPhenoList <- lapply(combineTPs, function(tp) get(paste0("phenoT", tp)))
phenoCombined <- do.call(rbind, combinedPhenoList)

combineSuffix <- paste0("T", paste(combineTPs, collapse = "T"))
write.csv(phenoCombined, 
          file = file.path(opt$outputPheno, 
                                          paste0("pheno", 
                                                 combineSuffix, ".csv")), 
          row.names = FALSE)

cat("Saved combined phenotype file for T1T2 at:", opt$outputPheno, "\n")
cat("=======================================================================\n")

# ----------- Merge Beta Matrix with Phenotype ----------
mergeBeta <- function(phenoFrame, betaMatrix, id1 = "SID", id2 = "Timepoint") {
        rownames(phenoFrame) <- paste0(phenoFrame[[id1]], 
                                       ".", phenoFrame[[id2]])
        matched <- intersect(rownames(phenoFrame), colnames(betaMatrix))
        phenoFrame <- phenoFrame[matched, ]
        betaMatrix <- betaMatrix[, matched]
        
        betaTranp <- as.data.frame(t(betaMatrix))
        mergedData <- cbind(phenoFrame, betaTranp)
        return(mergedData)
}

# Perform merge for each timepoint

mergedList <- list()
for (tp in timepoints) {
        cat("Processing merge for timepoint:", tp, "\n")
        
        phenoObj <- paste0("phenoT", tp)
        betaObj <- paste0("betaT", tp)
        
        if (!exists(phenoObj) || !exists(betaObj)) {
                cat("Warning: One or both objects not found for T", tp, "\n", sep = "")
                next
        }
        
        phenoTemp <- get(phenoObj)
        betaTemp <- get(betaObj)
        
        cat("  - pheno rows:", nrow(phenoTemp), "\n")
        cat("  - beta cols:", ncol(betaTemp), "\n")
        
        mergedTemp <- tryCatch({
                mergeBeta(phenoTemp, betaTemp)
        }, error = function(e) {
                cat("[ERROR] mergeBeta failed for timepoint", tp, ":\n", conditionMessage(e), "\n")
                return(NULL)
        })
        
        if (!is.null(mergedTemp)) {
                mergedList[[as.character(tp)]] <- mergedTemp
                save(mergedTemp, file = file.path(opt$outputRDataMerge,
                                                  paste0("phenoBetaT", tp, ".RData")))
                cat("Saved merged object for T", tp, "\n", sep = "")
        } else {
                cat("Skipping save for T", tp, " due to error\n")
        }
}


# ----------- Combine merged phenotype + beta matrix -----------
combined <- do.call(rbind, mergedList[as.character(combineTPs)])
save(combined, 
     file = file.path(opt$outputRDataMerge, 
                      paste0("phenoBeta", combineSuffix, ".RData")))

cat("Combined data saved for timepoints:", 
    paste(combineTPs, collapse = ", "), "\n")

cat("Merged data saved to:", opt$outputRDataMerge, "\n")
cat("=======================================================================\n")

# ----------- Close Logging -----------
cat("\nSession Info:\n")
print(sessionInfo())
# =============================================================================
sink(type = "message")
sink()
close(logCon)           
