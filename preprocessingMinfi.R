#!/usr/bin/env Rscript
# ==============================================================================
# DNAm Preprocessing Script (Minfi-based)
# Script Name: preprocessingMinfi.R
# Description: Preprocesses DNAm data using Minfi package. 
# ==============================================================================

# Usage Example (Full version):

# Rscript preprocessingMinfi.R \
#   --phenoFile data/preprocessingMinfi/pheno.csv \
#   --idatFolder data/idats/ \
#   --outputLogs logs/preprocessingMinfi/ \
#   --nSamples 100 \
#   --idColumns SID,Timepoint \
#   --arrayType IlluminaHumanMethylationEPICv2 \
#   --annotationVersion 20a1.hg38 \
#   --scriptLabel preprocessingMinfi \
#   --baseDataFolder rData \
#   --qcTiffPath figures/preprocessingMinfi/quality_control_MSet.tiff \
#   --tiffWidth 2000 \
#   --tiffHeight 1000 \
#   --tiffRes 150 \
#   --qcCutoff 10.5 \
#   --detPtype m+u \
#   --densityTiffPath figures/preprocessingMinfi/densityBeta_MSet.tiff \
#   --pdfReportPath reports/qc_report_RGSet.pdf \
#   --funnormSeed 123 \
#   --normMethods "funnorm;quantile" \
#   --pvalThreshold 0.01 \
#   --chrToRemove chrX,chrY \
#   --snpsToRemove SBE,CpG \
#   --mafThreshold 0.5 \
#   --crossReactivePath data/preprocessingMinfi/12864_2024_10027_MOESM8_ESM.csv \
#   --plotGroupVar Ethnicity \
#   --betaMPlotPath figures/preprocessingMinfi/densityBetaM_MSetF_Flt_Rxy_Ds_Rc.tiff
# ==============================================================================

# Usage Example (Default values with key parameters):

# Rscript preprocessingMinfi.R \
#   --phenoFile data/preprocessingMinfi/pheno.csv \
#   --idatFolder data/idats/ \
#   --nSamples 10 \ ## Testing with 10, NA for all
#   --pvalThreshold 0.01 \
#   --mafThreshold 0.5 \
#   --crossReactivePath data/preprocessingMinfi/12864_2024_10027_MOESM8_ESM.csv \
# ==============================================================================

# ----------- Libraries -----------
suppressPackageStartupMessages({
        library(minfi)
        library(IlluminaHumanMethylationEPICv2manifest)
        library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
        library(RColorBrewer)
        library(ggplot2)
        library(data.table)
        library(optparse)
})
# ==============================================================================
# DNAm Preprocessing Script — Input Arguments

#   --phenoFile            [FILE]   Path to phenotype/sample sheet (.csv or .tsv)
#   --idatFolder           [DIR]    Path to raw .idat files
#   --outputLogs           [DIR]    Directory to save logs and execution metadata
#   --nSamples             [INT]    Number of samples to use (for testing); NA means all
#   --idColumns            [STR]    Comma-separated column names to build sample IDs (e.g., "SID,Timepoint")
#   --arrayType            [STR]    Array type for annotation() (e.g., "IlluminaHumanMethylationEPICv2")
#   --annotationVersion    [STR]    Annotation string for array (e.g., "20a1.hg38")
#   --scriptLabel          [STR]    Script label for tagging folders/outputs (e.g., "preprocessingMinfi")
#   --baseDataFolder       [DIR]    Root folder to store .RData outputs (default: "rData")
#
#   --qcTiffPath           [FILE]   Output TIFF file path for sample QC plot
#   --tiffWidth            [INT]    Width of TIFF images (e.g., 2000)
#   --tiffHeight           [INT]    Height of TIFF images (e.g., 1000)
#   --tiffRes              [INT]    Resolution (dpi) for TIFF images (e.g., 150)
#
#   --qcCutoff             [NUM]    Quality control cutoff value for bad sample detection
#   --detPtype             [STR]    Detection p-value type used ("m+u", "both", etc.)
#   --densityTiffPath      [FILE]   Output TIFF for beta density plot from MSet
#   --pdfReportPath        [FILE]   Path to save PDF report from qcReport()
#
#   --funnormSeed          [INT]    Random seed for normalization
#   --normMethods          [STR]    Normalization method(s), separated by ";" (e.g., "funnorm;swan")
#
#   --pvalThreshold        [NUM]    Detection p-value cutoff for probe filtering (e.g., 0.01)
#   --chrToRemove          [STR]    Comma-separated list of chromosomes to remove (e.g., "chrX,chrY")
#   --snpsToRemove         [STR]    SNP categories to remove (e.g., "SBE,CpG")
#   --mafThreshold         [NUM]    Minor allele frequency threshold (e.g., 0.5)
#
#   --crossReactivePath    [FILE]   Path to file with list of cross-reactive probes to remove
#   --plotGroupVar         [STR]    Phenotype variable used for coloring density plots (e.g., "Ethnicity")
#   --betaMPlotPath        [FILE]   Output TIFF for Beta & M-value density plots after filtering
# ==============================================================================

# ----------- Command Line Arguments -----------
opt <- parse_args(OptionParser(option_list = list(
        make_option("--phenoFile", type = "character", help = "Path to phenotype CSV file", metavar = "FILE"),
        make_option("--idatFolder", type = "character", help = "Folder with IDAT files", metavar = "DIR"),
        make_option("--outputLogs", default = "logs/", help = "Directory for all output", metavar = "DIR"),
        make_option("--nSamples", type = "integer", default = NA, help = "Limit to first N samples [default: all]"),
        make_option("--idColumns", type = "character", default = "SID,Timepoint", help = "Comma-separated ID columns"),
        make_option("--arrayType", default = "IlluminaHumanMethylationEPICv2", help = "Array platform name"),
        make_option("--annotationVersion", default = "20a1.hg38", help = "Annotation version"),
        make_option("--scriptLabel", default = "preprocessingMinfi", help = "Label for output folders/logs"),
        make_option("--baseDataFolder", default = "rData", help = "Base folder for RData output"),
        make_option("--qcTiffPath", default = NULL, help = "Path for QC TIFF output"),
        make_option("--tiffWidth", type = "integer", default = 2000),
        make_option("--tiffHeight", type = "integer", default = 1000),
        make_option("--tiffRes", type = "integer", default = 150),
        make_option("--qcCutoff", type = "double", default = 10.5),
        make_option("--detPtype", default = "m+u", help = "Detection P-value type"),
        make_option("--densityTiffPath", default = NULL, help = "Density plot output TIFF"),
        make_option("--pdfReportPath", default = "reports/qc_report_RGSet.pdf", help = "Path to save QC PDF report"),
        make_option("--funnormSeed", type = "integer", default = 123, help = "Seed for normalization"),
        make_option("--normMethods", default = "funnorm", help = "Normalization methods separated by ; (e.g., funnorm;swan)"),
        make_option("--pvalThreshold", type = "double", default = 0.01),
        make_option("--chrToRemove", default = "chrX,chrY", help = "Chromosomes to remove"),
        make_option("--snpsToRemove", default = "SBE,CpG", help = "SNP positions to filter"),
        make_option("--mafThreshold", type = "double", default = 0.5),
        make_option("--crossReactivePath", type = "character", help = "Path to cross-reactive probe file"),
        make_option("--plotGroupVar", default = "Ethnicity", help = "Grouping variable for density plots"),
        make_option("--betaMPlotPath", default = NULL, help = "Final Beta/M-value plot output TIFF")
        
)))

# Set defaults for paths if NULL
if (is.null(opt$qcTiffPath)) {
        opt$qcTiffPath <- file.path("figures", opt$scriptLabel, 
                                    "quality_control_MSet.tiff")
}
if (is.null(opt$densityTiffPath)) {
        opt$densityTiffPath <- file.path("figures", opt$scriptLabel, 
                                         "densityBeta_MSet.tiff")
}
if (is.null(opt$betaMPlotPath)) {
        opt$betaMPlotPath <- file.path("figures", opt$scriptLabel, 
                                       "densityBetaM_MSetF_Flt_Rxy_Ds_Rc.tiff")
}

# Split comma/semicolon lists
opt$idColList       <- strsplit(opt$idColumns, ",")[[1]]
opt$chrToRemoveList <- strsplit(opt$chrToRemove, ",")[[1]]
opt$snpList         <- strsplit(opt$snpsToRemove, ",")[[1]]
opt$normMethodList  <- strsplit(opt$normMethods, ";")[[1]]
# ==============================================================================

# ----------- Logging Setup -----------
dir.create(opt$outputLogs, recursive = TRUE, showWarnings = FALSE)

logFilePath <- file.path(opt$outputLogs, paste0("log_", opt$scriptLabel, ".txt"))
logCon <- file(logFilePath, open = "wt")  

sink(logCon, split = TRUE)                     
sink(logCon, type = "message")                
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting Preprocessing ====\n")
cat("Start time: ", format(Sys.time()), "\n\n")
cat("Log file path: ", logFilePath, "\n\n")
cat("Pheno file: ", opt$phenoFile, "\n")
cat("IDAT folder: ", opt$idatFolder, "\n")
cat("Log directory: ", opt$outputLogs, "\n")
cat("Sample limit: ", ifelse(is.na(opt$nSamples), "All", opt$nSamples), "\n")
cat("ID columns: ", opt$idColumns, "\n")
cat("Array type: ", opt$arrayType, "\n")
cat("Annotation version: ", opt$annotationVersion, "\n")
cat("Script label: ", opt$scriptLabel, "\n")
cat("Base data folder: ", opt$baseDataFolder, "\n")
cat("Normalization methods: ", paste(opt$normMethodList, collapse = "; "), "\n")
cat("Chromosomes to remove: ", paste(opt$chrToRemoveList, collapse = ", "), "\n")
cat("SNPs to filter: ", paste(opt$snpList, collapse = ", "), "\n")
cat("MAF threshold: ", opt$mafThreshold, "\n")
cat("Pvalue threshold: ", opt$pvalThreshold, "\n")
# =============================================================================

# ----------- Prepare Subfolders -----------
objectDir  <- file.path(opt$baseDataFolder, opt$scriptLabel, "objects")
normDir    <- file.path(opt$baseDataFolder, opt$scriptLabel, "normObjects")
metricsDir <- file.path(opt$baseDataFolder, opt$scriptLabel, "metrics")
filterDir  <- file.path(opt$baseDataFolder, opt$scriptLabel, "filterObjects")

dir.create(objectDir, recursive = TRUE, showWarnings = FALSE)
dir.create(normDir, recursive = TRUE, showWarnings = FALSE)
dir.create(metricsDir, recursive = TRUE, showWarnings = FALSE)
dir.create(filterDir, recursive = TRUE, showWarnings = FALSE)

cat("=======================================================================\n")

# ----------- Read Phenotype File -----------
targets <- read.csv(opt$phenoFile)

if (!is.na(opt$nSamples) && opt$nSamples < nrow(targets)) {
        targets <- targets[1:opt$nSamples, ]
        cat("Subsetting to", opt$nSamples, "samples for testing.\n")
} else {
        cat("Using all", nrow(targets), "samples.\n")
}
cat("=======================================================================\n")

# ----------- Create Sample IDs Based on ID Columns -----------
idColList <- strsplit(opt$idColumns, ",")[[1]]
if (!all(idColList %in% colnames(targets))) {
        stop(paste("Some ID columns not found in phenotype file:", opt$idColumns))
}

targets$sampleID <- do.call(paste, c(targets[idColList], sep = "."))
cat("Sample IDs created using columns:", paste(idColList, collapse = ", "), "\n")
cat("=======================================================================\n")

# ----------- Load IDAT Files into RGSet -----------
rgSet <- read.metharray.exp(
        base = opt$idatFolder,
        targets = targets,
        extended = FALSE,
        recursive = FALSE,
        verbose = FALSE
)

# Assign custom sample names
sampleNames(rgSet) <- targets$sampleID
cat("RGSet loaded with", length(sampleNames(rgSet)), "samples.\n")
cat("=======================================================================\n")

# ----------- Apply Annotation -----------
annotation(rgSet) <- c(
        array = opt$arrayType,
        annotation = opt$annotationVersion
)
cat("Applied annotation: ", paste(annotation(rgSet), collapse = ", "), "\n")
cat("Manifest used:\n")
show(getManifest(rgSet))
cat("=======================================================================\n")

# ----------- Save RGSet -----------
rgSetPath <- file.path(objectDir, "RGSet.RData")
save(rgSet, file = rgSetPath)
cat("RGSet saved to: ", rgSetPath, "\n")
cat("=======================================================================\n")

# ----------- Preprocess Raw (create MSet) -----------
cat("Running preprocessRaw() to generate MSet...\n")
mSet <- preprocessRaw(rgSet)
cat("MSet created with", ncol(mSet), "samples and", nrow(mSet), "probes.\n")
cat("=======================================================================\n")

# Save MSet object
mSetPath <- file.path(objectDir, "MSet.RData")
save(mSet, file = mSetPath)
cat("MSet saved to:", mSetPath, "\n")
cat("=======================================================================\n")

# Display methylated and unmethylated intensity
cat("Preview of methylated intensities:\n")
print(head(getMeth(mSet)[, 1:3]))
cat("=======================================================================\n")
cat("Preview of unmethylated intensities:\n")
print(head(getUnmeth(mSet)[, 1:3]))
cat("=======================================================================\n")

# ----------- Ratio Conversion and Genome Mapping -----------
cat("Converting MSet to RatioSet and GSet...\n")
ratioSet <- ratioConvert(mSet, what = "both", keepCN = TRUE)
cat("RatioSet created.\n")
print(ratioSet)
cat("=======================================================================\n")
gSet <- mapToGenome(ratioSet)
cat("GSet created.\n")
print(gSet)
cat("=======================================================================\n")

# Save RatioSet and GSet
ratioSetPath <- file.path(objectDir, "RatioSet.RData")
gSetPath <- file.path(objectDir, "GSet.RData")
save(ratioSet, file = ratioSetPath)
save(gSet, file = gSetPath)
cat("=======================================================================\n")

# ----------- Extract Methylation Metrics -----------
cat("Extracting methylation level metrics from GSet...\n")

beta <- getBeta(gSet)
cat("Preview of beta values:\n")
print(head(beta[, 1:5]))
cat("=======================================================================\n")
m <- getM(gSet)
cat("Preview of M-values:\n")
print(head(m[, 1:5]))

cat("=======================================================================\n")

cn <- getCN(gSet)
cat("Preview of copy number values:\n")
print(head(cn[, 1:5]))

betaPath <- file.path(metricsDir, "beta_GSet.RData")
mPath <- file.path(metricsDir, "m_GSet.RData")
cnPath <- file.path(metricsDir, "cn_GSet.RData")

save(beta, file = betaPath)
save(m, file = mPath)
save(cn, file = cnPath)
cat("=======================================================================\n")

# ----------- Quality Control Plot (from MSet) -----------
cat("Running QC plotting from MSet object...\n")
qc <- getQC(mSet)

figureDir <- dirname(opt$qcTiffPath)
reportDir <- "reports"
if (!dir.exists(figureDir)) dir.create(figureDir, recursive = TRUE)
if (!dir.exists(reportDir)) dir.create(reportDir, recursive = TRUE)

tiff(filename = opt$qcTiffPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")
plotQC(qc, badSampleCutoff = opt$qcCutoff)
dev.off()

cat("QC plot saved to: ", opt$qcTiffPath, "\n")
cat("=======================================================================\n")

# ----------- Calculate Detection P-values -----------
cat("Calculating detection p-values...\n")

detP <- detectionP(rgSet, type = opt$detPtype)
cat("Detection p-values calculated using type: ", opt$detPtype, "\n")

cat("Preview of detection p-values:\n")
print(head(detP[, 1:5]))

detPpath <- file.path(metricsDir, "detP_RGSet.RData")
save(detP, file = detPpath)
cat("Detection p-values saved to: ", detPpath, "\n")
cat("=======================================================================\n")

# ----------- Density Plot of Beta Values from MSet -----------
cat("Generating density plot of Beta values...\n")

phenoData <- pData(mSet)

# Ensure output directory exists
dir.create(dirname(opt$densityTiffPath), 
           recursive = TRUE, 
           showWarnings = FALSE)

tiff(filename = opt$densityTiffPath,
     width = opt$tiffWidth,
     height = opt$tiffHeight,
     res = opt$tiffRes, type = "cairo")

densityPlot(mSet,
            sampGroups = phenoData[[opt$plotGroupVar]],
            pal = brewer.pal(8, "Dark2"),
            main = paste("Density Plot of Beta Values by", opt$plotGroupVar),
            add = TRUE,
            legend = TRUE)

dev.off()

cat("Density plot saved to: ", opt$densityTiffPath, "\n")
cat("=======================================================================\n")

# ----------- QC Report -----------
cat("Saving QC report PDF to: ", opt$pdfReportPath, "\n")
dir.create(dirname(opt$pdfReportPath), recursive = TRUE, 
           showWarnings = FALSE)

qcReport(rgSet, pdf = opt$pdfReportPath)
cat("=======================================================================\n")

cat("Running normalization methods: ", paste(opt$normMethodList, collapse = ", "), "\n")
normPaths <- c()

firstMethod <- TRUE  # Track whether to assign mSet
for (method in opt$normMethodList) {
        cat("  → Applying normalization:", method, "\n")
        set.seed(opt$funnormSeed)
        
        normObj <- switch(
                method,
                "funnorm"  = preprocessFunnorm(rgSet),
                "illumina" = preprocessIllumina(rgSet),
                "quantile" = preprocessQuantile(rgSet),
                "swan"     = preprocessSWAN(rgSet),
                stop(paste("Unknown normalization method:", method))
        )
        
        if (firstMethod) {
                mSetF <- normObj  # Assign only once
                firstMethod <- FALSE
        }
        
        normPath <- file.path(normDir, paste0("norm_", method, "_MSetF.RData"))
        save(normObj, file = normPath)
        normPaths <- c(normPaths, normPath)
        cat("Saved normalized object: ", normPath, "\n")
}

cat("=======================================================================\n")

# ----------- Probe Filtering Based on Detection P-values -----------
cat("Filtering probes with detection p-values ≥ ", 
    opt$pvalThreshold, "...\n")

# Recompute detection p-values
detP <- detectionP(rgSet)

# Align detection p-values with normalized probes
detP <- detP[match(featureNames(mSetF), rownames(detP)), ]

# Identify probes retained across all samples
keep <- rowSums(detP < opt$pvalThreshold) == ncol(mSetF)
cat("Probes retained: ", sum(keep), "/", length(keep), "\n")

mSetF_Flt <- mSetF[keep, ]
mSetFfltPath <- file.path(filterDir, "removProbes_MSetF_Flt.RData")
save(mSetF_Flt, file = mSetFfltPath)
cat("Filtered object saved to: ", mSetFfltPath, "\n")
cat("=======================================================================\n")

# ----------- Filter Probes on Sex Chromosomes -----------
cat("Removing probes on chromosomes: ", paste(opt$chrToRemoveList, 
                                              collapse = ", "), "\n")
# Identify probes to remove
ann <- getAnnotation(rgSet)
removeProbes <- ann$Name[ann$chr %in% opt$chrToRemoveList]
keepChr <- !(featureNames(mSetF_Flt) %in% removeProbes)

mSetF_Flt_Rxy <- mSetF_Flt[keepChr, ]

cat("Remaining probes after removing selected chromosomes:\n")
print(table(getAnnotation(mSetF_Flt_Rxy)$chr))

rxyPath <- file.path(filterDir, "removChrXY_MSetF_Flt_Rxy.RData")
save(mSetF_Flt_Rxy, file = rxyPath)
cat("Sex chromosome-filtered object saved to: ", rxyPath, "\n")
cat("=======================================================================\n")

# ----------- Remove Probes with SNPs -----------
cat("Removing probes with SNPs at: ", paste(opt$snpList, collapse = ", "), 
    " with MAF >=", opt$mafThreshold, "\n")

# Apply SNP probe filtering
mSetF_Flt_Rxy_Ds <- dropLociWithSnps(
        mSetF_Flt_Rxy,
        snps = opt$snpList,
        maf = opt$mafThreshold
)
cat("Remaining probes after SNP filtering: ", nrow(mSetF_Flt_Rxy_Ds), "\n")

snpPath <- file.path(filterDir, paste0("removSNPs_MAF", opt$mafThreshold, 
                                       "_MSetF_Flt_Rxy_Ds.RData"))
save(mSetF_Flt_Rxy_Ds, file = snpPath)
cat("SNP-filtered object saved to: ", snpPath, "\n")
cat("=======================================================================\n")

# ----------- Remove Cross-Reactive Probes -----------
cat("Loading cross-reactive probe list from:\n", opt$crossReactivePath, "\n")

xReactiveProbes <- read.csv(opt$crossReactivePath, stringsAsFactors = FALSE)

# Filter out cross-reactive probes
keepCr <- !(featureNames(mSetF_Flt_Rxy_Ds) %in% xReactiveProbes$ProbeID)
cat("Probes retained after cross-reactive filter: ", sum(keepCr), "\n")

mSetF_Flt_Rxy_Ds_Rc <- mSetF_Flt_Rxy_Ds[keepCr, ]
rcPath <- file.path(filterDir, "removCrossReactive_MSetF_Flt_Rxy_Ds_Rc.RData")
save(mSetF_Flt_Rxy_Ds_Rc, file = rcPath)
cat("Cross-reactive-filtered object saved to: ", rcPath, "\n")
cat("=======================================================================\n")

# ----------- Final DNAm Matrices from Filtered Data -----------
cat("Extracting final DNAm matrices (M, Beta, CN)...\n")

# M-values
m <- getM(mSetF_Flt_Rxy_Ds_Rc)
mOut <- file.path(metricsDir, "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(m, file = mOut)
cat("M-values saved to: ", mOut, "\n")

# Beta-values
beta <- getBeta(mSetF_Flt_Rxy_Ds_Rc)
betaOut <- file.path(metricsDir, "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(beta, file = betaOut)
cat("Beta-values saved to: ", betaOut, "\n")

# CN-values
cn <- getCN(mSetF_Flt_Rxy_Ds_Rc)
cnOut <- file.path(metricsDir, "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(cn, file = cnOut)
cat("CN-values saved to: ", cnOut, "\n")
cat("=======================================================================\n")

# ----------- Plot Density of Final Beta & M Values by Group Variable -----------
cat("Plotting final density plots for grouping variable: ", 
    opt$plotGroupVar, "\n")

groupFactor <- factor(targets[[opt$plotGroupVar]])

# Create TIFF output
tiff(opt$betaMPlotPath, width = opt$tiffWidth, height = opt$tiffHeight, 
     res = opt$tiffRes, type = "cairo")
par(mfrow = c(1, 2))

# Beta plot
densityPlot(beta,
            sampGroups = groupFactor,
            main = "Beta values",
            legend = FALSE,
            xlab = "Beta values")
legend("top", legend = levels(groupFactor), text.col = brewer.pal(8, "Dark2"))

# M-value plot
densityPlot(m,
            sampGroups = groupFactor,
            main = "M-values",
            legend = FALSE,
            xlab = "M values")
legend("topleft", legend = levels(groupFactor), text.col = brewer.pal(8, 
                                                                      "Dark2"))

dev.off()
cat("Density plots saved to: ", opt$betaMPlotPath, "\n")

cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")  
sink()                  
close(logCon)           

