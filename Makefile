# ===============================================
# Makefile: DNAm Pipeline
# ===============================================

# Directories
LOGS_DIR = logs
DATA_DIR = data
RDATA_DIR = rData
RESULTS_DIR = results
FIGURES_DIR = figures

# Variables
MODEL = profession

# Default target
all: \
	rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	data/svaEnmix/sva/summary_full_sva2.txt \
	data/preprocessingMinfiEwasWater/phenoLC.csv \
	data/preprocessingPheno/phenoT1.csv \
	data/preprocessingPheno/phenoT2.csv \
	data/preprocessingPheno/phenoT1T2.csv \
	data/methylationGLM_T1/$(MODEL)/annotatedGLM.csv \
	data/methylationGLMM_T1T2/$(MODEL)/annotatedLME.csv \
	reports/DNAm.pdf
	
# ----------------------------------------------------
# Group target: first3 (Steps 1to3 only)
# ----------------------------------------------------
.PHONY: f3
FIRST3 = \
  rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  data/preprocessingMinfiEwasWater/phenoLC.csv \
  data/svaEnmix/sva/summary_full_sva2.txt \
  data/preprocessingPheno/phenoT1.csv \
  data/preprocessingPheno/phenoT2.csv \
  data/preprocessingPheno/phenoT1T2.csv

f3: $(FIRST3) reports/DNAm.pdf

# ----------------------------------------------------
# Group target: f4 (Steps 1 to 4 only)
# ----------------------------------------------------
.PHONY: f4
FIRST4 = \
  rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  data/preprocessingMinfiEwasWater/phenoLC.csv \
  data/svaEnmix/sva/summary_full_sva2.txt \
  data/preprocessingPheno/phenoT1.csv \
  data/preprocessingPheno/phenoT2.csv \
  data/preprocessingPheno/phenoT1T2.csv \
  data/methylationGLM_T1/$(MODEL)/annotatedGLM.csv

f4: $(FIRST4) reports/DNAm.pdf

# ----------------------------------------------------
# Group target: f3lme (Steps 1 to 3 and LME only)
# ----------------------------------------------------
.PHONY: f3lme
F3LME = \
  rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
  data/preprocessingMinfiEwasWater/phenoLC.csv \
  data/svaEnmix/sva/summary_full_sva2.txt \
  data/preprocessingPheno/phenoT1.csv \
  data/preprocessingPheno/phenoT2.csv \
  data/preprocessingPheno/phenoT1T2.csv \
  data/methylationGLMM_T1T2/$(MODEL)/annotatedLME.csv

f3lme: $(F3LME) reports/DNAm.pdf

# ----------------------------------------------------
# Step 1: Minfi Preprocessing
# ----------------------------------------------------
rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
data/preprocessingMinfiEwasWater/phenoLC.csv: preprocessingMinfiEwasWater.R
	Rscript preprocessingMinfiEwasWater.R \
	  --SampleID SID \
	  --pvalThreshold 0.01 \
	  --mafThreshold 0.1 \
	  --plotGroupVar Timepoint \
	  --lcRef salivaEPIC \
	  --sexColumn Sex \
	  --phenoOrder "SID;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position"

# ----------------------------------------------------
# Step 2: Surrogate Variable Analysis
# ----------------------------------------------------
data/svaEnmix/sva/summary_full_sva2.txt: svaEnmix.R data/preprocessingMinfiEwasWater/phenoLC.csv
	Rscript svaEnmix.R \
	  --SampleID SID \
	  --SentrixIDColumn SentrixBarcode \
	  --SentrixPositionColumn Chipposition

# ----------------------------------------------------
# Step 3: Merge Phenotype
# ----------------------------------------------------
data/preprocessingPheno/phenoT1.csv \
data/preprocessingPheno/phenoT2.csv \
data/preprocessingPheno/phenoT1T2.csv: preprocessingPheno.R data/preprocessingMinfiEwasWater/phenoLC.csv
	Rscript preprocessingPheno.R \
	  --betaPath rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	  --mPath rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	  --cnPath rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	  --SampleID SID \
	  --timepoints 1,2 \
	  --combineTimepoints 1,2 \
	  --sexColumn Sex   
	  
# ----------------------------------------------------
# Step 4: GLM for T1
# ----------------------------------------------------
data/methylationGLM_T1/$(MODEL)/annotatedGLM.csv: methylationGLM_T1.R rData/preprocessingPheno/mergeData/phenoBetaT1.RData
	Rscript methylationGLM_T1.R \
	  --inputPheno rData/preprocessingPheno/mergeData/phenoBetaT1.RData \
	  --outputLogs logs/methylationGLM_T1/$(MODEL) \
	  --outputRData rData/methylationGLM_T1/models/$(MODEL) \
	  --outputPlots figures/methylationGLM_T1/$(MODEL) \
	  --phenotypes Profession \
	  --covariates Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells,BMI \
	  --factorVars Sex,Ethnicity,TraumaDefinition,Profession \
	  --cpgPrefix cg \
	  --cpgLimit NA \
	  --nCores 64 \
	  --plotWidth 2000 --plotHeight 1000 --plotDPI 150 \
	  --libPath ~/R/x86_64-pc-linux-gnu-library/4.4 \
	  --glmLibs glm2 \
	  --summaryPval NA \
	  --summaryResidualSD \
	  --saveSignificantCpGs \
	  --significantCpGDir preliminaryResults/cpgs/methylationGLM_T1/$(MODEL) \
	  --significantCpGPval 0.00001 \
	  --saveTxtSummaries \
	  --summaryTxtDir preliminaryResults/summary/methylationGLM_T1/glm/$(MODEL) \
	  --fdrThreshold  0.05 \
	  --annotationPackage IlluminaHumanMethylationEPICv2anno.20a1.hg38 \
	  --annotationCols Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group \
	  --annotatedGLMOut data/methylationGLM_T1/$(MODEL)

# ----------------------------------------------------
# Step 5: LME for T1 vs T2 (Longitudinal Analysis)
# ----------------------------------------------------
data/methylationGLMM_T1T2/$(MODEL)/annotatedLME.csv: methylationGLMM_T1T2.R rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData
	Rscript methylationGLMM_T1T2.R \
	  --inputPheno rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData \
	  --outputLogs logs/methylationGLMM_T1T2/$(MODEL) \
	  --outputRData rData/methylationGLMM_T1T2/models/$(MODEL) \
	  --outputPlots figures/methylationGLMM_T1T2/$(MODEL) \
	  --personVar person \
	  --timeVar Timepoint \
	  --phenotypes Profession \
	  --covariates Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells,BMI \
	  --factorVars Sex,Ethnicity,TraumaDefinition,Profession \
	  --lmeLibs lme4,lmerTest \
	  --libPath ~/R/x86_64-pc-linux-gnu-library/4.4 \
	  --cpgPrefix cg \
	  --cpgLimit NA \
	  --nCores 64 \
	  --interactionTerm Timepoint \
	  --saveSignificantInteractions \
	  --significantInteractionDir preliminaryResults/cpgs/methylationGLMM_T1T2/$(MODEL) \
	  --significantInteractionPval 0.00001 \
	  --saveTxtSummaries \
	  --summaryTxtDir preliminaryResults/summary/methylationGLMM_T1T2/$(MODEL) \
	  --fdrThreshold  0.05 \
	  --annotationPackage IlluminaHumanMethylationEPICv2anno.20a1.hg38 \
	  --annotationCols Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group \
	  --annotatedLMEOut data/methylationGLMM_T1T2/$(MODEL)

# ----------------------------------------------------
# Step 6: Final Report
# ----------------------------------------------------
REPORT_INPUTS = DNAm.Rmd \
                rData/preprocessingPheno/mergeData/phenoBetaT1.RData \
                rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData \
                data/preprocessingPheno/phenoT1T2.csv data/preprocessingMinfiEwasWater/phenoLC.csv

reports/DNAm.pdf: $(REPORT_INPUTS)
	mkdir -p reports logs
	Rscript -e "rmarkdown::render('DNAm.Rmd', output_file='reports/DNAm.pdf')" > logs/report.log 2>&1
	@echo "Report built: reports/DNAm.pdf (see logs/report.log for details)"


# ----------------------------------------------------
# Clean up outputs
# ----------------------------------------------------
clean:
	rm -rf data/preprocessingMinfiEwasWater/phenoLC.csv \
	       data/preprocessingPheno/* data/svaEnmix/* \
	       data/methylationGLM_T1/* \
	       data/methylationGLMM_T1T2/* \
=	       results/* figures/* logs/* preliminaryResults/* rData/* reports/*

# ----------------------------------------------------
# Status Check Target
# ----------------------------------------------------
status:
	@echo "===== Pipeline Status ====="
	@test -e rData/preprocessingMinfiEwasWater/objects/RGSet.RData && echo "? Step 1: preprocessingMinfiEwasWater done" || echo "? Step 1: preprocessingMinfiEwasWater outcome file missing"
	@test -e data/svaEnmix/sva/summary_full_sva2.txt && echo "? Step 2: SVA done" || echo "? Step 2: SVA outcome file missing"
	@test -e data/preprocessingPheno/phenoT1T2.csv && echo "? Step 3: preprocessingPheno done" || echo "? Step 3: preprocessingPheno outcome file missing"
	@test -e data/methylationGLM_T1/$(MODEL)/annotatedGLM.csv && echo "? Step 4: methylationGLM_T1 done" || echo "? Step 4: methylationGLM_T1 outcome file missing"
	@test -e data/methylationGLMM_T1T2/$(MODEL)/annotatedLME.csv && echo "? Step 5: methylationGLMM_T1T2 done" || echo "? Step 5: methylationGLMM_T1T2 outcome file missing"
	@test -e reports/DNAm.pdf && echo "? Step 6: Report generated" || echo "? Step 6: Report missing"
	@echo "============================"


.PHONY: all clean
