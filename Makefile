# ===============================================
# Makefile: DNAm Pipeline
# ===============================================
MAKEFLAGS += --output-sync

# Directories
LOGS_DIR = logs
DATA_DIR = data
RDATA_DIR = rData
RESULTS_DIR = results
FIGURES_DIR = figures

# Variables

MODEL ?= firefighter
MODELS = firefighter paramedic emergency

# ----------------------------------------------------
# Run all models 
# ----------------------------------------------------
ifeq ($(MODEL),firefighter)
  PHENOFILE = data/preprocessingMinfiEwasWater/phenoFIRE.csv
else ifeq ($(MODEL),paramedic)
  PHENOFILE = data/preprocessingMinfiEwasWater/phenoMED.csv
else ifeq ($(MODEL),emergency)
  PHENOFILE = data/preprocessingMinfiEwasWater/phenoEMD.csv
else
  PHENOFILE = data/preprocessingMinfiEwasWater/pheno.csv
endif

models:
	@echo "Running models in parallel: $(MODELS)"
	@$(MAKE) -j $(words $(MODELS)) $(addprefix run-, $(MODELS))

$(foreach m,$(MODELS),$(eval run-$(m): 
	@echo ">>> Starting full pipeline for: $(m)"
	$(MAKE) -j1 all MODEL=$(m)
	@echo "<<< Finished full pipeline for: $(m)"
))

default: all

all: \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	data/$(MODEL)/svaEnmix/summary_full_sva2.txt \
	data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1.csv \
	data/$(MODEL)/preprocessingPheno/phenoT2.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1T2.csv \
	data/$(MODEL)/methylationGLM_T1/annotatedGLM.csv \
	data/$(MODEL)/methylationGLMM_T1T2/annotatedLME.csv \
	reports/$(MODEL)/DNAm.pdf
	
# ----------------------------------------------------
# Group target: first3 (Steps 1to3 only)
# ----------------------------------------------------
.PHONY: f3
FIRST3 = \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	data/$(MODEL)/svaEnmix/summary_full_sva2.txt \
	data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1.csv \
	data/$(MODEL)/preprocessingPheno/phenoT2.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1T2.csv

f3: $(FIRST3) reports/$(MODEL)/DNAm.pdf

# ----------------------------------------------------
# Run f3 (Steps 1–3 only) for all models in parallel
# ----------------------------------------------------
f3_models:
	@echo "Running f3 (Steps 1–3) in parallel: $(MODELS)"
	@$(MAKE) -j $(words $(MODELS)) $(addprefix runf3-, $(MODELS))

$(foreach m,$(MODELS),$(eval runf3-$(m): 
	@echo ">>> Starting f3 for: $(m)"
	$(MAKE) -j1 f3 MODEL=$(m)
	@echo "<<< Finished f3 for: $(m)"
)

# ----------------------------------------------------
# Group target: f4 (Steps 1 to 4 only)
# ----------------------------------------------------
.PHONY: f4
FIRST4 = \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	data/$(MODEL)/svaEnmix/summary_full_sva2.txt \
	data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1.csv \
	data/$(MODEL)/preprocessingPheno/phenoT2.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1T2.csv \
	data/$(MODEL)/methylationGLM_T1/annotatedGLM.csv

f4: $(FIRST4) reports/$(MODEL)/DNAm.pdf

# ----------------------------------------------------
# Run f4 (Steps 1–4 only) for all models in parallel
# ----------------------------------------------------
f4_models:
	@echo "Running f4 (Steps 1–4) in parallel: $(MODELS)"
	@$(MAKE) -j $(words $(MODELS)) $(addprefix runf4-, $(MODELS))

$(foreach m,$(MODELS),$(eval runf4-$(m): 
	@echo ">>> Starting f4 for: $(m)"
	$(MAKE) -j1 f4 MODEL=$(m)
	@echo "<<< Finished f4 for: $(m)"
))
# ----------------------------------------------------
# Group target: f3lme (Steps 1 to 3 and LME only)
# ----------------------------------------------------
.PHONY: f3lme
F3LME = \
  rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	data/$(MODEL)/svaEnmix/summary_full_sva2.txt \
	data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1.csv \
	data/$(MODEL)/preprocessingPheno/phenoT2.csv \
	data/$(MODEL)/preprocessingPheno/phenoT1T2.csv \
	data/$(MODEL)/methylationGLMM_T1T2/annotatedLME.csv

f3lme: $(F3LME) reports/$(MODEL)/DNAm.pdf

# ----------------------------------------------------
# Run f3lme (Steps 1–3 + LME) for all models in parallel
# ----------------------------------------------------

f3lme_models:
	@echo "Running f3lme (Steps 1–3 + LME) in parallel: $(MODELS)"
	@$(MAKE) -j $(words $(MODELS)) $(addprefix runf3lme-, $(MODELS))

$(foreach m,$(MODELS),$(eval runf3lme-$(m):
	@echo ">>> Starting f3lme for: $(m)"
	$(MAKE) -j1 f3lme MODEL=$(m)
	@echo "<<< Finished f3lme for: $(m)"
))
# ----------------------------------------------------
# Step 1: Minfi Preprocessing
# ----------------------------------------------------
rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv: preprocessingMinfiEwasWater.R
	Rscript preprocessingMinfiEwasWater.R \
	  --phenoFile $(PHENOFILE) \
	  --outputLogs logs/$(MODEL) \
	  --SampleID SID \
	  --scriptLabel $(MODEL)/preprocessingMinfiEwasWater \
	  --pvalThreshold 0.01 \
	  --mafThreshold 0.1 \
	  --plotGroupVar Timepoint \
	  --lcRef salivaEPIC \
	  --sexColumn Sex \
	  --phenoOrder "SID;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position" \
	  --lcPhenoDir data/$(MODEL)/preprocessingMinfiEwasWater

# ----------------------------------------------------
# Step 2: Surrogate Variable Analysis
# ----------------------------------------------------
data/$(MODEL)/svaEnmix/summary_full_sva2.txt: svaEnmix.R data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv
	Rscript svaEnmix.R \
	  --phenoFile data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv \
	  --outputLogs logs/$(MODEL) \
	  --scriptLabel $(MODEL)/svaEnmix \
	  --SampleID SID \
	  --SentrixIDColumn SentrixBarcode \
	  --SentrixPositionColumn Chipposition

# ----------------------------------------------------
# Step 3: Merge Phenotype
# ----------------------------------------------------
data/$(MODEL)/preprocessingPheno/phenoT1.csv \
data/$(MODEL)/preprocessingPheno/phenoT2.csv \
data/$(MODEL)/preprocessingPheno/phenoT1T2.csv: preprocessingPheno.R data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv
	Rscript preprocessingPheno.R \
	  --phenoFile data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv \
	  --betaPath rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	  --mPath rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	  --cnPath rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData \
	  --SampleID SID \
	  --timepoints 1,2 \
	  --combineTimepoints 1,2 \
	  --outputPheno data/$(MODEL)/preprocessingPheno \
	  --outputRData rData/$(MODEL)/preprocessingPheno/metrics \
	  --outputRDataMerge rData/$(MODEL)/preprocessingPheno/mergeData \
	  --outputLogs logs/$(MODEL) \
	  --sexColumn Sex \
	  --outputDir data/$(MODEL)/preprocessingPheno   
	  
# ----------------------------------------------------
# Step 4: GLM for T1
# ----------------------------------------------------
data/$(MODEL)/methylationGLM_T1/annotatedGLM.csv: methylationGLM_T1.R rData/$(MODEL)/preprocessingPheno/mergeData/phenoBetaT1.RData
	Rscript methylationGLM_T1.R \
	  --inputPheno rData/$(MODEL)/preprocessingPheno/mergeData/phenoBetaT1.RData \
	  --outputLogs logs/$(MODEL) \
	  --outputRData rData/$(MODEL)/methylationGLM_T1/models/ \
	  --outputPlots figures/$(MODEL)/methylationGLM_T1/ \
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
	  --significantCpGDir preliminaryResults/$(MODEL)/cpgs/methylationGLM_T1 \
	  --significantCpGPval 0.00001 \
	  --saveTxtSummaries \
	  --summaryTxtDir preliminaryResults/$(MODEL)/summary/methylationGLM_T1/glm \
	  --fdrThreshold  0.05 \
	  --annotationPackage IlluminaHumanMethylationEPICv2anno.20a1.hg38 \
	  --annotationCols Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group \
	  --annotatedGLMOut data/$(MODEL)/methylationGLM_T1

# ----------------------------------------------------
# Step 5: LME for T1 vs T2 (Longitudinal Analysis)
# ----------------------------------------------------
data/$(MODEL)/methylationGLMM_T1T2/annotatedLME.csv: methylationGLMM_T1T2.R rData/$(MODEL)/preprocessingPheno/mergeData/phenoBetaT1T2.RData
	Rscript methylationGLMM_T1T2.R \
	  --inputPheno rData/$(MODEL)/preprocessingPheno/mergeData/phenoBetaT1T2.RData \
	  --outputLogs logs/$(MODEL) \
	  --outputRData rData/$(MODEL)/methylationGLMM_T1T2/models \
	  --outputPlots figures/$(MODEL)/methylationGLMM_T1T2 \
	  --personVar person \
	  --timeVar Timepoint \
	  --phenotypes Timepoint \
	  --covariates Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells,BMI \
	  --factorVars Sex,Ethnicity,TraumaDefinition,Timepoint \
	  --lmeLibs lme4,lmerTest \
	  --libPath ~/R/x86_64-pc-linux-gnu-library/4.4 \
	  --cpgPrefix cg \
	  --cpgLimit NA \
	  --nCores 64 \
	  --saveSignificantInteractions \
	  --significantInteractionDir preliminaryResults/$(MODEL)/cpgs/methylationGLMM_T1T2 \
	  --significantInteractionPval 0.00001 \
	  --saveTxtSummaries \
	  --summaryTxtDir preliminaryResults/$(MODEL)/summary/methylationGLMM_T1T2 \
	  --fdrThreshold  0.05 \
	  --annotationPackage IlluminaHumanMethylationEPICv2anno.20a1.hg38 \
	  --annotationCols Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group \
	  --annotatedLMEOut data/$(MODEL)/methylationGLMM_T1T2

# ----------------------------------------------------
# Step 6: Final Report
# ----------------------------------------------------
REPORT_INPUTS = DNAm.Rmd \
                rData/$(MODEL)/preprocessingPheno/mergeData/phenoBetaT1.RData \
                rData/$(MODEL)/preprocessingPheno/mergeData/phenoBetaT1T2.RData \
                data/$(MODEL)/preprocessingPheno/phenoT1T2.csv data/$(MODEL)/preprocessingMinfiEwasWater/phenoLC.csv

reports/$(MODEL)/DNAm.pdf: $(REPORT_INPUTS)
	mkdir -p reports/$(MODEL) logs/$(MODEL)
	Rscript -e "rmarkdown::render('DNAm.Rmd', output_file='reports/$(MODEL)/DNAm.pdf')" > logs/$(MODEL)/report.log 2>&1
	@echo "Report built: reports/$(MODEL)/DNAm.pdf (see logs/$(MODEL)/report.log for details)"

# ----------------------------------------------------
# Clean up outputs
# ----------------------------------------------------
clean:
ifeq ($(MODEL),all)
	@echo "===== Cleaning all models ====="
	@for m in $(MODELS); do \
	  echo "Cleaning model $$m ..."; \
	  rm -rf \
	    data/$$m \
	    rData/$$m \
	    reports/$$m \
	    logs/$$m \
	    figures/$$m \
	    results/$$m \
	    preliminaryResults/$$m; \
	done
	@echo "Clean completed for all models. Kept data/preprocessingMinfiEwasWater/."
else
	@echo "Cleaning outputs for: $(MODEL)"
	rm -rf \
		data/$(MODEL) \
		rData/$(MODEL) \
		reports/$(MODEL) \
		logs/$(MODEL) \
		figures/$(MODEL) \
		results/$(MODEL) \
		preliminaryResults/$(MODEL)
	@echo "Clean completed for: $(MODEL)"
endif

# ----------------------------------------------------
# Status Check Target
# ----------------------------------------------------
status:
	@echo "===== Pipeline Status ====="
	@test -e rData/$(MODEL)/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData && echo "? Step 1: preprocessingMinfiEwasWater done" || echo "? Step 1: preprocessingMinfiEwasWater outcome file missing"
	@test -e data/$(MODEL)/svaEnmix/summary_full_sva2.txt && echo "? Step 2: SVA done" || echo "? Step 2: SVA outcome file missing"
	@test -e data/$(MODEL)/preprocessingPheno/phenoT1T2.csv && echo "? Step 3: preprocessingPheno done" || echo "? Step 3: preprocessingPheno outcome file missing"
	@test -e data/$(MODEL)/methylationGLM_T1/annotatedGLM.csv && echo "? Step 4: methylationGLM_T1 done" || echo "? Step 4: methylationGLM_T1 outcome file missing"
	@test -e data/$(MODEL)/methylationGLMM_T1T2/annotatedLME.csv && echo "? Step 5: methylationGLMM_T1T2 done" || echo "? Step 5: methylationGLMM_T1T2 outcome file missing"
	@test -e reports/$(MODEL)/DNAm.pdf && echo "? Step 6: Report generated" || echo "? Step 6: Report missing"
	@echo "============================"

.PHONY: all clean status models f3 f4 f3lme f3_models f4_models f3lme_models
