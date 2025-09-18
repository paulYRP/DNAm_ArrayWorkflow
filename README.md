--------------
# DNAm EPIC Array Analysis Workflow
--------------

This repository contains a reproducible pipeline for DNA methylation (DNAm) analysis using Illumina MethylationEPIC v2.0 arrays. It includes preprocessing, quality control, phenotype merging, and statistical modeling (GLM and GLMM) for mental health phenotypes using both `minfi`, `watermelon` and `ewastools` frameworks. The pipeline is modular, HPC-compatible, and fully parameterized via command-line arguments. 

## News

### 21/08/2025
- It is compatible with all types of tissues using [ewastools](https://hhhh5.github.io/ewastools/articles/exemplary_ewas.html) 
- `ctrlsva` added from [ENmix](https://www.bioconductor.org/packages/devel/bioc/vignettes/ENmix/inst/doc/ENmix.html) 
- `adjusted_funnorm` added from [wateRmelon](https://www.bioconductor.org/packages/release/bioc/vignettes/wateRmelon/inst/doc/wateRmelon.html)
- `make -j1 f3` added to execute only the first three steps. Output:
  - Metrics (beta, M)
  - Objects (RGSet, MSet ...)
  - CSV file with cell estimation
  - Surrogate variable report
  - CSV file and ZIP file compatible with [ClockFundation](https://dnamage.clockfoundation.org/)
- `make -j1 f4` added to execute only the first 4 steps. Output:
  - annotatedGLM.csv
  ```bash
  | IlmnID             | Phenotype1P.Value | Phenotype2P.Value | Phenotype3P.Value | Phenotype4P.Value | Phenotype5P.Value | Phenotype6P.Value | Phenotype7P.Value | Name               | chr   | pos       | UCSC_RefGene_Group                           | UCSC_RefGene_Name         | Relation_to_Island | GencodeV41_Group                     |
  |--------------------|------------------------|----------------------|---------------------|-------------------------|--------------------------|------------------------|--------------------------|--------------------|-------|-----------|-----------------------------------------------|----------------------------|---------------------|--------------------------------------|
  | cgXXXXXXXX_TC21    |                        |                      |                     |                         |                          |                        |                          | cgXXXXXXXX_TC21    | chrX  | ######### | TSS1500;Exon1;5UTR;...                      | RBL2;RBL2;...              | Shore / OpenSea     | exon_1;TSS1500;...                    |
  ``` 
- `make -j1 f3lme` added to execute only the first 3 + LMER. Output:
  - annotatedLME.csv
  ```bash
  | IlmnID             | Phenotype1_Timepoint3_P.Value | Phenotype2_Timepoint3_P.Value | Phenotype3_Timepoint3_P.Value | Phenotype4_Timepoint3_P.Value | Phenotype5_Timepoint3_P.Value | Phenotype6_Timepoint3_P.Value | Phenotype7_Timepoint3_P.Value | Name               | chr   | pos       | UCSC_RefGene_Group                           | UCSC_RefGene_Name         | Relation_to_Island | GencodeV41_Group                     |
  |--------------------|------------------------|----------------------|---------------------|-------------------------|--------------------------|------------------------|--------------------------|--------------------|-------|-----------|-----------------------------------------------|----------------------------|---------------------|--------------------------------------|
  | cgXXXXXXXX_TC21    |                        |                      |                     |                         |                          |                        |                          | cgXXXXXXXX_TC21    | chrX  | ######### | TSS1500;Exon1;5UTR;...                      | RBL2;RBL2;...              | Shore / OpenSea     | exon_1;TSS1500;...                    |
  ```
## Articles/Tutorial:
- [**A Novel Longitudinal Epigenome-Wide Study of Posttraumatic**](https://github.com/n10962646/DNAm_ArrayWorkflow/blob/main/A%20Novel%20Longitudinal%20Epigenome-Wide%20Study%20of%20Posttraumatic.pdf)
- [**DNA Methylation Tutorial**](https://n10962646.github.io/2025CGPHNeurogenomicsWorkshop/tutorial.html)

--------------
## Hardware Requirements

- HPC or workstation with sufficient RAM for large array datasets (≥ 8TB recommended)
- For visualization, interactive usage of RStudio (desktop/server) is recommended

## Software Requirements

- OS: Linux (tested on Rocky Linux 9.2 via Aqua HPC)
- R version: ≥4.4.0 (tested on 4.4.1)
- R packages (CRAN + Bioconductor):

--------------
## Project Structure
The repository is organized to facilitate reproducible analysis, modular development, and HPC compatibility. Below is an overview of the directory layout and the role of each file and folder.

```bash
├── data/
│   ├── preprocessingMinfiEwasWater/
│   │   ├── idats/                          # Folder containing raw IDAT files
│   │   ├── pheno.csv                       # Phenotype file with sample metadata
│   │   └── 12864_2024_10027_MOESM8_ESM.csv # Optional: Cross-reactivity comparison reference
│
├── rData/                                  # Processed R objects (MSet, Beta, CN matrices, etc.)
├── logs/                                   # Logging output from script runs
├── figures/                                # QC plots and visualization output
├── preliminaryResults/                     # Output tables and summary statistics
├── reports/                                # Optional rendered reports (e.g., HTML, PDF)
│
├── preprocessingMinfiEwasWater.R           # Minfi-Ewastool-WateRmelon-based preprocessing pipeline
├── preprocessingPheno.R                    # Phenotype merge and factor conversion
│
├── methylationGLM_T1.R                     # GLM analysis per CpG at T1
├── methylationGLMM_T1T2.R                  # GLMM longitudinal analysis (e.g., T1 vs T2)
├── epigeneticAge_T1T2.R                    # Epigenetic clock analysis
│
├── Makefile                                # Rule-based automation for pipeline steps
└── pipeline.pbs                            # PBS job script for Aqua HPC execution
```
--------------
## Running the Pipeline with Makefile and PBS

This project supports automated execution of all major steps using a Makefile. For high-performance computing (HPC) environments (QUT's Aqua cluster), the pipeline is compatible with PBS Pro job submission.

**Run Locally via Makefile**:
Use GNU Make to control execution of individual scripts or entire workflows.

```bash
make all
make status
make clean
```

**Submit via PBS (on Aqua HPC)**: 
To run the full pipeline in a high-performance environment as Aqua (QUT). 

```bash
qsub pipeline.pbs
```
