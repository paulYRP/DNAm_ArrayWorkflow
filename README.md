# EMERGES_DNAm

--------------
# DNAm EPIC Array Analysis Workflow

This repository contains a reproducible pipeline for DNA methylation (DNAm) analysis using Illumina MethylationEPIC v2.0 arrays. It includes preprocessing, quality control, phenotype merging, and statistical modeling (GLM and GLMM) for mental health phenotypes using both `minfi` and `ewastools` frameworks. The pipeline is modular, HPC-compatible, and fully parameterized via command-line arguments.

--------------
## Hardware Requirements

- HPC or workstation with sufficient RAM for large array datasets (≥ 5TB recommended)
- For visualization, interactive usage of RStudio (desktop/server) is recommended

## Software Requirements

- OS: Linux (tested on Rocky Linux 9.2 via Aqua HPC)
- R version: ≥4.3.0 (tested on 4.4.1)
- R packages (CRAN + Bioconductor):

--------------
## Project Structure
The repository is organized to facilitate reproducible analysis, modular development, and HPC compatibility. Below is an overview of the directory layout and the role of each file and folder.

```bash
├── data/
│   ├── preprocessingMinfi/
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
├── preprocessingMinfi.R                    # Minfi-based preprocessing pipeline
├── preprocessingEwastools.R                # Ewastools-based preprocessing pipeline
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

**Run Locally via Makefile**
Use GNU Make to control execution of individual scripts or entire workflows.

```bash
make all
make status
make clean
```

**Submit via PBS (on Aqua HPC)**
To run the full pipeline in a high-performance environment like Aqua (QUT). 

```bash
qsub pipeline.pbs
```

--------------
## References

```r
Bates, D., Mächler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects Models Usinglme4. Journal of Statistical Software, 67(1). https://doi.org/10.18637/jss.v067.i01 
Fisher, J. (2017). BeadSorted.Saliva.EPIC. In https://bioconductor.org/packages/release/data/experiment/vignettes/BeadSorted.Saliva.EPIC/inst/doc/BeadSorted.Saliva.EPIC.html
Hansen, K. D., & Fortin, J.-P. (2025). The minfi User’s Guide. In https://bioconductor.org/packages/devel/bioc/vignettes/minfi/inst/doc/minfi.html
Heiss, J. (2013). Recommended Work Flow. In https://hhhh5.github.io/ewastools/articles/exemplary_ewas.html
Horvath, S. (2013). DNA methylation age of human tissues and cell types. Genome Biol, 14(10), R115. https://doi.org/10.1186/gb-2013-14-10-r115 
Maksimovic, J., Phipson, B., & Oshlack, A. (2017). A cross-package Bioconductor workflow for analysing methylation array data. F1000Research, 5. https://f1000research.com/articles/5-1281 
Marschner, I. C. (2011). glm2: Fitting Generalized Linear Models with Convergence Problems. R Journal, 3(2), 12-15. https://journal.r-project.org/archive/2011/RJ-2011-012/RJ-2011-012.pdf 
Pelegri , D., & Gonzalez, J. R. (2015). Chronological and gestational DNAm age estimation using different methylation-based clocks. In https://bioconductor.org/packages/release/bioc/vignettes/methylclock/inst/doc/methylclock.html
Smialowska, A., Vasquez, L., Ringnér, M., Elsässer, S., Navarro Luzón, C., Nordlund, J., Metzger, A., Contreras‐López, O., Westholm, J., Van Hoef, V., Dethlefsen, O., & Ewels, P. (2025). DNA Methylation: Array Workflow. In https://nbis-workshop-epigenomics.readthedocs.io/en/latest/content/tutorials/methylationArray/Array_Tutorial.html#gene-ontology-testing
```



