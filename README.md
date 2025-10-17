--------------
# DNAm EPIC Array Analysis Workflow
--------------

This repository contains a reproducible pipeline for DNA methylation (DNAm) analysis using Illumina MethylationEPIC v2.0 arrays. It includes preprocessing, quality control, phenotype merging, and statistical modeling (GLM and GLMM) for mental health phenotypes using both `minfi`, `watermelon` and `ewastools` frameworks. The pipeline is modular, HPC-compatible, and fully parameterized via command-line arguments. 

--------------
## Articles:
- [**A Pilot Epigenome-Wide Study of Posttraumatic Growth: Identifying Novel Candidates for Future Research**](https://www.mdpi.com/2075-4655/9/4/39)

--------------
## Tutorials:
- [**DNA Methylation Tutorial**](https://n10962646.github.io/2025CGPHNeurogenomicsWorkshop/tutorial.html)
- [**Getting Started**](https://github.com/paulYRP/DNAm_ArrayWorkflow/wiki/Getting-Started)
- [**Requirements**](https://github.com/paulYRP/DNAm_ArrayWorkflow/wiki/Requirements)

--------------
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
- `make -j1 f3lme` added to execute only the first 3 + LMER. Output:
  - annotatedLME.csv
- `DNA.pdf`, automic report generated with a summary of all steps.
- Interaction added to `methylationGLM_T1`.
- `methylationGLM_T1` and `methylationGLMM_T1T2` extract categorical and numerical coefficients to the annotatedGLM.csv and annotatedLME.csv. 
