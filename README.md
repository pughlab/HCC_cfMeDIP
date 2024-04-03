# Plasma cell-free DNA methylomes for hepatocellular carcinoma (HCC) detection and monitoring of recurrence after liver resection or transplantation

## Introduction
The detection of hepatocellular carcinoma (HCC) via mutational profiling of cell-free DNA (cfDNA) in blood plasma is hampered by mutational diversity and the challenge of acquiring tumor tissue for targeted gene panels. Alternatively, DNA methylation patterns can illuminate biological processes without requiring prior mutation-specific knowledge. Our investigation evaluated the effectiveness of cell-free methylated DNA immunoprecipitation and sequencing (cfMeDIP-Seq) in identifying HCC and monitoring its recurrence after surgery.
This research represents the first application of cfMeDIP-Seq to HCC for both its identification and the tracking of residual disease post-operatively. We collected baseline plasma samples from HCC patients undergoing curative interventions or liver transplants, as well as from healthy individuals, including living liver donors. Utilizing cfMeDIP-Seq and machine learning, we identified HCC-specific differentially methylated regions (DMRs) and developed a HCC methylation score (HMS). This score accurately distinguished HCC with high specificity and sensitivity, and importantly, predicted the likelihood of postoperative recurrence. Analysis of post-surgery follow-up samples using cfMeDIP-Seq revealed that changes in HMS were indicative of future recurrence, demonstrating cfMeDIP-Seq’s utility for both precise HCC detection and recurrence prediction without the need for tumor-specific mutations.  
Our findings demonstrate the potential of cfMeDIP-Seq for sensitive and specific HCC detection irrespective of tumor type, as well as for predicting recurrent disease post-surgery. The comprehensive analysis of tumor-agnostic cfDNA methylomes accurately detected HCC and forecasted recurrence following liver resection or transplantation. This innovative approach may address current challenges in minimally invasive HCC diagnosis using liquid biopsy, offering significant implications for screening, early detection, treatment decisions, and disease progression monitoring and therapy response assessment.

## Contents
The repository contains bash scripts, data processing scripts, and R scripts designed to reproduce figures for the manuscript.

## Setup
R scripts were executed in either R version 3.5.0 or 4.0.1, depending on the specific package dependencies (listeded in Table 1), utilizing RStudio versions (2022.12.0+353 ~ 2023.12.0+369).  

The R packages required are detailed within each script.  

Bash scripts were run on the high performance computing (HPC) cluster at the Princess Margaret Genomics Centre (University Health Network).

## Workflow
![Graphical_abstract](https://github.com/pughlab/HCC_cfMeDIP/assets/109993615/91b31a5c-1920-4214-99c9-5d5c28981fb4)

## Data Analysis
### List of Tools
|Tool         | Function                                   | Version    | Running_Platform  | Language | Alternative_tool| Link                                                             |
|--------------:|:---------------------------------------------------------|:------------|:--------------------------|:------:|:----------------------------|:------------------------------------------------------------------|
| FastQC       | FASTQ QC                                   | 0.11.5     | H4H_shell        | Java     |                            | http://www.bioinformatics.babraham.ac.uk/projects/fastqc        |
| MultiQC      | ❶ FASTQ QC<br>❷ count reads                       | 1.7        | H4H_shell        | Python   |                            | https://multiqc.info                                             |
| fastp        | ❶ IFSD QC3<br>❷ trim                              | 0.23.1     | H4H_shell        | C++      |❶ Trim Galore<br>❷ Trimmomatic<br>❸ Cutadapt| https://github.com/OpenGene/fastp                  |
| bowtie2      | align                                      | 2.4.5      | H4H_shell        | C++      |BWA-mem                     | https://bowtie-bio.sourceforge.net/bowtie2/index.shtml          |
| qualimap     | BAM QC                                     | 2.2        | H4H_shell        | Java     |                            | http://qualimap.conesalab.org                                   |
| SAMtools     | sort - index                               | 1.14       | H4H_shell        | C        |Picard                      | http://www.htslib.org                                           |
| sambamba     | deduplicate                                | 0.7.0      | H4H_shell        | D        |❶ SAMtools<br>❷ Picard| https://lomereiter.github.io/sambamba                           |
| MEDIPS       | ❶ QC2<br>❷ export wig<br>❸ DMRs/ROIs(edgeR)</br>| 1.50.0     | H4H_R/3.5.0      | R        | DESeq2                      | https://doi.org/doi:10.18129/B9.bioc.MEDIPS                     |
| MeDEStrand   |absolute_methylation_level| 0.0.0.9000 |H4H_R/3.5.0      | R        |   | https://github.com/jxu1234/MeDEStrand                           |
| sva/ComBat_seq | reduce batch effect                       | 3.46.0     |H4H_R/4.0.1      | R        |                            | https://github.com/zhangyuqing/ComBat-seq                       |
| edgeR        |❶ TMM normalization<br>❷ CPM counts| 3.28.0     |H4H_R/4.0.1      | R        |                            | https://bioconductor.org/packages/release/bioc/html/edgeR.html  |
| limma        |❶ vCounts<br>❷ DMRs(CPM-trend)| 3.42.0 |❶ H4H_R/4.0.1<br>❷ RStudio | R  |                            | https://bioconductor.org/packages/release/bioc/html/limma.html  |
| FactoMineR/factoextra |PCA                              | 2.8/1.0.7  |❶ H4H_R/4.0.1<br>❷ RStudio | R  |                            | https://rpkgs.datanovia.com/factoextra/index.html               |
| randomForest |❶ build & train models<br>❷ classification prediction| 4.7.1.1 |❶ H4H_R/3.5.0<br>❷ RStudio | R  |                            | https://www.stat.berkeley.edu/users/breiman/RandomForests       |
| glmnet       |regularized linear modeling| 4.1.7      |❶ H4H_R/3.5.0<br>❷ RStudio | R  |                            | https://glmnet.stanford.edu/index.html                          |
| caret        |model training configuration| 6.0        |❶ H4H_R/3.5.0<br>❷ RStudio | R  |                            | https://topepo.github.io/caret                                  |
| pROC         |AUROC                                      | 1.18.0     | RStudio           | R        |                            | https://xrobin.github.io/pROC                                   |
| ggplot2      |                                            | 3.4.2      | RStudio           | R        |                            | https://ggplot2.tidyverse.org                                   |


❶ ❷ ❸ ❹

### H4H - R Compatibility Chart

|      Package | Version    | H4H/R3.5 | H4H/R4.0 | H4H/R4.2 | --mem  | Running Time |
| -----------: | :--------- | :------: | :------: | :------: | ------ | ------------ |
|       MEDIPS | 1.50.0     |   ✅    |    ✅     |  ❌  | ≥ 180G | ~1 day       |
|   MeDEStrand | 0.0.0.9000 |    ✅    |  ❌  |  ❌  | ≥ 180G | 1~2 days     |
|        ChAMP | 2.28.0     |  ❌  |    ✅     |    NA    | ≥ 300G | 1~2 days     |
|       DESeq2 | 1.38.2     |    ✅     |  ❌  |    ✅     | ≥ 180G | ~1 day       |
|        limma | 3.42.0     |    ✅     |    ✅     |    ✅     |        | 1~2 days     |
|        edgeR | 3.28.0     |    ✅     |    ✅    |  ❌  |        | 1~2 days     |
|   FactoMineR | 2.8        |  ❌  |    ✅    |    ✅    |        |              |
|   factoextra | 1.0.7      |  ❌  |    ✅     |    ✅     |        |              |
|        dplyr | 1.1.1      |    ✅    |  ❌  |    ✅     |        |              |
|        tidyr | 1.3.0      |    ✅     |  ❌  |    ✅     |        |              |
|          sva | 3.46.0     |  ❌ |    ✅     |  ❌ | ≥ 500G | ~1 day       |
|        caret | 6.0        |    ✅     |  ❌  |  ❌ | ≥ 80G  | ~2 days      |
| randomForest | 4.7.1.1    |    ✅    |    ✅     | ❌  | ≥ 80G  | ~2 days      |
|       glmnet | 4.1.7      |    ✅     |  ❌  |  ❌  | ≥ 80G  | ~2 days      |
