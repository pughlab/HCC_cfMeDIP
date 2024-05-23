# Plasma cell-free DNA methylomes for hepatocellular carcinoma (HCC) detection and monitoring of recurrence after liver resection or transplantation

## Introduction
The detection of hepatocellular carcinoma (HCC) via mutational profiling of cell-free DNA (cfDNA) in blood plasma is hampered by mutational diversity and the challenge of acquiring tumor tissue for targeted gene panels. Alternatively, DNA methylation patterns can illuminate biological processes without requiring prior mutation-specific knowledge. Our investigation evaluated the effectiveness of cell-free methylated DNA immunoprecipitation and sequencing (cfMeDIP-Seq) in identifying HCC and monitoring its recurrence after surgery.
This research represents the first application of cfMeDIP-Seq to HCC for both its identification and the tracking of residual disease post-operatively. We collected baseline plasma samples from HCC patients undergoing curative interventions or liver transplants, as well as from healthy individuals, including living liver donors. Utilizing cfMeDIP-Seq and machine learning, we identified HCC-specific differentially methylated regions (DMRs) and developed a HCC methylation score (HMS). This score accurately distinguished HCC with high specificity and sensitivity, and importantly, predicted the likelihood of postoperative recurrence. Analysis of post-surgery follow-up samples using cfMeDIP-Seq revealed that changes in HMS were indicative of future recurrence, demonstrating cfMeDIP-Seq’s utility for both precise HCC detection and recurrence prediction without the need for tumor-specific mutations.  
Our findings demonstrate the potential of cfMeDIP-Seq for sensitive and specific HCC detection irrespective of tumor type, as well as for predicting recurrent disease post-surgery. The comprehensive analysis of tumor-agnostic cfDNA methylomes accurately detected HCC and forecasted recurrence following liver resection or transplantation. This innovative approach may address current challenges in minimally invasive HCC diagnosis using liquid biopsy, offering significant implications for screening, early detection, treatment decisions, and disease progression monitoring and therapy response assessment.

## Contents
The repository contains bash scripts, data processing scripts, and R scripts designed to reproduce figures for the manuscript "Plasma cell-free DNA methylomes for hepatocellular carcinoma detection and monitoring after liver resection or transplantation"

## Setup
R scripts were executed in either R version 3.5.0 or 4.0.1, depending on the specific package dependencies (listeded in Table 1), utilizing RStudio versions (2022.12.0+353 ~ 2023.12.0+369).  
The R packages required are detailed within each script.  
Bash scripts were run on the high performance computing (HPC) cluster at the Princess Margaret Genomics Centre (University Health Network).

## HCC cfMeDIP Study Workflow
![Graphical_abstract](https://github.com/pughlab/HCC_cfMeDIP/assets/109993615/91b31a5c-1920-4214-99c9-5d5c28981fb4)

## HCC cfMeDIP Data Analysis Scheme
![HCC_cfMeDIP_data_analysis_shceme_v2](https://github.com/pughlab/HCC_cfMeDIP/assets/109993615/5f2a9a0a-8a34-48d5-9475-938b8a6c256d)

## List of Dependencies
|Tool         | Function                                   | Version    | Running_Platform  | Language | Alternative_tool| Link                                                             |
|--------------:|:----------------------------------------------------------------------|:------------|:--------------------------|:------:|:----------------------------|:------------------------------------------------------------------|
| FastQC       | FASTQ QC                                   | 0.11.5     | H4H_shell        | Java     |                            | http://www.bioinformatics.babraham.ac.uk/projects/fastqc        |
| MultiQC      | ❶ FASTQ QC<br>❷ count reads                       | 1.7        | H4H_shell        | Python   |                            | https://multiqc.info                                             |
| fastp        | ❶ ISD QC3<br>❷ trim                              | 0.23.1     | H4H_shell        | C++      |❶ Trim Galore<br>❷ Trimmomatic<br>❸ Cutadapt| https://github.com/OpenGene/fastp                  |
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
| ggplot2      |data visualization | 3.4.2      | RStudio           | R        |                            | https://ggplot2.tidyverse.org                                   |

### Scripts List
| Step | Script                       | Analysis                                                                                                           | Panels Generated                                           |
| :--- | :--------------------------- | :----------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------- |
| 1    | step_1_fastq_to_bam.sh       |                                                                                                                    |                                                            |
| 1s   | step_1s_ISD.sh               | QC                                                                                                                 | Figure S2B/2C                                              |
| 2    | step_2_bam_to_wig.R          |                                                                                                                    |                                                            |
| 3    | step_3_wig_to_txt.sh         |                                                                                                                    |                                                            |
| 4    | step_4_matrices_processing.R |                                                                                                                    |                                                            |
| 5    | step_5_medips_qc.R           | QC                                                                                                                 | Figure S2A                                                 |
| 6    | step_6_hcc_classifiy_hms.R   | kFold modeling<br>AUROC<br>HMS<br>PCA<br>Confusion matrix<br>Simply moving average<br>Postoperative HMS trajectory | Figure 2<br>Figure 3<br>Figure 4<br>Figure 5<br>Figure S3C |
| 7    | step_7_hcc_subtyping.R       | HCC subgrouping                                                                                                    | Figure S3A/B                                               |
| 7s   | step_7s_one_vs_each_hcc.R    | OnevsEach.hcc() required by hcc_subgroup.R                                                                         |                                                            |

### Supporting Files List

|  Name                            |  Size      |  Script using this file  |  Description                            |
|:---------------------------------|:-----------|:-------------------------|:----------------------------------------|
|  hg19_chr1_22_coord.rds          |  94.4 MB   |  -                       |                                         |
|  hg19_chr1_22_m_coord.rds        |  94.4 MB   |  step_4_matrices_processing.R   |                                         |
|  hg19_chr1_22_x_y_coord.rds<br>  |  101.4 MB  |  -                       |                                         |
|  hg19_chr1_22_m_x_y_coord.rds    |  101.4 MB  |  -                       |                                         |
|  black_bin_v2.RData              |  7.8 MB    |  step_6_hcc_classifiy_hms.R     |                                         |
|  n236_cpm.rds                    |  1.31 GB   |  step_6_hcc_classifiy_hms.R     |  cfMeDIP signals (summed CPM) analysis  |
|  n236_lcpm.rds                   |  1.34 GB   |  step_6_hcc_classifiy_hms.R     |                                         |
|  vCount_n236.rds              |  1.43 GB   |  step_6_hcc_classifiy_hms.R     |                                         |





### FastQ to BAM Processing Summary
|  File Type                   |  Storage Path                  |  Operation/Function                         |
|:-----------------------------|:-------------------------------|:--------------------------------------------|
|  **Merged FastQ Files**      |  `../bam`                      |  `cat` for merging L001 and L002            |
|  **FastQC Pre-Processing**   |  `../fastqc_pre`               |  `fastqc` for initial QC                    |
|  **Fastp Reports**           |  `../fastqc_pre/fastp_report`  |  `fastp` for quality and adapter trimming   |
|  **FastQC Post-Processing**  |  `../fastqc_post`              |  `fastqc` for QC after `fastp`              |
|  **Aligned SAM Files**       |  Temporary storage             |  `bowtie2` for alignment                    |
|  **Unsorted BAM Files**      |  Temporary storage             |  `samtools view` for SAM to BAM conversion  |
|  **Sorted BAM Files**        |  `../bam`                      |  `samtools sort` for sorting BAM            |
|  **BAM Index Files**         |  `../bam`                      |  `samtools index` for indexing BAM          |
|  **Pre-QC Reports**          |  `../fastq_bam_QC`             |  `samtools flagstat` and `qualimap`         |
|  **Deduplicated BAM Files**  |  `../bam`                      |  `sambamba markdup` for deduplication       |
|  **Post-QC Reports**         |  `../fastq_bam_QC`             |  `samtools flagstat` and `qualimap`         |     


Please note that temporary files, such as unsorted BAM and aligned SAM files, are removed after they are no longer needed to save storage space.



### MEDIPS_QC Results List
This table clearly delineates the output files from the MEDIPS QC analysis, where to find them, and which specific function in the MEDIPS package or custom script is responsible for generating each file. This comprehensive overview aids in understanding the workflow and accessing specific results for further analysis or reporting.

|  Result File Description                                      |  File Pattern/Type                 |  Storage Path                         |  Generating Function                                                                        |
|:--------------------------------------------------------------|:-----------------------------------|:--------------------------------------|:--------------------------------------------------------------------------------------------|
|  Saturation Plots                                             |  `_SaturationPlot.pdf`             |  `outdir/pdf/`                        |  `MEDIPS.plotSaturation()`                                                                  |
|  Sequence Coverage Pie Charts                                 |  `_SeqCoveragePlot_Pie.pdf`        |  `outdir/pdf/`                        |  `MEDIPS.plotSeqCoverage()` (pie)                                                           |
|  Sequence Coverage Histograms                                 |  `_SeqCoveragePlot_Hist.pdf`       |  `outdir/pdf/`                        |  `MEDIPS.plotSeqCoverage()` (hist)                                                          |
|  Sequence Coverage Pie Charts (Including Non-Unique Reads)    |  `_SeqCoveragePlot_Pie_uniqF.pdf`  |  `outdir/pdf/`                        |  `MEDIPS.plotSeqCoverage()` (pie)                                                           |
|  Sequence Coverage Histograms (Including Non-Unique Reads)    |  `_SeqCoveragePlotHist_uniqF.pdf`  |  `outdir/pdf/`                        |  `MEDIPS.plotSeqCoverage()` (hist)                                                          |
|  Saturation Metrics (**FigureS2A, BAM QC**)                   |  `_saturation.txt`                 |  `outdir/QCstats/`                    |  `MEDIPS.saturation()`                                                                      |
|  Coverage Metrics                                             |  `_coverage.txt`                   |  `outdir/QCstats/`                    |  `MEDIPS.seqCoverage()`                                                                     |
|  CpG Enrichment Analysis (**FigureS2A, BAM QC**)              |  `_enrichment.txt`                 |  `outdir/QCstats/`                    |  `MEDIPS.CpGenrich()`                                                                       |
|  Comprehensive QC Statistics                                  |  `_QCstats.txt`                    |  `outdir/QCstats/`                    |  Custom aggregation of MEDIPS data                                                          |
|  Window Coordinates per Chromosome                            |  `_window_per_chr.csv`             |  `outdir/` or specified subdirectory  |  Custom script using MEDIPS data                                                            |
|  PNG versions of Saturation, Coverage, and Calibration Plots  |  Various PNG files                 |  `outdir/png/`                        |  `MEDIPS.plotSaturation()`<br>`MEDIPS.plotSeqCoverage()`<br>`MEDIPS.plotCalibrationPlot()`  |  


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




### Key Arguments

| Tools                | <span style="font-size: 15px;">Arguments</span>                                                                                                                                           | Support Data                                                                                                         | Dependencies                                                                       |
| :------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| fastp                | -w 10<br>-g<br>-f 5<br>-F 5<br>--**adapor_F==**<br>--**adapor_R==**<br>--correction<br>--html=$.html                                                                                      | Adaptor Sequences Table                                                                                              |                                                                                    |
| Bowtie2              | -p 10 <br>**--minins 50**<br>**--maxins 500**<br>-x hg19                                                                                                                                  | 􀄵hg19.1.bt2<br>􀄵hg19.2.bt2<br>􀄵hg19.3.bt2<br>􀄵hg19.4.bt2<br>􀄵hg19.rev.1.bt2<br>􀄵hg19.rev.2.bt2                 |                                                                                    |
| samtools view        | -@ 10 <br>-bS <br>**-f 2**<br>**-F 4**                                                                                                                                                    |                                                                                                                      |                                                                                    |
| qualimap             | -c <br>--java-mem-size=30G <br>-gff hg19.gtf <br>-outformat pdf <br>-nt 10                                                                                                                | hg19.gtf                                                                                                             |                                                                                    |
| sambamba markdup     | -r <br>-p <br>-t 10 <br>**--overflow-list-size 600000**                                                                                                                                   |                                                                                                                      |                                                                                    |
| MEDIPS/QC            | - BSgenome="BSgenome.Hsapiens.UCSC.hg19"<br>- uniq <- 1<br>- **extend <- 300**<br>- shift <- 0<br>-**ws <- 300**<br>- paired <- TRUE<br>- chr.select <- paste0("chr",c(1:22,"X","Y","M")) | NA                                                                                                                   | BSgenome.Hsapiens.UCSC.hg19                                                        |
| MEDIPS/wig_to_counts | - BSgenome="BSgenome.Hsapiens.UCSC.hg19"<br>- **uniq <- 1e-3**<br>- extend <- 300<br>- shift <- 0<br>- ws <- 300<br>- paired <- TRUE<br>- chr.select <- paste0("chr",c(1:22,"M"))         | 􀄵hg19_chr1_22_m_coord.rds<br>􀄵hg19_all_chr_coord.rds<br>􀄵hg19_chr1_22_x_coord.rds<br>􀄵hg19_chr1_22_x_y_coord.rds | 􀄵BSgenome.Hsapiens.UCSC.hg19<br>􀄵edgeR                                           |
| MeDEStrand           | - BSgenome="BSgenome.Hsapiens.UCSC.hg19"<br>- uniq <- 1<br>- extend <- 300<br>- shift <- 0<br>- ws <- 300<br>- paired <- TRUE<br>- chr.select <- paste0("chr",c(1:22,"M"))                |                                                                                                                      | 􀄵BSgenome.Hsapiens.UCSC.hg19<br>􀄵BSgenome<br>􀄵GenomicRanges<br>􀄵MEDIPSData<br> |
