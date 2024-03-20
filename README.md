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
