# MicroarrayQC_R
Comprehensive microarray analysis of GSE7305 dataset: QC, normalization, and differential expression in endometriosis.

## Overview
This repository contains an R workflow for **quality control, normalization, and differential expression analysis** of the publicly available **GSE7305 microarray dataset**, focusing on **normal vs endometriosis endometrial tissues**.

The workflow includes:  
1. Reading raw CEL files.  
2. Quality control before and after normalization using `arrayQualityMetrics`.  
3. RMA normalization with `affy`.  
4. Filtering low-intensity probes.  
5. Assigning sample groups (normal vs disease).  
6. Preparing data for downstream differential expression analysis.

---

## Repository Structure
MicroarrayQC_R/
│
├── scripts/ # R scripts for QC, normalization, and analysis
│ ├── 01_read_data.R
│ ├── 02_QC_pre_norm.R
│ ├── 03_RMA_normalization.R
│ ├── 04_QC_post_norm.R
│ ├── 05_filter_low_intensity.R
│ └── 06_define_groups.R
│
├── QC_Raw/ # QC reports before normalization
├── QC_Normalized/ # QC reports after normalization
├── README.md
