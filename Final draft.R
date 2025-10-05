# -----------------------------
Load packages
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))
install.packages(c("dplyr","matrixStats"))

library(affy)
library(arrayQualityMetrics)
library(matrixStats)
library(GEOquery)

-----------------------------
2. Read CEL files
-----------------------------

  raw_data <- ReadAffy(celfile.path = "D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_RAW")

list.files("D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_RAW")

library(affy)
library(R.utils)

library(R.utils)


gz_files <- list.files("D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_RAW",
                       pattern = "\\.CEL.gz$", full.names = TRUE)
sapply(gz_files, gunzip, overwrite = TRUE) 

library(affy)

raw_data <- ReadAffy(celfile.path = "D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_RAW")

-----------------------------
  # 3. QC before normalization
  # -----------------------------
library(arrayQualityMetrics)

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "QC_Raw",
                    force = TRUE,
                    do.logtransform = TRUE)


# -----------------------------
# 4. RMA normalization
# -----------------------------
normalized_data <- rma(raw_data)

# -----------------------------
# 5. QC after normalization
# -----------------------------
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "QC_Normalized",
                    force = TRUE)

# -----------------------------
# 6. Filter low-intensity probes
# -----------------------------
expr_matrix <- exprs(normalized_data)
row_median <- rowMedians(expr_matrix)
(threshold <- 3.5)
filtered_data <- expr_matrix[row_median > threshold, ]
dim(filtered_data)  

# -----------------------------
# 7. Load phenotype data
# -----------------------------
gse <- getGEO(filename = "D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_series_matrix.txt")
pheno <- pData(phenoData(gse))

# -----------------------------
# 8. Assign groups (normal vs disease)
# -----------------------------

unique(pheno$source_name_ch1)
groups <- ifelse(grepl(" - E-", pheno$source_name_ch1), "disease", "normal")
groups <- factor(groups, levels = c("normal", "disease"))
table(groups)


table(groups)

nrow(exprs(normalized_data))
dim(filtered_data)
