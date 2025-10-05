if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"))
install.packages("dplyr")
install.packages("matrixStats")


library(affy)

raw_data <- ReadAffy(celfile.path = "D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_RAW")


library(arrayQualityMetrics)

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "QC_Raw",
                    force = TRUE,
                    do.logtransform = TRUE)
Normalization

rma(raw_data)


arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "QC_Normalized",
                    force = TRUE)

normalized_data <- rma(raw_data)


library(affy)
normalized_data <- rma(raw_data)


arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "QC_Normalized",
                    force = TRUE)

if (!requireNamespace("arrayQualityMetrics", quietly = TRUE)) {
  BiocManager::install("arrayQualityMetrics")
}
library(arrayQualityMetrics)


arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "QC_Normalized",
                    force = TRUE)



expr_matrix <- exprs(normalized_data)
row_median <- rowMedians(expr_matrix)
threshold <- 3.5
filtered_data <- expr_matrix[row_median > threshold, ]



expr_matrix <- exprs(normalized_data)


library(matrixStats)

row_median <- rowMedians(expr_matrix)


threshold <- 3.5    

filtered_data <- expr_matrix[row_median > threshold, ]

dim(filtered_data)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO(filename = "D:/My Files/AI and Bioinfo/R_Bioinformatics_Class/Data/GSE7305_series_matrix.txt")
pheno <- pData(phenoData(gse))


(N = normal, E = disease)
groups <- ifelse(grepl("^-?E-", pheno$source_name_ch1) | grepl("-E-", pheno$source_name_ch1), "disease", "normal")
groups <- factor(groups, levels = c("normal", "disease"))


table(groups)
