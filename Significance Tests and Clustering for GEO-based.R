#Under the codes we runned the command line in R.

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()


BiocManager::install("GEOquery")  #GSE1987
library(GEOquery)

#question a

lung_c <- getGEO("GSE1987", AnnotGPL = TRUE)
lung_c1 <- lung_c[[1]]
matrix_dataset = exprs(lung_c1)
ExpInfo = lung_c1@phenoData@data
ExpInfo$description 

#We would learn the renal and colon metastasis
#remove the 2 Renal Metastasis 1 Colon metastasis
#[16] "Renal Metastasis. Male." 
#[18] "Colon Metastasis. Female."
#[20] "Renal Metastasis. Male."



#question b

gse_lung <- lung_c1@assayData$exprs
columns_to_remove <- c(2, 4, 6)
gse_lung <- gse_lung[, -columns_to_remove]

#sd = 0 genes, we removed
zero_sd_genes <- which(apply(gse_lung, 1, sd) == 0)
gse_lung_filtered <- gse_lung[-zero_sd_genes, ]

#sd=0 genes, we removed the featureData
feature_lung <- lung_c1@featureData
feature_lung <- feature_lung[-zero_sd_genes, ]

log_gse_lung_filtered <- log2(gse_lung_filtered)
head(log_gse_lung_filtered)

#t-test
p_value = NULL
for (i in 1:nrow(log_gse_lung_filtered)) {
  p_value[i] <- t.test(log_gse_lung_filtered[i, c(1:24, 34)], log_gse_lung_filtered[i, 25:33])$p.value
}

p.val1 <- length(which(p_value<0.01))
p.val1 #output 1534



#question c

#Benjamini-Hochberg correction
corrected_p_values <- p.adjust(p_value, method = "BH")

# Find significantly changed genes using a corrected p-value cutoff of 0.05
significant_genes <- which(corrected_p_values < 0.05)

sig_genes <- length(significant_genes)
sig_genes #1063



#question e
#3 genes-most significant p-values, Find the indices of the smallest three values
most_3 <- order(corrected_p_values)[1:3]
print(most_3) #4346-1183-1922

a <- feature_lung@data$`Gene symbol`
a[most_3] #most common 3 genes; "SPP1"  "SPP1"  "DDX11"

#Select columns for tumor (1:24 and 34) and normal tissue (25:33)

#Correlation test

top3_genes <- log_gse_lung_filtered[most_3, ]

c_top3_genes <- data.frame(t(top3_genes[, c(1:24, 34)]))
colnames(c_top3_genes) <- rownames(top3_genes)
n_top3_genes <- data.frame(t(top3_genes[, 25:33]))
colnames(n_top3_genes) <- rownames(top3_genes)

#Cancer samples correlation tests
p_values_cancer <- matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  for (j in 1:3) {
    p_values_cancer[i, j] <- cor.test(c_top3_genes[, i], c_top3_genes[, j], method = 'pearson')$p.value
  }
}

#Normal samples correlation tests
p_values_normal <- matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  for (j in 1:3) {
    p_values_normal[i, j] <- cor.test(n_top3_genes[, i], n_top3_genes[, j], method = 'pearson')$p.value
  }
}

#Print the results asthe pvalues-cancer and pvalues-normal
print("Cancer Samples Correlation P-Values:")
print(p_values_cancer)
print("Normal Samples Correlation P-Values:")
print(p_values_normal)

#Cancer samples Spearman correlation tests
p_values_cancer_spearman <- matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  for (j in 1:3) {
    p_values_cancer_spearman[i, j] <- cor.test(c_top3_genes[, i], c_top3_genes[, j], method = 'spearman')$p.value
  }
}

#Normal samples Spearman correlation tests
p_values_normal_spearman <- matrix(NA, nrow = 3, ncol = 3)

for (i in 1:3) {
  for (j in 1:3) {
    p_values_normal_spearman[i, j] <- cor.test(n_top3_genes[, i], n_top3_genes[, j], method = 'spearman')$p.value
  }
}

#Print the results for Spearman correlation
print("Cancer Samples Spearman Correlation P-Values:")
print(p_values_cancer_spearman)
print("Normal Samples Spearman Correlation P-Values:")
print(p_values_normal_spearman)



#question f

#Benjamini-Hochberg correction
#corrected_p_values <- p.adjust(p_value, method = "BH")

#Find significantly changed genes using a corrected p-value cutoff of 0.05
significant_genes <- which(corrected_p_values < 0.01)
significant_genes
sig_genes <- length(significant_genes)
sig_genes #533
a[significant_genes]

write.table(a[significant_genes], file = "significant_genes.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



#question g

#Load required libraries
install.packages(c("gplots", "RColorBrewer"))
library(gplots)

# Read the data
data <- log_gse_lung_filtered[significant_genes,]

#Perform hierarchical clustering using Spearman correlation
sg <- as.dist(1 - cor(data, method = "spearman"))
ss <- as.dist(1 - cor(t(data), method = "spearman"))

hc_sg <- hclust(sg, method = "complete")
hc_ss <- hclust(ss, method = "complete")

dend_sg <- as.dendrogram(hc_sg)
dend_ss <- as.dendrogram(hc_ss)

heatmap(data, Colv = dend_sg, Rowv = dend_ss, scale = "row")

#Perform hierarchical clustering using pearson correlation
pg <- as.dist(1 - cor(data, method = "pearson"))
ps <- as.dist(1 - cor(t(data), method = "pearson"))

hc_pg <- hclust(pg, method = "complete")
hc_ps <- hclust(ps, method = "complete")

dend_pg <- as.dendrogram(hc_pg)
dend_ps <- as.dendrogram(hc_ps)

heatmap(data, Colv = dend_pg, Rowv = dend_ps , scale = "row")
              
