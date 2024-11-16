library(DESeq2)
library(tidyverse)
library(ggplot)
library(ggplot2)
library(plotMA)
library(dplyr)

count_data_2001 <- read.csv("H0351.2001/RNAseqCounts.csv", header = TRUE, row.names = 1)
count_data_2002 <- read.csv("H0351.2002/RNAseqCounts.csv", header = TRUE, row.names = 1)
metadata_2001 <- read.csv("H0351.2001/SampleAnnot.csv", header = TRUE)
metadata_2002 <- read.csv("H0351.2002/SampleAnnot.csv", header = TRUE)
dds_2001 <- DESeqDataSetFromMatrix(countData = count_data_2001,
colData = metadata_2001,
design = ~ RNAseq_sample_name)
count_data_2001 <- round(count_data_2001)
count_data_2002 <- round(count_data_2002)
dds_2001 <- DESeqDataSetFromMatrix(countData = count_data_2001,
colData = metadata_2001,
design = ~ RNAseq_sample_name)
dds_2002 <- DESeqDataSetFromMatrix(countData = count_data_2002,
colData = metadata_2002,
design = ~ RNAseq_sample_name)
dds_2001 <- DESeq(dds_2001)
str(metadata_2001)
head(metadata_2001)
# Create a design formula based on the main_structure variable
design_formula <- ~ main_structure
# Create DESeq2 object
dds_2001 <- DESeqDataSetFromMatrix(countData = count_data_2001,
colData = metadata_2001,
design = design_formula)
dds_2001 <- DESeq(dds_2001)
res_2001 <- results(dds_2001)
head(res_2001)
write.csv(as.data.frame(res_2001), file = "DESeq2_results_2001.csv", row.names = TRUE)
dds_2002 <- DESeqDataSetFromMatrix(countData = count_data_2002,
colData = metadata_2002,
design = design_formula)
dds_2002 <- DESeq(dds_2002)
res_2002 <- results(dds_2002)
head(res_2002)
write.csv(as.data.frame(res_2002), file = "DESeq2_results_2002.csv", row.names = TRUE)
plotMA(res, main = "DESeq2", ylim = c(-5, 5))
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
geom_point() +
theme_minimal() +
labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
install.packages("ggplot2")
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
geom_point() +
theme_minimal() +
labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
geom_point() +
theme_minimal() +
labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
ggplot(res_2001, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
geom_point() +
theme_minimal() +
labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted p-value")
rlang::last_trace()
res_2001_significant <- res_2001 %>% filter(padj < 0.05)
library(DESeq2)
res_2001_df <- as.data.frame(res_2001)
res_2002_df <- as.data.frame(res_2002)
res_2001_significant <- res_2001_df %>% filter(padj < 0.05)
res_2002_significant <- res_2002_df %>% filter(padj < 0.05)
res_2001_significant$source <- "2001"
res_2002_significant$source <- "2002"
combined_significant <- bind_rows(res_2001_significant, res_2002_significant)
ggplot(combined_significant, aes(x = baseMean, y = log2FoldChange, color = source)) +
geom_point(alpha = 0.5) +
scale_x_log10() +
theme_minimal() +
labs(title = "Significant Genes",
x = "Base Mean (log scale)",
y = "Log2 Fold Change",
color = "Source")
res_2001_significant$source <- "H0351.2001"
res_2002_significant$source <- "H0351.2002"
combined_significant <- bind_rows(res_2001_significant, res_2002_significant)
ggplot(combined_significant, aes(x = baseMean, y = log2FoldChange, color = source)) +
geom_point(alpha = 0.5) +
scale_x_log10() +
theme_minimal() +
labs(title = "Significant Genes",
x = "Base Mean (log scale)",
y = "Log2 Fold Change",
color = "Source")

res_2001$gene <- rownames(res_2001)
res_2002$gene <- rownames(res_2002)
res_2001_df <- as.data.frame(res_2001)
res_2002_df <- as.data.frame(res_2002)
upregulated_genes <- res_2001_df %>%
filter(Log2FoldChange > 0)
upregulated_genes <- res_2001_df %>%
filter(log2FoldChange > 0)
# Filter downregulated genes (Log2FoldChange < 0)
downregulated_genes <- res_2001_df %>%
filter(log2FoldChange < 0)
upregulated_file <- 'upregulated_genes.csv'
write_csv(upregulated_genes, upregulated_file)
upregulated_file <- 'upregulated_genes.csv'
write.csv(upregulated_genes, upregulated_file)
# Save the downregulated genes to a new Excel file
downregulated_file <- 'downregulated_genes.csv'
write.csv(downregulated_genes, downregulated_file)
upregulated_genes <- res_2002_df %>%
filter(log2FoldChange > 0)
# Filter downregulated genes (Log2FoldChange < 0)
downregulated_genes <- res_2002_df %>%
filter(log2FoldChange < 0)
upregulated_file <- 'upregulated_genes.csv'
write.csv(upregulated_genes, upregulated_file)
# Save the downregulated genes to a new Excel file
downregulated_file <- 'downregulated_genes.csv'
write.csv(downregulated_genes, downregulated_file)
