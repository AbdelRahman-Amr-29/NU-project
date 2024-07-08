########################### Information about the data
## GSE134878
## Suammary: RNA-seq evaluation of post-mortem human cerebelllum from 33 patients with
#diagnosed Essential tremor, compared to 22 age-matched control patients.
#Two samples were under-sequenced and therefore removed from the final analysis.
#The raw data has been included in this submission.
## DOI: 10.1016/j.neulet.2019.134540

##################################################
## setting working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE134878")

## calling libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(annotables)
library(org.Hs.eg.db)
library(purrr)
library(ggplot2)
library(DESeq2)
library(airway)
library(pheatmap)
library(readxl)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(ggrepel)
library(openxlsx)
library(plotly)

## read the samples
raw_counts <- read_excel("GSE134878_raw_read_counts.xls")

# remove C5 and D12 columns because they are outliers
raw_counts <- raw_counts[, c(-15, -42)]

# make raw_counts data frame
raw_counts <-  data.frame(raw_counts)
class(raw_counts)

##### check the duplicates of gene column
x <- duplicated(raw_counts$gene)
sum(x)

# aggregation
raw_counts <- aggregate(raw_counts, list(raw_counts$gene),FUN=mean)
raw_counts <- raw_counts[, -2]
colnames(raw_counts)[1] <- "gene"

# make the column gene row names 
rownames(raw_counts) <- raw_counts$gene

# remove the gene column
raw_counts <- raw_counts[, -1]

# Remove low counts
raw_counts <- raw_counts[rowSums(raw_counts) >=10, ]

################## Get metadata ##########################
gse <- getGEO(GEO = "GSE134878", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
metadata <- metadata[c(-29, -47),]
head(metadata)

# get needed columns and modify it
colData <- metadata %>%
  dplyr::select(1, 11) %>%
  dplyr::rename(diagnosis = characteristics_ch1.1) %>%
  dplyr::rename(sample_name = title) %>%
  dplyr::mutate(diagnosis = gsub("diagnosis: Control", "control", diagnosis)) %>%
  dplyr::mutate(diagnosis = gsub("diagnosis: Essential Tremor", "ET", diagnosis)) 

head(colData)

# Remove row names from a colData frame
rownames(colData) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData) <- colData$sample_name

# remove the sample_name column from the meatadata_modified
colData <- subset(colData, select = -sample_name)

all(colnames(raw_counts) %in% rownames(colData))
all(rownames(colData) %in% colnames(raw_counts))

# know the right order indexing
data_idx <- match(rownames(colData), colnames(raw_counts))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
raw_counts  <- raw_counts[ , data_idx]

# check if that they are in the same order
all(rownames(colData) == colnames(raw_counts))

# download the modified metadata to use it in the analysis
#write.csv(colData, "colData.csv", row.names = TRUE)

#################### preparing for the DESeq2 analysis

#check the row names in colData match with column names with raw_counts
all(colnames(raw_counts) %in% rownames(colData))

#check the order
all(colnames(raw_counts) == rownames(colData))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_counts, is.numeric)

# Convert numeric columns to integer format
raw_counts[, numeric_columns] <- lapply(raw_counts[, numeric_columns], as.integer)

# Convert 'diagnosis' column to factor
colData$diagnosis <- factor(colData$diagnosis)


########## Data exploration ########
# ######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘A1’
ggplot(raw_counts) +
  geom_histogram(aes(x = A1), stat = "bin", bins = 200) +
  xlab("Raw raw_countsression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘B7’
ggplot(raw_counts) +
  geom_histogram(aes(x = B7), stat = "bin", bins = 200) +
  xlab("Raw raw_countsression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for all the samples
exp.data <- apply(raw_counts, 2, as.numeric)

# Plot histogram of processed RNA-Seq data
hist(log2(exp.data + 1), main = "processed RNA-Seq Histogram")

# Plot boxplot of processed RNA-Seq data
boxplot(log2(exp.data + 1), main = "RNA-Seq Box Plot", col = seq(1:15))

# Mean vs Variance
mean_counts <- apply(raw_counts[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(raw_counts[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")


########### Begin with the DESeq2 analysis 
# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = colData,
                              design = ~ diagnosis)


# set the factor level
dds$diagnosis <-  relevel(dds$diagnosis, ref = "control")

##### Run DESeq
dds <- DESeq(dds)

res <- results(dds, alpha = 0.25)
summary(res)

# Remove rows with null values
res <- res[complete.cases(res), ]

# Convert results to a data frame
res_df <- as.data.frame(res)

# Extract significant genes based on p-value < 0.025 and padj < 0.25
sig_genes <- res[res$pvalue < 0.025 & res$padj < 0.25, ]
summary(sig_genes)
sig_genes <- as.data.frame(sig_genes)

########## Adjusting DEGs For Plotting #############

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Extract the transformed data for significant genes
vst_data <- assay(vsd)

# Add a column to identify significant genes
res_df$significant <- ifelse(res_df$pvalue < 0.025 & res_df$padj < 0.25, 
                             ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"), 
                             "Not Significant")

# Define top genes based on the highest -log10(padj) values
top_genes <- rownames(head(res_df[order(res_df$padj), ], 10))

# Create a new column for labeling top genes
res_df$label <- ifelse(rownames(res_df) %in% top_genes, rownames(res_df), NA)


############### Create the volcano plot ##############
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), size = 1.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue"),
                     name = "Gene Regulation") +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add text labels for the top significant genes
p + geom_text_repel(aes(label = label),
                    size = 3,
                    box.padding = 0.3,
                    point.padding = 0.3,
                    segment.color = 'grey50')


############# 2D PCA ##############
# Perform PCA manually
pca_res <- prcomp(t(vst_data))

# Create a data frame with the PCA results
pca_data <- as.data.frame(pca_res$x)
pca_data$diagnosis <- colData(dds)$diagnosis

# Plot the first two principal components (PC1 vs. PC2)
pca_plot1 <- ggplot(pca_data, aes(PC1, PC2, color = diagnosis)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * pca_res$sdev[1]^2 / sum(pca_res$sdev^2), 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "% variance")) +
  ggtitle("PCA of RNA-seq Data: PC1 vs. PC2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

print(pca_plot1)

# Plot the second and third principal components (PC2 vs. PC3)
pca_plot2 <- ggplot(pca_data, aes(PC2, PC3, color = diagnosis)) +
  geom_point(size = 3) +
  xlab(paste0("PC2: ", round(100 * pca_res$sdev[2]^2 / sum(pca_res$sdev^2), 1), "% variance")) +
  ylab(paste0("PC3: ", round(100 * pca_res$sdev[3]^2 / sum(pca_res$sdev^2), 1), "% variance")) +
  ggtitle("PCA of RNA-seq Data: PC2 vs. PC3") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

print(pca_plot2)

# Scree plot to visualize the explained variance of each principal component
scree_data <- data.frame(
  PC = paste0("PC", 1:length(pca_res$sdev)),
  Variance = round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
)

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity") +
  xlab("Principal Components") +
  ylab("Percentage of Variance Explained") +
  ggtitle("Scree Plot") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

print(scree_plot)


######### 3D PCA ################

# Perform PCA
pca_res <- prcomp(t(vst_data))

# Create a data frame with PCA results and sample information
pca_data <- data.frame(PC1 = pca_res$x[,1], 
                       PC2 = pca_res$x[,2], 
                       PC3 = pca_res$x[,3], 
                       diagnosis = colData(dds)$diagnosis)

# Plot the 3D PCA
fig <- plot_ly(pca_data, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~diagnosis, 
               colors = c("blue", "red"),
               text = ~paste("Sample: ", rownames(pca_data)),
               marker = list(size = 5)) %>%
  layout(title = '3D PCA of RNA-seq Data',
         scene = list(xaxis = list(title = paste0("PC1: ", round(100 * summary(pca_res)$importance[2,1], 1), "% variance")),
                      yaxis = list(title = paste0("PC2: ", round(100 * summary(pca_res)$importance[2,2], 1), "% variance")),
                      zaxis = list(title = paste0("PC3: ", round(100 * summary(pca_res)$importance[2,3], 1), "% variance"))))

# Display the plot
fig

############# Heat Map ##############
# Extract the top 50 DEGs based on adjusted p-values
top50_genes <- head(rownames(res[order(res$padj), ]), 50)

# Subset the transformed data for these top 50 genes
top50_data <- assay(vsd)[top50_genes, ]

# Ensure the data is numeric
top50_data <- as.matrix(top50_data)

# Create annotation for columns
annotation_col <- data.frame(diagnosis = colData(dds)$diagnosis)
rownames(annotation_col) <- colnames(top50_data)

# Define colors for the annotation
ann_colors <- list(
  diagnosis = c("ET" = "red", "control" = "blue")
)

# Create a heatmap with a different color palette and sorted samples
pheatmap(top50_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50),  # New color palette
         main = "Heatmap of Top 50 DEGs",
         fontsize_row = 8, 
         fontsize_col = 10, 
         fontsize = 12,
         legend = TRUE,
         annotation_legend = TRUE,
         border_color = NA)

############### Functional enrichment analysis using clusterprofiler ################
ego_CC <- enrichGO(gene = degs$gene,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = ,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

# ## Output results from CC term
cluster_summary_CC <- data.frame(ego_CC)

#write.csv(ego_CC, "CC_results.csv", row.names = TRUE)

# MF
ego_MF <- enrichGO(gene = degs$gene,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

#save results
cluster_summary_MF <- data.frame(ego_MF)

# BP
ego_BP <- enrichGO(gene = degs$gene,
                   keyType = "SYMBOL",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.1,
                   qvalueCutoff = 0.5,
                   readable = TRUE)

#save results
cluster_summary_BP <- data.frame(ego_BP)

write.csv(ego_BP, "BP_results.csv", row.names = TRUE)


# Convert gene symbols to ENTREZ IDs
gene_info <- bitr(degs$gene, fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = org.Hs.eg.db)
gene_entrez <- gene_info$ENTREZID

# Perform KEGG pathway enrichment analysis
kk <- enrichKEGG(gene         = gene_entrez,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

# Save the results to a CSV file
kk_data_frame <- as.data.frame(kk)
# write.csv(kk_data_frame, file = "KEGG_pathway_results.csv", row.names = FALSE)

dotplot(kk, showCategory=20)

# Visualize a specific pathway
pathview(gene.data  = gene_entrez,
         pathway.id = "hsa04933",
         species    = "hsa")



# versions
packageVersion("DESeq2")
packageVersion("GEOquery")
packageVersion("ggplot2")
packageVersion("clusterProfiler")
