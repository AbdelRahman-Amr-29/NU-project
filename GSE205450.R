################### Information about the data#######################
# GSE205450
# Summary:    We analyzed transcriptome differences in postmortem caudate and putamen
#from controls (n=40) and PD (n=35) as confirmed by autopsy.
#Further, we analyzed region specific differences in caudate and putamen
#associated with clinical variables.
# NOTE: For differential gene expression analysis between controls and PD,
#all 35 PD and 40 control specimens (caudate and putamen) obtained from NIH NeuroBiobank
#were included and processed for bulk RNA-seq as INDEPENDENT biological replicates.
# DOI: https://doi.org/10.1038/s41467-023-39652-6
#############################################################################

# Setting Working Directory
setwd("D:/Nile University/OneDrive - Nile University/Desktop/PD2")

# calling libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(org.Hs.eg.db)
library(purrr)
library(ggplot2)
library(DESeq2)
library(airway)
library(pheatmap)
library(ggrepel)
library(plotly)
library(pheatmap)

# Read txt file & convert it to table
raw_data <- read.table("D:/Nile University/OneDrive - Nile University/Desktop/DS Data/GSE205450_counts.table.txt", header = TRUE, sep = "\t")

# check the duplicates of symbols column
x <- duplicated(raw_data$Gene_symbol)
sum(x)

# aggregation
raw_data <- aggregate(raw_data, list(raw_data$Gene_symbol),FUN=mean)
raw_data <- raw_data[, -2]
colnames(raw_data)[1] <- "symbol"

# check the duplicates of symbols column
x <- duplicated(raw_data$symbol)
sum(x)

# Set the row names of the dataframe to the values in the Gene_column column
rownames(raw_data) <- raw_data$symbol

# remove the symbol column
raw_data <- raw_data[, -1]

############### impute the missing values mean ##############
# Calculate the proportion of zeros in each row
prop_zeros <- rowSums(raw_data == 0) / ncol(raw_data)

# Identify rows with less than 40% zeros
rows_to_fill <- which(prop_zeros < 0.4)
imputed_genes <- rows_to_fill

# Calculate the row means for these rows
row_means <- rowMeans(raw_data[rows_to_fill, ], na.rm = TRUE)

# Replace the zeros in these rows with the row means and remove other rows
raw_data[rows_to_fill, ] <- sweep(raw_data[rows_to_fill, ], 2, row_means, FUN = function(x, y) ifelse(x == 0, y, x))
raw_data <- raw_data[imputed_genes, ]

# Round the imputed values
raw_data <- round(raw_data)

################################################
### get metadata
gse <- getGEO(GEO = "GSE205450", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))

# get needed columns and modify it
colData <- metadata %>%
  dplyr::select(1, 10, 11) %>%
  dplyr::rename(samples_name = title) %>%
  dplyr::rename(sample_condition = characteristics_ch1) %>%
  dplyr::rename(region = characteristics_ch1.1) %>% 
  dplyr::mutate(sample_condition = gsub("disease: ", "",sample_condition)) %>% 
  dplyr::mutate(region = gsub("region: ", "", region)) 

# Remove row names from a colData frame
rownames(colData) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData) <- colData$samples_name

# remove the first column from the meatadata_modified
colData <- colData[, -1]

# Add "X" to the beginning of each column name
rownames(colData) <- paste0("X", rownames(colData))

all(colnames(raw_data) %in% rownames(colData))
all(rownames(colData) %in% colnames(raw_data))

# check if that they are in the same order
all(rownames(colData) == colnames(raw_data))

#######################################
### we divide the samples to CAU and PUT
## let's begin with CAU

head(colData)

meta_CAU_data <- colData %>%
  dplyr::filter(region == "CAU")

# download the modified metadata to use it in the analysis
# write.csv(meta_CAU_data, "meta_CAU_data.csv", row.names = TRUE)

# extract the columns that matches with the rows of modified metadata
raw_CAU_data <- raw_data[, names(raw_data) %in% rownames(meta_CAU_data)]

########## preparing for the DESeq2 analysis for stage 6 data

#check the row names in colData match with column names with data_5_ordered
all(colnames(raw_CAU_data) %in% rownames(meta_CAU_data))

#check the order
all(colnames(raw_CAU_data) == rownames(meta_CAU_data))

# Convert 'sample_condition' column to factor
meta_CAU_data$sample_condition <- factor(meta_CAU_data$sample_condition)






########## Data exploration ########
exp.data <- apply(raw_CAU_data, 2, as.numeric)

# Plot histogram of processed RNA-Seq data
hist(log2(exp.data + 1), main = "processed RNA-Seq Histogram")

# Plot boxplot of processed RNA-Seq data
boxplot(log2(exp.data + 1), main = "RNA-Seq Box Plot", col = seq(1:15))

# Mean vs Variance
mean_counts <- apply(raw_CAU_data[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(raw_CAU_data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

#let’s plot a histogram of the counts for a single sample, ‘X2589_CTRL_CAU’:
ggplot(raw_CAU_data) +
geom_histogram(aes(x = X2589_CTRL_CAU), stat = "bin", bins = 200) +
xlab("Raw expression counts") +
ylab("Number of genes")
 
#let’s plot a histogram of the counts for a single sample, ‘X2653_PD_CAU’:
ggplot(raw_CAU_data) +
geom_histogram(aes(x = X2653_PD_CAU), stat = "bin", bins = 200) +
xlab("Raw expression counts") +
ylab("Number of genes")

########### Begin with the DESeq2 analysis 
# Construct a DESeqDataSet object 
ddsCAU <- DESeqDataSetFromMatrix(countData = raw_CAU_data,
                                 colData = meta_CAU_data,
                                 design = ~ sample_condition)


# set the factor level
ddsCAU$sample_condition <-  relevel(ddsCAU$sample_condition, ref = 'Control')

# Run DESeq
ddsCAU <- DESeq(ddsCAU)

resCAU <- results(ddsCAU, alpha = 0.05)

#check results
summary(resCAU)


# Remove rows with null values
resCAU <- resCAU[complete.cases(resCAU), ]

# Convert results to a data frame
resCAU.df <- as.data.frame(resCAU)

# Extract significant genes based on padj & LFC
sig_genes_CAU <- resCAU[resCAU$padj < 0.05 & abs(resCAU$log2FoldChange)>(1), ]
summary(sig_genes_CAU)
sig_genes_CAU <- as.data.frame(sig_genes_CAU)




########## Adjusting DEGs For Plotting #############

# Perform variance stabilizing transformation
vsd <- vst(ddsCAU, blind = FALSE)

# Extract the transformed data for significant genes
vst_data <- assay(vsd)

# Add a column to identify significant genes based on new criteria
resCAU.df$significant <- ifelse(resCAU.df$padj < 0.05 & abs(resCAU.df$log2FoldChange)>1 , 
                                ifelse(resCAU.df$log2FoldChange > 1, "Upregulated", 
                                       ifelse(resCAU.df$log2FoldChange < -1, "Downregulated", "Not Significant")), 
                                "Not Significant")

# Define top genes based on the highest -log10(padj) values
top_genes <- rownames(head(resCAU.df[order(resCAU.df$padj), ], 10))

# Create a new column for labeling top genes
resCAU.df$label <- ifelse(rownames(resCAU.df) %in% top_genes, rownames(resCAU.df), NA)

# Create the volcano plot
p <- ggplot(resCAU.df, aes(x = log2FoldChange, y = -log10(padj))) +
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

# Perform PCA manually
pca_res <- prcomp(t(vst_data))

# Create a data frame with the PCA results
pca_data <- as.data.frame(pca_res$x)
pca_data$sample_condition <- colData(ddsCAU)$sample_condition

# Plot the first two principal components (PC1 vs. PC2)
pca_plot1 <- ggplot(pca_data, aes(PC1, PC2, color = sample_condition)) +
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
pca_plot2 <- ggplot(pca_data, aes(PC2, PC3, color = sample_condition)) +
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
                       sample_condition = colData(ddsCAU)$sample_condition)

# Plot the 3D PCA
fig <- plot_ly(pca_data, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~sample_condition, 
               colors = c("blue", "red"),
               text = ~paste("Sample: ", rownames(pca_data)),
               marker = list(size = 5)) %>%
  layout(title = '3D PCA of RNA-seq Data',
         scene = list(xaxis = list(title = paste0("PC1: ", round(100 * summary(pca_res)$importance[2,1], 1), "% variance")),
                      yaxis = list(title = paste0("PC2: ", round(100 * summary(pca_res)$importance[2,2], 1), "% variance")),
                      zaxis = list(title = paste0("PC3: ", round(100 * summary(pca_res)$importance[2,3], 1), "% variance"))))

# Display the plot
fig


######################## Heat Map ##########################

# Sort significant genes by adjusted p-value (padj) in ascending order
sig_genes_CAU <- sig_genes_CAU[order(sig_genes_CAU$padj), ]

# Select the top 50 significant genes based on the smallest padj
top_50_sig_genes <- head(rownames(sig_genes_CAU), 50)

# Subset the original counts data to only include the top 50 significant genes
top_50_sig_gene_counts <- counts(ddsCAU, normalized=TRUE)[top_50_sig_genes, ]

top_50_sig_gene_counts <- as.data.frame(top_50_sig_gene_counts)

# Perform a log transformation for better visualization
log_top_50_sig_gene_counts <- log2(top_50_sig_gene_counts + 1)

# Check the structure of colData
print(colData(ddsCAU))

# Ensure the correct column exists and extract sample conditions
if ("sample_condition" %in% colnames(colData(ddsCAU))) {
  sample_conditions <- as.data.frame(colData(ddsCAU)$sample_condition)
  colnames(sample_conditions) <- "Condition"
  rownames(sample_conditions) <- colnames(log_top_50_sig_gene_counts)
} else {
  stop("Sample condition column not found in colData(ddsCAU). Please check the column name.")
}

# Add region "CAU" to the sample conditions
sample_conditions$Region <- "CAU"

# Print the sample_conditions to check correctness
print(sample_conditions)

# Dynamically define colors for the conditions
condition_levels <- levels(factor(sample_conditions$Condition))
annotation_colors <- list(
  Condition = setNames(c("blue", "red")[1:length(condition_levels)], condition_levels),
  Region = c(CAU = "orange") # You can choose any color you prefer
)

# Generate the heatmap with sample condition annotations
pheatmap(log_top_50_sig_gene_counts, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE, 
         show_rownames=TRUE, 
         show_colnames=TRUE, 
         scale="row",
         annotation_col=sample_conditions,
         annotation_colors=annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# Save Results
# write.csv(resCAU.degs, "all_CAU_genes.csv", row.names = TRUE)
# write.csv(sig_CAU, "CAU_significant_genes.csv", row.names = TRUE)


##########################################################################
--------------------------------------------------------------------------
############################## PUT #######################################

#######################################
### we divide the samples to CAU and PUT
## NOW with PUT

meta_PUT_data <- colData %>%
  dplyr::filter(region == "PUT")

# download the modified metadata to use it in the analysis
# write.csv(meta_PUT_data, "meta_PUT_data.csv", row.names = TRUE)

# extract the columns that matches with the rows of modified metadata
raw_PUT_data <- raw_data[, names(raw_data) %in% rownames(meta_PUT_data)]

#check the row names in colData match with column names with data_5_ordered
all(colnames(raw_PUT_data) %in% rownames(meta_PUT_data))

#check the order
all(colnames(raw_PUT_data) == rownames(meta_PUT_data))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_PUT_data, is.numeric)

# Convert numeric columns to integer format
raw_PUT_data[, numeric_columns] <- lapply(raw_PUT_data[, numeric_columns], as.integer)

# Convert 'sample_condition' column to factor
meta_PUT_data$sample_condition <- factor(meta_PUT_data$sample_condition)

######### Check the distribution of RNA-seq counts 
#let’s plot a histogram of the counts for a single sample, ‘X2590_CTRL_PUT’:
ggplot(raw_PUT_data) +
  geom_histogram(aes(x = X2590_CTRL_PUT), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

########## Data exploration ########
exp.data <- apply(raw_PUT_data, 2, as.numeric)

# Plot histogram of processed RNA-Seq data
hist(log2(exp.data + 1), main = "processed RNA-Seq Histogram")

# Plot boxplot of processed RNA-Seq data
boxplot(log2(exp.data + 1), main = "RNA-Seq Box Plot", col = seq(1:15))

# Mean vs Variance
mean_counts <- apply(raw_PUT_data[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(raw_PUT_data[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

########### Begin with the DESeq2 analysis 
# Construct a DESeqDataSet object 
ddsPUT <- DESeqDataSetFromMatrix(countData = raw_PUT_data,
                                 colData = meta_PUT_data,
                                 design = ~ sample_condition)


# set the factor level
ddsPUT$sample_condition <-  relevel(ddsPUT$sample_condition, ref = 'Control')

# Run DESeq
ddsPUT <- DESeq(ddsPUT)
resPUT <- results(ddsPUT, alpha = 0.05)

#check results
summary(resPUT)

# Remove rows with null values
resPUT <- resPUT[complete.cases(resPUT), ]

# Convert results to a data frame
resPUT.df <- as.data.frame(resPUT)

# Extract significant genes based on padj & LFC
sig_genes_PUT <- resPUT[resPUT$padj < 0.05 & abs(resPUT$log2FoldChange)>(1), ]
summary(sig_genes_PUT)
sig_genes_PUT <- as.data.frame(sig_genes_PUT)


########## Adjusting DEGs For Plotting #############

# Perform variance stabilizing transformation
vsd <- vst(ddsPUT, blind = FALSE)

# Extract the transformed data for significant genes
vst_data <- assay(vsd)

# Add a column to identify significant genes based on new criteria
resPUT.df$significant <- ifelse(resPUT.df$padj < 0.05 & abs(resPUT.df$log2FoldChange)>1 , 
                                ifelse(resPUT.df$log2FoldChange > 1, "Upregulated", 
                                       ifelse(resPUT.df$log2FoldChange < -1, "Downregulated", "Not Significant")), 
                                "Not Significant")

# Define top genes based on the highest -log10(padj) values
top_genes <- rownames(head(resPUT.df[order(resPUT.df$padj), ], 10))

# Create a new column for labeling top genes
resPUT.df$label <- ifelse(rownames(resPUT.df) %in% top_genes, rownames(resPUT.df), NA)

# Create the volcano plot
p <- ggplot(resPUT.df, aes(x = log2FoldChange, y = -log10(padj))) +
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

# Perform PCA manually
pca_res <- prcomp(t(vst_data))

# Create a data frame with the PCA results
pca_data <- as.data.frame(pca_res$x)
pca_data$sample_condition <- colData(ddsPUT)$sample_condition

# Plot the first two principal components (PC1 vs. PC2)
pca_plot1 <- ggplot(pca_data, aes(PC1, PC2, color = sample_condition)) +
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
pca_plot2 <- ggplot(pca_data, aes(PC2, PC3, color = sample_condition)) +
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
                       sample_condition = colData(ddsPUT)$sample_condition)

# Plot the 3D PCA
fig <- plot_ly(pca_data, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~sample_condition, 
               colors = c("blue", "red"),
               text = ~paste("Sample: ", rownames(pca_data)),
               marker = list(size = 5)) %>%
  layout(title = '3D PCA of RNA-seq Data',
         scene = list(xaxis = list(title = paste0("PC1: ", round(100 * summary(pca_res)$importance[2,1], 1), "% variance")),
                      yaxis = list(title = paste0("PC2: ", round(100 * summary(pca_res)$importance[2,2], 1), "% variance")),
                      zaxis = list(title = paste0("PC3: ", round(100 * summary(pca_res)$importance[2,3], 1), "% variance"))))

# Display the plot
fig


######################## Heat Map ##########################

# Sort significant genes by adjusted p-value (padj)
sig_genes_PUT <- sig_genes_PUT[order(sig_genes_PUT$padj), ]

# Select the top 50 significant genes based on padj
top_50_sig_genes <- head(rownames(sig_genes_PUT), 50)

# Subset the original counts data to only include the top 50 significant genes
top_50_sig_gene_counts <- counts(ddsPUT, normalized=TRUE)[top_50_sig_genes, ]

# Perform a log transformation for better visualization
log_top_50_sig_gene_counts <- log2(top_50_sig_gene_counts + 1)

# Check the structure of colData
print(colData(ddsPUT))

# Ensure the correct column exists and extract sample conditions
if ("sample_condition" %in% colnames(colData(ddsPUT))) {
  sample_conditions <- as.data.frame(colData(ddsPUT)$sample_condition)
  colnames(sample_conditions) <- "Condition"
  rownames(sample_conditions) <- colnames(log_top_50_sig_gene_counts)
} else {
  stop("Sample condition column not found in colData(ddsPUT). Please check the column name.")
}

# Add region "PUT" to the sample conditions
sample_conditions$Region <- "PUT"

# Print the sample_conditions to check correctness
print(sample_conditions)

# Dynamically define colors for the conditions
condition_levels <- levels(factor(sample_conditions$Condition))
annotation_colors <- list(
  Condition = setNames(c("blue", "red")[1:length(condition_levels)], condition_levels),
  Region = c(PUT = "lightblue") # You can choose any color you prefer
)

# Generate the heatmap with sample condition annotations
pheatmap(log_top_50_sig_gene_counts, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE, 
         show_rownames=TRUE, 
         show_colnames=TRUE, 
         scale="row",
         annotation_col=sample_conditions,
         annotation_colors=annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# save results
write.csv(resPUT.degs, "all_PUT_genes.csv", row.names = TRUE)
write.csv(sig_PUT, "PUT_significant_genes.csv", row.names = TRUE)
