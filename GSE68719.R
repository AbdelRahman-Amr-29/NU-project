############# Information about the Data ##########
#
# GSE68719
#
# Summary: Parkinson disease (PD) is a neurodegenerative disease characterized by the accumulation
# of alpha-synuclein (SNCA) and other proteins in aggregates termed “Lewy Bodies” within neurons.
# PD has both genetic and environmental risk factors, and while processes leading to aberrant protein
# aggregation are unknown, past work points to abnormal levels of SNCA and other proteins.
# Although several genome-wide studies have been performed for PD, these have focused on
# DNA sequence variants by genome-wide association studies (GWAS) and on RNA levels
# (microarray transcriptomics), while genome-wide proteomics analysis has been lacking.
# After appropriate filters, proteomics identified 3,558 unique proteins and 283 of these (7.9%)
# were significantly different between PD and controls (q-value<0.05).
# RNA-sequencing identified 17,580 protein-coding genes and 1,095 of
# these (6.2%) were significantly different (FDR p-value<0.05), but only 166 of
# the FDR significant protein-coding genes (0.94%) were present among the 3,558 proteins
# characterized. Of these 166, eight genes (4.8%) were significant in both studies,
# with the same direction of effect. Functional enrichment analysis of the proteomics
# results strongly supports mitochondrial-related pathways, while comparable analysis of
# the RNA-sequencing results implicates protein folding pathways and metallothioneins.
# Ten of the implicated genes or proteins co-localized to GWAS loci.
# Evidence implicating SNCA was stronger in proteomics than in RNA-sequencing analyses.
# Notably, differentially expressed protein-coding genes were more likely to not be characterized
# in the proteomics analysis, which lessens the ability to compare across platforms.
# Combining multiple genome-wide platforms offers novel insights into the pathological
# processes responsible for this disease by identifying pathways implicated across methodologies.
#
# DOI: 10.1186/s12920-016-0164-y
#######################################################################

# Setting The Working Directory
setwd("D:/Nile University/OneDrive - Nile University/Desktop/PD_NEW")

## calling libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
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

### Upload the Data
raw_counts <- read_excel("D:/Nile University/OneDrive - Nile University/Desktop/PD_NEW/GSE68719_mlpd_PCG_DESeq2_norm_counts.xlsx")


##### check the duplicates of gene column
x <- duplicated(raw_counts$symbol)
sum(x)

# aggregation
raw_counts <- aggregate(raw_counts, list(raw_counts$symbol),FUN=mean)
raw_counts <- raw_counts[, c(-2,-3)]
colnames(raw_counts)[1] <- "gene"

##### check the duplicates of gene column
x <- duplicated(raw_counts$symbol)
sum(x)

# make the column gene row names 
rownames(raw_counts) <- raw_counts$gene

# remove the gene column
raw_counts <- raw_counts[, -1]

# remove rows that have value less than 10 
# raw_counts <- raw_counts[rowSums(raw_counts) >=10, ]
# view(counts(dds_sn))


# Calculate the proportion of zeros in each row
prop_zeros <- rowSums(raw_counts == 0) / ncol(raw_counts)

# Identify rows with less than 40% zeros
rows_to_fill <- which(prop_zeros < 0.4)
imputed_genes <- rows_to_fill


# Calculate the row means for these rows
row_means <- rowMeans(raw_counts[rows_to_fill, ], na.rm = TRUE)

# Replace the zeros in these rows with the row means and remove other rows
raw_counts[rows_to_fill, ] <- sweep(raw_counts[rows_to_fill, ], 2, row_means, FUN = function(x, y) ifelse(x == 0, y, x))
raw_counts <- raw_counts[imputed_genes, ]

# Round the imputed values
raw_counts <- round(raw_counts)

#---------------------------------------------------------------------------#

# get metadata
gse <- getGEO(GEO = "GSE68719", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
colData <- metadata %>%
  dplyr::select(1)

# Save the File To Modify it
# write.csv(colData, "colData.csv")

# Upload The Metadata After modifying it
colData <- read_excel("colData.xlsx")

# ---> Dataframe
colData <- data.frame(colData)
class(colData)

# make the column gene row names 
rownames(colData) <- colData$title

#-------------------------------------------------------------------------#

#check the row names in colData_amg match with column names with raw_counts
all(colnames(raw_counts) %in% rownames(colData))
all(rownames(colData) %in% colnames(raw_counts))

# check if that they are in the same order
all(rownames(colData) == colnames(raw_counts))

# know the right order indexing
data_idx <- match(rownames(colData), colnames(raw_counts))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
raw_counts  <- raw_counts[ , data_idx]

# check if that they are in the same order
all(rownames(colData) == colnames(raw_counts))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_counts, is.numeric)

# Convert numeric columns to integer format
raw_counts[, numeric_columns] <- lapply(raw_counts[, numeric_columns], as.integer)

# Convert 'condition' column to factor
colData$condition <- factor(colData$condition)
class(colData$condition)

#-------------------------------------------------------------------#
######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘GSM3900032’
ggplot(raw_counts) +
  geom_histogram(aes(x = 	C_0002), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘GSM3900040’
ggplot(raw_counts) +
  geom_histogram(aes(x = P_0003), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

########## Data exploration ########
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
                                 design = ~ condition)

# set the factor level
dds$condition <-  relevel(dds$condition, ref = "control")

##### Run DESeq
dds <- DESeq(dds)

res <- results(dds)
res

# Explore results
summary(res)

res <- results(dds, alpha = 0.05)
summary(res)

# Remove rows with null values
res <- res[complete.cases(res), ]

# Convert results to a data frame
res.df <- as.data.frame(res)

# Extract significant genes based on padj & LFC
sig_genes <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
summary(sig_genes)
sig_genes <- as.data.frame(sig_genes)

# write.csv(sig_genes,"new_sig_genes.csv")


########## Adjusting DEGs For Plotting #############

# Perform variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Extract the transformed data for significant genes
vst_data <- assay(vsd)

# Add a column to identify significant genes based on new criteria
res.df$significant <- ifelse(res.df$padj < 0.05 & abs(res.df$log2FoldChange)>1 , 
                                ifelse(res.df$log2FoldChange > 1, "Upregulated", 
                                       ifelse(res.df$log2FoldChange < -1, "Downregulated", "Not Significant")), 
                                "Not Significant")

# Define top genes based on the highest -log10(padj) values
top_genes <- rownames(head(res.df[order(res.df$padj), ], 10))

# Create a new column for labeling top genes
res.df$label <- ifelse(rownames(res.df) %in% top_genes, rownames(res.df), NA)

# Create the volcano plot
p <- ggplot(res.df, aes(x = log2FoldChange, y = -log10(padj))) +
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
pca_data$condition <- colData(dds)$condition

# Plot the first two principal components (PC1 vs. PC2)
pca_plot1 <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
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
pca_plot2 <- ggplot(pca_data, aes(PC2, PC3, color = condition)) +
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
                       condition = colData(dds)$condition)

# Plot the 3D PCA
fig <- plot_ly(pca_data, 
               x = ~PC1, 
               y = ~PC2, 
               z = ~PC3, 
               color = ~condition, 
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
sig_genes <- sig_genes[order(sig_genes$padj), ]

# Select the top 50 significant genes based on padj
top_50_sig_genes <- head(rownames(sig_genes), 50)

# Subset the original counts data to only include the top 50 significant genes
top_50_sig_gene_counts <- counts(dds, normalized=TRUE)[top_50_sig_genes, ]

# Perform a log transformation for better visualization
log_top_50_sig_gene_counts <- log2(top_50_sig_gene_counts + 1)

# Check the structure of colData
print(colData(dds))

# Ensure the correct column exists and extract sample conditions
if ("condition" %in% colnames(colData(dds))) {
  conditions <- as.data.frame(colData(dds)$condition)
  colnames(conditions) <- "Condition"
  rownames(conditions) <- colnames(log_top_50_sig_gene_counts)
} else {
  stop("Sample condition column not found in colData(dds). Please check the column name.")
}


# Print the conditions to check correctness
print(conditions)

# Dynamically define colors for the conditions
condition_levels <- levels(factor(conditions$Condition))
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
         annotation_col=conditions,
         annotation_colors=annotation_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50))