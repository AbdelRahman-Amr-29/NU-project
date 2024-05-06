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

## read the samples
raw_counts <- read_excel("GSE134878_raw_read_counts.xls")
View(raw_counts)

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

#####################################################
# Remove rows containing zeros
raw_counts <- raw_counts[apply(raw_counts, 1, function(row) !any(row == 0)), ]

# Check if rows containing zeros are removed
zeros_removed <- any(apply(raw_counts, 1, function(row) any(row == 0)))

# Print the result
if (!zeros_removed) {
  print("Rows containing zeros have been removed.")
} else {
  print("Rows containing zeros could not be removed completely.")
}

#####################################################
# Remove rows with NA values
raw_counts <- na.omit(raw_counts)

#####################################################
# get metadata
gse <- getGEO(GEO = "GSE134878", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
metadata <- metadata[c(-29, -47),]
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(1, 11) %>%
  dplyr::rename(diagnosis = characteristics_ch1.1) %>%
  dplyr::rename(sample_name = title) %>%
  dplyr::mutate(diagnosis = gsub("diagnosis: Control", "control", diagnosis)) %>%
  dplyr::mutate(diagnosis = gsub("diagnosis: Essential Tremor", "ET", diagnosis)) 

head(metadata_modified)

# Remove row names from a metadata_modified frame
rownames(metadata_modified) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(metadata_modified) <- metadata_modified$sample_name

# remove the sample_name column from the meatadata_modified
metadata_modified <- subset(metadata_modified, select = -sample_name)

all(colnames(raw_counts) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(raw_counts))

# know the right order indexing
data_idx <- match(rownames(metadata_modified), colnames(raw_counts))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
raw_counts  <- raw_counts[ , data_idx]

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(raw_counts))

# download the modified metadata to use it in the analysis
#write.csv(metadata_modified, "colData.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis
# remove rows form the raw_counts that don't have gene name
raw_counts <- raw_counts[c(-1:-21), ]

# read in the metadata we have adjusted lately
colData <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE134878/colData.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData) <- colData$X

# remove the first column from the meatadata_modified
colData <- subset(colData, select = -X)

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

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘A1’
ggplot(raw_counts) +
  geom_histogram(aes(x = A1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘B7’
ggplot(raw_counts) +
  geom_histogram(aes(x = B7), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

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

view(counts(dds))

# remove rows that have value less than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

view(counts(dds))

# set the factor level
dds$diagnosis <-  relevel(dds$diagnosis, ref = "control")

##### Run DESeq
dds <- DESeq(dds)

res <- results(dds)
res

# Explore results
summary(res)

res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

res0.25 <- results(dds, alpha = 0.25)
summary(res0.25)
