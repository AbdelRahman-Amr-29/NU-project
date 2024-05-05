## setting the working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE197345")

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
raw_counts <- read_excel("raw_counts.xlsx")
View(raw_counts)

# make raw_counts data frame
raw_counts <-  data.frame(raw_counts)
class(raw_counts)

##### check the duplicates of gene column
x <- duplicated(raw_counts$gene)
sum(x)

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
gse <- getGEO(GEO = "GSE197345", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(1, 10, 12, 11) %>%
  dplyr::rename(diagnosis = characteristics_ch1) %>%
  dplyr::rename(age = characteristics_ch1.2) %>%
  dplyr::rename(gender = characteristics_ch1.1) %>% 
  dplyr::mutate(diagnosis = gsub("diagnosis: Control", "control", diagnosis)) %>%
  dplyr::mutate(diagnosis = gsub("diagnosis: Essential Tremor", "ET", diagnosis)) %>%
  dplyr::mutate(age = gsub("age: ", "", age)) %>%
  dplyr::mutate(gender = gsub("gender: ", "", gender))
  
head(metadata_modified)

# Remove row names from a metadata_modified frame
rownames(metadata_modified) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(metadata_modified) <- metadata_modified$title

# remove the first column from the meatadata_modified
metadata_modified <- metadata_modified[, -1]

all(colnames(raw_counts) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(raw_counts))

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(raw_counts))

# download the modified metadata to use it in the analysis
write.csv(metadata_modified, "colData.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis
# read in the metadata we have adjusted lately
colData <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE197345/colData.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData) <- colData$X

# remove the first column from the meatadata_modified
colData <- colData[, -1]

#check the row names in colData match with column names with raw_counts
all(colnames(raw_counts) %in% rownames(colData))

#check the order
all(colnames(raw_counts) == rownames(colData))

# Convert 'diagnosis' column to factor
colData$diagnosis <- factor(colData$diagnosis)

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘A_1’
ggplot(raw_counts) +
  geom_histogram(aes(x = A_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘C_5’
ggplot(raw_counts) +
  geom_histogram(aes(x = C_5), stat = "bin", bins = 200) +
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
