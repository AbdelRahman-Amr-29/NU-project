################### Information about the data
## GSE205450
## Summary:	We analyzed transcriptome differences in postmortem caudate and putamen 
#from controls (n=40) and PD (n=35) as confirmed by autopsy.
#Further, we analyzed region specific differences in caudate and putamen
#associated with clinical variables.
# NOTE: For differential gene expression analysis between controls and PD,
#all 35 PD and 40 control specimens (caudate and putamen) obtained from NIH NeuroBiobank
#were included and processed for bulk RNA-seq as INDEPENDENT biological replicates.
## DOI: https://doi.org/10.1038/s41467-023-39652-6

#############################################################################
## setting working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE205450")

# calling libraries
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

# Read txt file & convert it to table
raw_data <- read.table("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE205450/GSE205450_counts.table.txt", header = TRUE, sep = "\t")

# check the duplicates of symbols column
x <- duplicated(raw_data$Gene_symbol)
sum(x)

# aggregation
raw_data <- aggregate(raw_data, list(raw_data$Gene_symbol),FUN=mean)
raw_data <- raw_data[, -2]
colnames(raw_data)[1] <- "symbol"

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
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(1, 10, 11) %>%
  dplyr::rename(samples_name = title) %>%
  dplyr::rename(sample_condition = characteristics_ch1) %>%
  dplyr::rename(region = characteristics_ch1.1) %>% 
  dplyr::mutate(sample_condition = gsub("disease: ", "",sample_condition)) %>% 
  dplyr::mutate(region = gsub("region: ", "", region)) 

# Remove row names from a metadata_modified frame
rownames(metadata_modified) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(metadata_modified) <- metadata_modified$samples_name

# remove the first column from the meatadata_modified
metadata_modified <- metadata_modified[, -1]

# Add "X" to the beginning of each column name
rownames(metadata_modified) <- paste0("X", rownames(metadata_modified))

all(colnames(raw_data) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(raw_data))

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(raw_data))

#######################################
### we divide the samples to CAU and PUT
## let's begin with CAU

head(metadata_modified)

meta_CAU_data <- metadata_modified %>%
  dplyr::filter(region == "CAU")

# download the modified metadata to use it in the analysis
# write.csv(meta_CAU_data, "colDataCAU.csv", row.names = TRUE)

# extract the columns that matches with the rows of modified metadata
raw_CAU_data <- raw_data[, names(raw_data) %in% rownames(meta_CAU_data)]

########## preparing for the DESeq2 analysis for stage 6 data
# read in the metadata we have adjusted lately
colDataCAU <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE205450/colDataCAU.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colDataCAU) <- colDataCAU$X

# remove the first column from the meatadata_modified
colDataCAU <- colDataCAU[, -1]

#check the row names in colData match with column names with data_5_ordered
all(colnames(raw_CAU_data) %in% rownames(colDataCAU))

#check the order
all(colnames(raw_CAU_data) == rownames(colDataCAU))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_CAU_data, is.numeric)

# Convert 'sample_condition' column to factor
colDataCAU$sample_condition <- factor(colDataCAU$sample_condition)

######### Check the distribution of RNA-seq counts 
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

########### Begin with the DESeq2 analysis 
# Construct a DESeqDataSet object 
ddsCAU <- DESeqDataSetFromMatrix(countData = raw_CAU_data,
                               colData = colDataCAU,
                               design = ~ sample_condition)

view(counts(ddsCAU))


# remove rows that have value less than 10 
keep <- rowSums(counts(ddsCAU)) >= 10
ddsCAU <- ddsCAU[keep, ]

# set the factor level
ddsCAU$sample_condition <-  relevel(ddsCAU$sample_condition, ref = 'Control')

# Run DESeq
ddsCAU <- DESeq(ddsCAU)
resCAU <- results(ddsCAU)

resCAU

#check results
summary(resCAU)

resCAU_0.05 <- results(ddsCAU, alpha = 0.05)
summary(resCAU_0.05)

#######################################
### we divide the samples to CAU and PUT
## NOW with PUT

head(metadata_modified)

meta_PUT_data <- metadata_modified %>%
  dplyr::filter(region == "PUT")

# download the modified metadata to use it in the analysis
# write.csv(meta_PUT_data, "colDataPUT.csv", row.names = TRUE)

# extract the columns that matches with the rows of modified metadata
raw_PUT_data <- raw_data[, names(raw_data) %in% rownames(meta_PUT_data)]

########## preparing for the DESeq2 analysis for stage 6 data
# read in the metadata we have adjusted lately
colDataPUT <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE205450/colDataPUT.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colDataPUT) <- colDataPUT$X

# remove the first column from the meatadata_modified
colDataPUT <- colDataPUT[, -1]

#check the row names in colData match with column names with data_5_ordered
all(colnames(raw_PUT_data) %in% rownames(colDataPUT))

#check the order
all(colnames(raw_PUT_data) == rownames(colDataPUT))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_PUT_data, is.numeric)

# Convert numeric columns to integer format
raw_PUT_data[, numeric_columns] <- lapply(raw_PUT_data[, numeric_columns], as.integer)

# Convert 'sample_condition' column to factor
colDataPUT$sample_condition <- factor(colDataPUT$sample_condition)

######### Check the distribution of RNA-seq counts 
#let’s plot a histogram of the counts for a single sample, ‘X2590_CTRL_PUT’:
ggplot(raw_PUT_data) +
  geom_histogram(aes(x = X2590_CTRL_PUT), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#let’s plot a histogram of the counts for a single sample, ‘X2654_PD_PUT’:
ggplot(raw_PUT_data) +
  geom_histogram(aes(x = X2654_PD_PUT), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

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
                                 colData = colDataPUT,
                                 design = ~ sample_condition)

view(counts(ddsPUT))


# remove rows that have value less than 10 
keep <- rowSums(counts(ddsPUT)) >= 10
ddsPUT <- ddsPUT[keep, ]

# set the factor level
ddsPUT$sample_condition <-  relevel(ddsPUT$sample_condition, ref = 'Control')

# Run DESeq
ddsPUT <- DESeq(ddsPUT)
resPUT <- results(ddsPUT)

resPUT

#check results
summary(resPUT)

resPUT_0.05 <- results(ddsPUT, alpha = 0.05)
summary(resPUT_0.05)
