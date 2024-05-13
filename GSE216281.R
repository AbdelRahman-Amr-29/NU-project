############### Information about the data
## Summary:	To investigate transcriptomic changes associated with Parkinson's disease development 
#we compared frontal cortex gene expression across four Braak Lewy body stage groups.
#Additionally we investigated sex-specific gene expression differences in neuropathologically
#healthy donors and Parkinson's disease patients at Braak Lewy body stage 5.
#We investigated sex-specific gene expression differences in neuropathologically
#healthy donors and Parkinson's disease patients at Braak Lewy body stage 5.
##DOI: https://doi.org/10.1007/s00401-023-02597-7

###################################################################
# Setting working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281")

# installing packages
#BiocManager::install("biomaRt")
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install('airway')
#install.packages("pheatmap")

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

#################################################
# Read txt file & convert it to table
raw_data <- read.table("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/GSE216281_raw_data.txt", header = TRUE, sep = " ")

# Bind the row names with the "gene" column
raw_data <- cbind(ensembl_gene_id = rownames(raw_data), raw_data)

# Remove row names from a raw_data frame
rownames(raw_data) <- NULL

# extract first column
ids <- raw_data[ , 1, drop = FALSE]
head(ids)

# download the first column from the data file contains ids
#write.csv(ids, "ensembl_gene_ids.csv", row.names = FALSE)

# read the ensembl_gene_ids file to use it in the next step
ensembl_ids <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/ensembl_gene_ids.csv", header = TRUE)

########################################################
# convert gene ids to gene names
keytypes(org.Hs.eg.db)
columns(org.Hs.eg.db)

mapper <- mapIds(org.Hs.eg.db,
                 keys = ensembl_ids$ensembl_gene_id,
                 keytype = "ENSEMBL",
                 column = "SYMBOL")


mapper.df <- as.data.frame(mapper)
mapper.df <- cbind(rownames(mapper.df), mapper.df)
names(mapper.df) <- c("ensembl_gene_id","symbol")
raw_data <- merge(mapper.df, raw_data, by = "ensembl_gene_id", all.x = TRUE)
raw_data <- raw_data[, -1]
raw_data <- raw_data[ ! is.na(raw_data$symbol),]

########################################################
# check the duplicates of symbols column
x <- duplicated(raw_data$symbol)  
sum(x)

####remove  duplication by aggregation
raw_data.agg <- aggregate(raw_data, list(raw_data$symbol),FUN=mean)
raw_data.agg <- raw_data.agg[, -2]
colnames(raw_data.agg)[1] <- "symbol"

# Set the 'symbol' column as row names
row.names(raw_data.agg) <- raw_data.agg$symbol

# remove symbol column
raw_data.agg <- raw_data.agg[ , -1]

#####################################################
# Remove rows containing zeros
raw_data.agg <- raw_data.agg[apply(raw_data.agg, 1, function(row) !any(row == 0)), ]

# Check if rows containing zeros are removed
zeros_removed <- any(apply(raw_data.agg, 1, function(row) any(row == 0)))

# Print the result
if (!zeros_removed) {
  print("Rows containing zeros have been removed.")
} else {
  print("Rows containing zeros could not be removed completely.")
}

########################################################
### get metadata
gse <- getGEO(GEO = "GSE216281", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(21, 1, 12) %>%
  dplyr::rename(samples_name = description.1) %>%
  dplyr::rename(sample_condition = characteristics_ch1.2) %>%
  dplyr::mutate(sample_condition = gsub("braak lewy_body_stage: ", "", sample_condition)) 

# Remove row names from a metadata_modified frame
rownames(metadata_modified) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(metadata_modified) <- metadata_modified$samples_name

# remove the first column from the meatadata_modified
metadata_modified <- metadata_modified[, -1]

all(colnames(raw_data.agg) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(raw_data.agg))

#######################################
# know the right order indexing
data_idx <- match(rownames(metadata_modified), colnames(raw_data.agg))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
data_ordered  <- raw_data.agg[ , data_idx]

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(data_ordered))

###################################
# adjusting metadata (colData) for DESeq2 analysis
colData <- metadata_modified %>% 
  mutate(sample_condition = gsub(0, "control", sample_condition)) %>%
  mutate(sample_condition = gsub(1, "PD", sample_condition)) %>% 
  mutate(sample_condition = gsub(2, "PD", sample_condition)) %>% 
  mutate(sample_condition = gsub(3, "PD", sample_condition)) %>% 
  mutate(sample_condition = gsub(4, "PD", sample_condition)) %>% 
  mutate(sample_condition = gsub(5, "PD", sample_condition)) %>% 
  mutate(sample_condition = gsub(6, "PD", sample_condition))

# download the modified metadata to use it in the analysis
#write.csv(colData, "colData_all.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis
# read in the metadata we have adjusted lately
colData <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData_all.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData) <- colData$X

# remove the first column from the meatadata_modified
colData <- colData[, -1]

#check the row names in colData match with column names with data_5_ordered
all(colnames(data_ordered) %in% rownames(colData))

#check the order
all(colnames(data_ordered) == rownames(colData))

# Convert numeric columns to integer
numeric_columns <- sapply(data_ordered, is.numeric)

# Convert numeric columns to integer format
data_ordered[, numeric_columns] <- lapply(data_ordered[, numeric_columns], as.integer)

# Convert 'sample_condition' column to factor
colData$sample_condition <- factor(colData$sample_condition)
class(colData$sample_condition)

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘sample_R1’
ggplot(data_ordered) +
  geom_histogram(aes(x = sample_R1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘sample_R6’
ggplot(data_ordered) +
  geom_histogram(aes(x = sample_R6), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Mean vs Variance
mean_counts <- apply(data_ordered[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data_ordered[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

########### Begin with the DESeq2 analysis 
# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = data_ordered,
                              colData = colData,
                              design = ~ sample_condition)

view(counts(dds))

# remove rows that have value less than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# set the factor level
dds$sample_condition <-  relevel(dds$sample_condition, ref = "control")

# Run DESeq
dds <- DESeq(dds)

res <- results(dds)

res

# Explore results
summary(res)

res0.05 <- results(dds, alpha = 0.05)
summary(res0.05)

##########################################################
### sample condition 1-4 comparing with control 0
meta_1_data <- metadata_modified %>%
  filter(sample_condition == 0 | sample_condition == 1 | sample_condition == 2
         | sample_condition == 3 | sample_condition == 4)

# change the 0 to control and 1,2,3 and 4 to PD
colData_1 <- meta_1_data %>% 
  dplyr::mutate(sample_condition = gsub(0, "control", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(1, "PD", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(2, "PD", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(3, "PD", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(4, "PD", sample_condition))

# download the modified metadata to use it in the analysis
#write.csv(colData_1, "colData1.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis for stage 5
# read in the metadata we have adjusted lately
colData_1 <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData1.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData_1) <- colData_1$X

# remove the first column from the meatadata_modified
colData_1 <- colData_1[, -1]

# keep only the columns that matches with the rows of modified metadata
raw_data1 <- raw_data.agg[, names(raw_data.agg) %in% rownames(colData_1)]

#check the row names in colData match with column names with raw_data1
all(colnames(raw_data1) %in% rownames(colData_1))

########################################
# know the right order indexing
data_idx <- match(rownames(colData_1), colnames(raw_data1))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
raw_data1_ordered  <- raw_data1[ , data_idx]

# check if that they are in the same order
all(rownames(colData_1) == colnames(raw_data1_ordered))

#check the order
all(colnames(raw_data1_ordered) == rownames(colData_1))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_data1_ordered, is.numeric)

# Convert numeric columns to integer format
raw_data1_ordered[, numeric_columns] <- lapply(raw_data1_ordered[, numeric_columns], as.integer)

# Convert 'braak_lewy_body_stage' column to factor
colData_1$sample_condition <- factor(colData_1$sample_condition)
class(colData_1$sample_condition)

######################################################################
######### Check the distribution of RNA-seq counts 
#let’s plot a histogram of the counts for a single sample, ‘sample_R6’:
ggplot(raw_data1_ordered) +
  geom_histogram(aes(x = sample_R6), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

#let’s plot a histogram of the counts for a single sample, ‘sample_R10’:
ggplot(raw_data1_ordered) +
  geom_histogram(aes(x = sample_R10), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Mean vs Variance
mean_counts <- apply(raw_data1_ordered[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(raw_data1_ordered[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

########### Begin the DESeq2 analysis 
# Construct a DESeqDataSet object 
dds1 <- DESeqDataSetFromMatrix(countData = raw_data1_ordered,
                               colData = colData_1,
                               design = ~ sample_condition)

view(counts(dds1))


# remove rows that have value less than 10 
keep <- rowSums(counts(dds1)) >= 10
dds1 <- dds1[keep, ]

# set the factor level
dds1$sample_condition <-  relevel(dds1$sample_condition, ref = "control")

# Run DESeq
dds1 <- DESeq(dds1)
res1 <- results(dds1)

res1

#check results
summary(res1)

res1_0.05 <- results(dds1, alpha = 0.05)
summary(res1_0.05)

#######################################
### sample condition 5 comparing with control 0
meta_5_data <- metadata_modified %>%
  filter(sample_condition == 0 | sample_condition == 5)

# change the 0 to control and 5 to PD stage 5 in the column
colData_5 <- meta_5_data %>% 
  dplyr::mutate(sample_condition = gsub(0, "control", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(5, "PD", sample_condition))

# download the modified metadata to use it in the analysis
#write.csv(colData_5, "colData5.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis for stage 5
# read in the metadata we have adjusted lately
colData_5 <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData5.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData_5) <- colData_5$X

# remove the first column from the meatadata_modified
colData_5 <- colData_5[, -1]

# keep only the columns that matches with the rows of colData_5
data_5 <- raw_data.agg[, names(raw_data.agg) %in% rownames(colData_5)]

# know the right order indexing
data_idx <- match(rownames(colData_5), colnames(data_5))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
data_5_ordered  <- data_5[ , data_idx]

#check the row names in colData match with column names with data_5_ordered
all(colnames(data_5_ordered) %in% rownames(colData_5))

#check the order
all(colnames(data_5_ordered) == rownames(colData_5))

# Convert numeric columns to integer
numeric_columns <- sapply(data_5_ordered, is.numeric)

# Convert numeric columns to integer format
data_5_ordered[, numeric_columns] <- lapply(data_5_ordered[, numeric_columns], as.integer)

# Convert 'braak_lewy_body_stage' column to factor
colData_5$sample_condition <- factor(colData_5$sample_condition)
class(colData_5$sample_condition)

######### Check the distribution of RNA-seq counts 
#let’s plot a histogram of the counts for a single sample, ‘sample_R1’:
ggplot(data_5_ordered) +
  geom_histogram(aes(x = sample_R1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Mean vs Variance
mean_counts <- apply(data_5_ordered[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data_5_ordered[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

########### Begin the DESeq2 analysis 
# Construct a DESeqDataSet object 
dds5 <- DESeqDataSetFromMatrix(countData = data_5_ordered,
                              colData = colData_5,
                              design = ~ sample_condition)

view(counts(dds5))


# remove rows that have value less than 10 
keep <- rowSums(counts(dds5)) >= 10
dds5 <- dds5[keep, ]

# set the factor level
dds5$sample_condition <-  relevel(dds5$sample_condition, ref = 'control')

# Run DESeq
dds5 <- DESeq(dds5)
res5 <- results(dds5)

res5

#check results
summary(res5)

res5_0.05 <- results(dds5, alpha = 0.05)
summary(res5_0.05)

#####################################################
### sample condition 6 comparing with control 0
meta_6_data <- metadata_modified %>%
  filter(sample_condition == 0 | sample_condition == 6)

# change the 0 to control and 5 to PD stage 5 in the column
colData_6 <- meta_6_data %>% 
  dplyr::mutate(sample_condition = gsub(0, "control", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(6, "PD", sample_condition))

# download the modified metadata to use it in the analysis
#write.csv(colData_6, "colData6.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis for stage 6 data
# read in the metadata we have adjusted lately
colData_6 <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData6.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData_6) <- colData_6$X

# remove the first column from the meatadata_modified
colData_6 <- colData_6[, -1]

# keep only the columns that matches with the rows of colData_5
data_6 <- raw_data.agg[, names(raw_data.agg) %in% rownames(colData_6)]

# know the right order indexing
data_idx <- match(rownames(colData_6), colnames(data_6))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
data_6_ordered  <- data_6[ , data_idx]

#check the row names in colData match with column names with data_5_ordered
all(colnames(data_6_ordered) %in% rownames(colData_6))

#check the order
all(colnames(data_6_ordered) == rownames(colData_6))

# Convert numeric columns to integer
numeric_columns <- sapply(data_6_ordered, is.numeric)

# Convert numeric columns to integer format
data_6_ordered[, numeric_columns] <- lapply(data_6_ordered[, numeric_columns], as.integer)

# Convert 'braak_lewy_body_stage' column to factor
colData_6$sample_condition <- factor(colData_6$sample_condition)
class(colData_6$sample_condition)

######### Check the distribution of RNA-seq counts 
#let’s plot a histogram of the counts for a single sample, ‘sample_R6’:
ggplot(data_6_ordered) +
  geom_histogram(aes(x = sample_R6), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Mean vs Variance
mean_counts <- apply(data_6_ordered[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(data_6_ordered[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

########### Begin the DESeq2 analysis 
# Construct a DESeqDataSet object 
dds6 <- DESeqDataSetFromMatrix(countData = data_6_ordered,
                               colData = colData_6,
                               design = ~ sample_condition)

view(counts(dds6))


# remove rows that have value less than 10 
keep <- rowSums(counts(dds6)) >= 10
dds6 <- dds6[keep, ]

# set the factor level
dds6$sample_condition <-  relevel(dds6$sample_condition, ref = 'control')

# Run DESeq
dds6 <- DESeq(dds6)
res6 <- results(dds6)

res6

#check results
summary(res6)

res6_0.05 <- results(dds6, alpha = 0.05)
summary(res6_0.05)
