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
write.csv(ids, "ensembl_gene_ids.csv", row.names = FALSE)

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
combined_data <- merge(mapper.df, raw_data, by = "ensembl_gene_id", all.x = TRUE)
combined_data <- combined_data[, -1]
combined_data <- combined_data[ ! is.na(combined_data$symbol),]

########################################################
# check the duplicates of symbols column
x <- duplicated(combined_data$symbol)  
sum(x)

### yes .. why ? transcripts?  solutions : aggregation
# do not run for now 
# final_data <- combined_data[-dim(combined_data)[2]]
#final_data=apply(final_data,2, as.numeric)

####remove  duplication by aggregation
final_data.agg= aggregate(combined_data, list(combined_data$symbol),FUN=mean)
final_data.agg <- final_data.agg[, -2]
final_data_agg <- colnames(final_data.agg)[1] <- "symbol"

# Set the 'symbol' column as row names
row.names(final_data.agg) <- final_data.agg$symbol

# remove symbol column
final_data.agg <- final_data.agg[ , -1]

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

all(colnames(final_data.agg) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(final_data.agg))

########################################
# know the right order indexing
data_idx <- match(rownames(metadata_modified), colnames(final_data.agg))
data_idx

# Reorder the counts data frame to have the sample names in the same order as the metadata data frame
data_ordered  <- final_data.agg[ , data_idx]

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(data_ordered))

#######################################
### we have many conditions so we will begin with sample condition 5 comparing with control 0
meta_5_data <- metadata_modified %>%
  filter(sample_condition == 0 | sample_condition == 5)

# change the 0 to control and 5 to PD stage 5 in the column
colData_5 <- meta_5_data %>% 
  dplyr::mutate(sample_condition = gsub(0, "control", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(5, "PD", sample_condition))

# download the modified metadata to use it in the analysis
write.csv(colData_5, "colData5.csv", row.names = TRUE)

#extract only columns matching with metadata in the same order as row names of meta data 
data_5_ordered <- data_ordered %>% 
  dplyr::select(1, 2, 3, 4, 5, 6, 9, 10, 13, 14,
                17, 18, 19, 22, 23, 26, 27, 29,
                33, 34, 37, 38, 43, 44, 46, 54, 
                55,60, 65, 66, 67, 68, 71, 73, 
                74, 75, 76, 77, 80, 81, 82, 83) 

########## preparing for the DESeq2 analysis for stage 5
# read in the metadata we have adjusted lately
colData_5 <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData5.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData_5) <- colData_5$X

# remove the first column from the meatadata_modified
colData_5 <- colData_5[, -1]

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

########### Begin with the DESeq2 analysis 
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
### we have many conditions so we will begin with sample condition 6 comparing with control 0
meta_6_data <- metadata_modified %>%
  filter(sample_condition == 0 | sample_condition == 6)

# change the 0 to control and 5 to PD stage 5 in the column
colData_6 <- meta_6_data %>% 
  dplyr::mutate(sample_condition = gsub(0, "control", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub(6, "PD", sample_condition))

# download the modified metadata to use it in the analysis
write.csv(colData_6, "colData6.csv", row.names = TRUE)

#extract only columns matching with metadata in the same order as row names of meta data 
data_6_ordered <- data_ordered %>% 
  dplyr::select(5, 6, 7, 11, 16, 21, 22, 31, 34, 35,
                37, 38, 39, 40, 42, 43, 45, 46, 47, 48,
                49, 50, 52, 53, 54, 57, 58, 59, 60, 62,
                63, 65, 66, 67, 70, 71, 72, 73, 74, 75, 
                76, 77, 80, 81, 82, 83)

########## preparing for the DESeq2 analysis for stage 6 data
# read in the metadata we have adjusted lately
colData_6 <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData6.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData_6) <- colData_6$X

# remove the first column from the meatadata_modified
colData_6 <- colData_6[, -1]

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

########### Begin with the DESeq2 analysis 
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
