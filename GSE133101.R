############### Information about the data #####################
#Summary:	total RNA sequencing of Parkinson's disease 3 regions of human brain
#Overall design:	ribo-depleted RNA libraries of Parkinson's disease
#and control samples of 3 brain regions: Amygdala (AMG), Substantia nigra (SN)
#and medil temporal gyrus (MTG)
#DOI: 10.15252/emmm.201911942
#################################################################

## setting the working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE133101")

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
raw_counts <- read_excel("amygdala_raw.xlsx")
View(raw_counts)

# make raw_counts data frame
raw_counts <-  data.frame(raw_counts)
class(raw_counts)

##### check the duplicates of gene column
x <- duplicated(raw_counts$symbol)
sum(x)

# aggregation
raw_counts <- aggregate(raw_counts, list(raw_counts$symbol),FUN=mean)
raw_counts <- raw_counts[, -2]
colnames(raw_counts)[1] <- "gene"

# make the column gene row names 
rownames(raw_counts) <- raw_counts$gene

# remove the gene column
raw_counts <- raw_counts[, -1]

# Remove rows where all values are zero
raw_counts <- raw_counts[rowSums(raw_counts != 0, na.rm = TRUE) > 1, ]

# Remove rows with NA values
raw_counts <- na.omit(raw_counts)

#####################################################
# get metadata
gse <- getGEO(GEO = "GSE133101", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(1, 10, 11) %>%
  dplyr::rename(sample_condition = characteristics_ch1) %>%
  dplyr::rename(tissue = characteristics_ch1.1) %>% 
  dplyr::mutate(sample_condition = gsub("sample_condition: Control", "control", sample_condition)) %>%
  dplyr::mutate(sample_condition = gsub("sample_condition: Parkinson's disease", "PD", sample_condition)) %>%
  dplyr::mutate(tissue = gsub("tissue: ", "", tissue))  

metadata_amg <- metadata_modified %>%
  dplyr::filter(tissue == "amygdala (AMG)")

head(metadata_amg)

#check the row names in colData_amg match with column names with raw_counts
all(colnames(raw_counts) %in% rownames(metadata_amg))
all(rownames(metadata_amg) %in% colnames(raw_counts))

# check if that they are in the same order
all(rownames(metadata_amg) == colnames(raw_counts))

# download the modified metadata to use it in the analysis
#write.csv(metadata_amg, "colData_amg_amg.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis
# read in the metadata we have adjusted lately
colData_amg <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE133101/colData_amg.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData_amg) <- colData_amg$X

# remove the first column from the meatadata_modified
colData_amg <- colData_amg[, -1]

#check the row names in colData_amg match with column names with raw_counts
all(colnames(raw_counts) %in% rownames(colData_amg))

#check the order
all(colnames(raw_counts) == rownames(colData_amg))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_counts, is.numeric)

# Convert numeric columns to integer format
raw_counts[, numeric_columns] <- lapply(raw_counts[, numeric_columns], as.integer)

# Convert 'sample_condition' column to factor
colData_amg$sample_condition <- factor(colData_amg$sample_condition)
class(colData_amg$sample_condition)

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘GSM3900032’
ggplot(raw_counts) +
  geom_histogram(aes(x = GSM3900032), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘GSM3900040’
ggplot(raw_counts) +
  geom_histogram(aes(x = GSM3900040), stat = "bin", bins = 200) +
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
dds_amg <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = colData_amg,
                              design = ~ sample_condition)

view(counts(dds_amg))

# remove rows that have value less than 10 
keep <- rowSums(counts(dds_amg)) >= 10
dds_amg <- dds_amg[keep, ]

view(counts(dds_amg))

# set the factor level
dds_amg$sample_condition <-  relevel(dds_amg$sample_condition, ref = "control")

##### Run DESeq
dds_amg <- DESeq(dds_amg)

res_amg <- results(dds_amg)
res_amg

# Explore results
summary(res_amg)

res_amg0.05 <- results(dds_amg, alpha = 0.05)
summary(res_amg0.05)

res_amg0.25 <- results(dds_amg, alpha = 0.25)
summary(res_amg0.25)
