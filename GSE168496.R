################################ information about the data
## GSE168496
## Summary:	Purpose: Parkinson’s disease (PD) is one of the most common neurological movement
#disorders. In the incidence and the clinical features of PD, significant gender differences
#have been observed. While males have a higher disease prevalence and are more frequently
#affected by muscle rigidity, females present more often with disabling tremors.
#The molecular mechanisms behind these disease-associated gender differences are still
#largely unknown, and a detailed understanding of the factors involved may
#open new avenues for pharmacological disease modification.
##DOI: 10.1038/s41531-023-00446-8
##########################################################################

# Setting working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE168496")

# installing packages
#BiocManager::install("biomaRt")
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install('airway')
#install.packages("pheatmap")
#install.packages("AnnotationDbi")
#BiocManager::install('EnsDb.Hsapiens.v86')

# calling libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(purrr)
library(ggplot2)
library(DESeq2)
library(airway)
library(pheatmap)
library(biomaRt)

#################################################
# Read txt file & convert it to table
raw_data <- read.table("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE168496/GSE168496_all_samples_preprocessed_data.tsv", header = TRUE, sep = "\t")

# Remove "X" prefix from column names
colnames(raw_data) <- sub("^X", "", colnames(raw_data))

# extract first column
ensembl_trans_id <- raw_data[ , 1, drop = FALSE]
head(ensembl_trans_id)

# download the first column from the data file contains ids
#write.csv(ensembl_trans_id, "ensembl_trans_ids.csv", row.names = FALSE)

# read theensembl_trans_id# read the ensembl_gene_ids file to use it in the next step
ensembl_trans_id <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE168496/ensembl_trans_ids.csv", header = TRUE)

########################################################
# convert gene ids to gene names
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)

gene_symbols <- getBM(attributes = c("ensembl_transcript_id_version", "external_gene_name"),
      filters = "ensembl_transcript_id_version",
      values = ensembl_trans_id$id,
      mart = ensembl.con)

# Rename the columns
names(gene_symbols) <- c("id","symbol")
gene_symbols <- gene_symbols[gene_symbols$symbol != "", , drop = FALSE]
class(gene_symbols)

########################################################
# merging the gene symbols data frame to raw_data
raw_data_merged <- merge(gene_symbols, raw_data, by = "id", all.x = TRUE)
raw_data_merged <- raw_data_merged[, -1]

########################################################
# check the duplicates of symbols column
x <- duplicated(raw_data_merged$symbol)  
sum(x)

####remove  duplication by aggregation
raw_data_agg <- aggregate(raw_data_merged, list(raw_data_merged$symbol),FUN=mean)
raw_data_agg <- raw_data_agg[, -2]
colnames(raw_data_agg)[1] <- "symbol"

# Set the 'symbol' column as row names
row.names(raw_data_agg) <- raw_data_agg$symbol

# remove symbol column
raw_data_agg <- raw_data_agg[ , -1]

# Remove rows where all values are zero
raw_data_agg <- raw_data_agg[rowSums(raw_data_agg != 0, na.rm = TRUE) > 1, ]

########################################################
### get metadata
gse <- getGEO(GEO = "GSE168496", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(1, 10) %>%
  dplyr::rename(samples_name = title) %>%
  dplyr::rename(sample_condition = characteristics_ch1) %>%
  dplyr::mutate(samples_name = gsub("-", ".", samples_name)) %>%
  dplyr::mutate(sample_condition = gsub("disease state: Non-demented control", "control", sample_condition)) %>% 
  dplyr::mutate(sample_condition = gsub("disease state: Parkinson's disease", "PD", sample_condition))

# Remove row names from a metadata_modified frame
rownames(metadata_modified) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(metadata_modified) <- metadata_modified$samples_name

all(colnames(raw_data_agg) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(raw_data_agg))

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(raw_data_agg))

# download the modified metadata to use it in the analysis
#write.csv(metadata_modified, "colData.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis
# read in the metadata we have adjusted lately
colData <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE168496/colData.csv")

# Set the row names of the dataframe to the values in the sample_name column
rownames(colData) <- colData$X

# remove the first column from the meatadata_modified
colData <- colData[, -1]

#check the row names in colData match with column names with data_5_ordered
all(colnames(raw_data_agg) %in% rownames(colData))

#check the order
all(colnames(raw_data_agg) == rownames(colData))

# Convert numeric columns to integer
numeric_columns <- sapply(raw_data_agg, is.numeric)

# Convert numeric columns to integer format
raw_data_agg[, numeric_columns] <- lapply(raw_data_agg[, numeric_columns], as.integer)

# Remove rows where all values are zero
raw_data_agg <- raw_data_agg[rowSums(raw_data_agg != 0, na.rm = TRUE) > 1, ]

# Convert 'sample_condition' column to factor
colData$sample_condition <- factor(colData$sample_condition)
class(colData$sample_condition)

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘1996.105’
ggplot(raw_data_agg) +
  geom_histogram(aes(x = "1996.105"), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘sample_R6’
ggplot(raw_data_agg) +
  geom_histogram(aes(x = sample_R6), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Mean vs Variance
mean_counts <- apply(raw_data_agg[,6:8], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(raw_data_agg[,6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

# plot the mean and variance values against each others
ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

########### Begin with the DESeq2 analysis 
# Construct a DESeqDataSet object 
dds <- DESeqDataSetFromMatrix(countData = raw_data_agg,
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
