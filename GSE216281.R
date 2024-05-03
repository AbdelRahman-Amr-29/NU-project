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
raw_data <- read.table("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/GSE216281_raw_counts.txt", header = TRUE, sep = " ")

# Bind the row names with the "gene" column
data <- cbind(ensembl_gene_id = rownames(raw_data), raw_data)

# Remove row names from a raw_data frame
rownames(data) <- NULL

# extract first column
ids <- data[ , 1, drop = FALSE]
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
combined_data <- merge(mapper.df, data, by = "ensembl_gene_id", all.x = TRUE)
combined_data <- combined_data[, -1]
combined_data <- combined_data[ ! is.na(combined_data$symbol),]

########################################################
# check the duplicates of symbols column
x <- duplicated(combined_data$symbol)  
sum(x)

####remove  duplication by aggregation
final_data.agg= aggregate(combined_data, list(combined_data$symbol),FUN=mean)
final_data.agg <- final_data.agg[, -2]
final_data_agg <- colnames(final_data.agg)[1] <- "symbol"

# Set the 'symbol' column as row names
row.names(final_data.agg) <- final_data.agg$symbol

# remove symbol column
final_data.agg <- final_data.agg[ , -1]

########################################################
# get metadata
gse <- getGEO(GEO = "GSE216281", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(21, 1, 2, 11, 12, 13) %>%
  dplyr::rename(samples_name = description.1) %>%
  dplyr::rename(gender = characteristics_ch1.1) %>%
  dplyr::rename(sample_condition = characteristics_ch1.2) %>%
  dplyr::rename(age_of_death = characteristics_ch1.3) %>%
  dplyr::mutate(gender = gsub("Sex: ", "", gender)) %>%
  dplyr::mutate(sample_condition = gsub("braak lewy_body_stage: ", "", sample_condition)) %>%
  dplyr::mutate(age_of_death = gsub("age at_death: ", "", age_of_death))

head(metadata_modified)

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
write.csv(colData, "colData.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis
# read in the metadata we have adjusted lately
colData <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE216281/colData.csv")

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

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘sample_R1’
ggplot(data_ordered) +
  geom_histogram(aes(x = sample_R1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘sample_R2’
ggplot(data_ordered) +
  geom_histogram(aes(x = sample_R2), stat = "bin", bins = 200) +
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

# Contrasts
resultsNames(dds)

#MA plot
plotMA(res)

### Transform counts for data visualization
rld <- vst(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="sample_condition")

plotPCA(rld, intgroup="gender")

plotPCA(rld, intgroup="age_of_death")

#### Hierarchical Clustering
### Extract the rlog matrix from the object
rld_mat <- assay(rld) 

### Compute pairwise correlation values
rld_cor <- cor(rld_mat)

## check the output of cor(), make note of the row names and column names
head(rld_cor)

head(colData)

### Plot heatmap using the correlation matrix and the metadata object
pheatmap(rld_cor, annotation = colData)
