######################### Information about the data #######################
#Summary:	Adipose-derived human mesenchymal stem cells (hADSCs) transplantation
#has recently emerged as a promising method in the treatment of Parkinson's disease (PD),
#however, the mechanism underlying has not been fully illustrated.
#In this study, we discovered that hADSCs protected the dopaminergic (DA) neurons
#in a 6-hydroxydopamine(6-OHDA) induced PD mice model. Using a transwell co-culture system,
#we reported that, in 6-OHDA brain slice cultures, hADSCs significantly promoted host
#DA neuronal viability. Within the analysis of hADSCs' exocrine proteins through RNA-seq,
#Human protein cytokine arrays and label-free quantitative proteomics, we identified
#Pentraxin3 (PTX3) as a key extracellular factor in hADSCs secretion environment.
#Moreover, we found that human recombinant Pentraxin3 (rhPTX3) treatment could rescue
#the physiological behaviour of the PD mice in-vivo, as well as prevent DA neurons from
#death and increase the neuronal terminals in the Ventral tegmental area(VTA)+ substantia nigra
#pars compacta (SNc) and striatum (STR) on the PD brain slices in-vitro.
#Furthermore, within testing on the pro-apoptotic markers of PD mice brain following the treatment
#of rhPTX3, we found that rhPTX3 can prevent the apoptosis and degeneration of DA neurons.
#Overall, the current study investigated that PTX3, a hADSCs secreted protein, played a
#potential role in protecting the DA neurons from apoptosis and degeneration in PD progression
#and improving the motor performances in PD mice to give a possible mechanism of how hADSCs works 
#in the cell replacement therapy in PD. Importantly, our study also provided a potential
#translational implications for the development of PTX3-based therapeutics in PD.
#DOI: 10.1096/fj.202100408RR

###########################################################################################

## setting working directory
setwd("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE163176")

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
raw_counts <- read_excel("GSE163176_raw_counts.xlsx")
View(raw_counts)

# make raw_counts data frame
raw_counts <-  data.frame(raw_counts)
class(raw_counts)

##### check the duplicates of gene column
x <- duplicated(raw_counts$gene_symbol)
sum(x)

# aggregation
raw_counts <- aggregate(raw_counts, list(raw_counts$gene_symbol),FUN=mean)
raw_counts <- raw_counts[, -2]
colnames(raw_counts)[1] <- "gene"

# remove rows form the raw_counts that don't have gene name
raw_counts <- raw_counts[c(-1:-76), ]

# make the column gene row names 
rownames(raw_counts) <- raw_counts$gene

# remove the gene column
raw_counts <- raw_counts[, -1]

# Remove rows where all values are zero
raw_counts <- raw_counts[rowSums(raw_counts != 0, na.rm = TRUE) > 1, ]

#####################################################
# get metadata
gse <- getGEO(GEO = "GSE163176", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# get needed columns and modify it
metadata_modified <- metadata %>%
  dplyr::select(1, 11) %>%
  dplyr::rename(treatment = characteristics_ch1.1) %>%
  dplyr::mutate(treatment = gsub("treatment: single culture", "single_culture", treatment)) %>%
  dplyr::mutate(treatment = gsub("treatment: co-culture with 6-OHDA-induced brain slices", "co_cultured", treatment)) 

head(metadata_modified)

# Remove row names from a metadata_modified frame
rownames(metadata_modified) <- NULL

# Set the row names of the dataframe to the values in the sample_name column
rownames(metadata_modified) <- metadata_modified$title

all(colnames(raw_counts) %in% rownames(metadata_modified))
all(rownames(metadata_modified) %in% colnames(raw_counts))

# check if that they are in the same order
all(rownames(metadata_modified) == colnames(raw_counts))

# download the modified metadata to use it in the analysis
#write.csv(metadata_modified, "colData.csv", row.names = TRUE)

########## preparing for the DESeq2 analysis

# read in the metadata we have adjusted lately
colData <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU-Project/GSE163176/colData.csv")

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

# Remove rows where all values are zero
raw_counts <- raw_counts[rowSums(raw_counts != 0, na.rm = TRUE) > 1, ]

# Convert 'diagnosis' column to factor
colData$treatment <- factor(colData$treatment)
class(colData$treatment)

######### Check the distribution of RNA-seq counts 
# plot a histogram of the counts for a single sample, ‘ASC_1’
ggplot(raw_counts) +
  geom_histogram(aes(x = ASC_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# plot a histogram of the counts for a single sample, ‘ASC_PD_1’
ggplot(raw_counts) +
  geom_histogram(aes(x = ASC_PD_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")

# Mean vs Variance
mean_counts <- apply(raw_counts, 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(raw_counts, 1, var)
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
                              design = ~ treatment)

view(counts(dds))

# remove rows that have value less than 10 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

view(counts(dds))

# set the factor level
dds$treatment <-  relevel(dds$treatment, ref = "single_culture")

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
