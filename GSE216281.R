# setting the working directory
getwd()
setwd("C:/Users/arahm/OneDrive/Desktop/NU/tutorials/GSE216281/GSE216281_counts.txt")

# installing packages
BiocManager::install("biomaRt")
install.packages("devtools")
devtools::install_github("stephenturner/annotables")
BiocManager::install("org.Hs.eg.db")

# calling libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(ggplot2)
library(biomaRt)
library(annotables)
library(org.Hs.eg.db)

#################################################
# Read txt file & convert it to table
raw_data <- read.table("C:/Users/arahm/OneDrive/Desktop/NU/tutorials/GSE216281/GSE216281_counts.txt/GSE216281_counts.txt", header = TRUE, sep = " ")

# Bind the row names with the "gene" column
data <- cbind(ensembl_gene_id = rownames(raw_data), raw_data)

# Remove row names from a raw_data frame
rownames(data) <- NULL

# extract first column
ids <- data[ , 1, drop = FALSE]
head(ids)


# download the first column from the data file contains ids
write.csv(ids, "ensembl_gene_ids.csv", row.names = FALSE)


########################################################
# ways of get the gene symbol but will not used 
# input list of ensembl ID's
#ensembl_ids <- read.delim("C:/Users/arahm/OneDrive/Desktop/NU/tutorials/GSE216281/GSE216281_counts.txt/ensembl_gene_ids.csv")

#listEnsembl()
#ensembl <- useEnsembl(biomart = "genes")
#datasets <- listDatasets(ensembl)

#ensembl.con <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#attr <- listAttributes(ensembl.con)
#filters <- listFilters(ensembl.con)


#gene_name_id_ensembl <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
      #filters = "ensembl_gene_id",
      #values = ensembl_ids$ensembl_gene_id,
      #mart = ensembl.con)
########################################################
#annot_way <- grch38 %>%
  #filter(ensgene %in% ensembl_ids$ensembl_gene_id)
########################################################
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

### yes .. why ? transcripts?  solutions : aggregation
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


#########################################################
# get medtadata
gse <- getGEO(GEO = "GSE216281", GSEMatrix = TRUE)

gse

metadata <- pData(phenoData(gse[[1]]))

rownames(metadata) <- NULL

metadata_subset <- metadata %>%
  select(21, 2, 1, 43, 44, 47)
  head()

# make description.1 column (samples) row names
row_name_column <- metadata_subset$description.1
rownames(metadata_subset) <- row_name_column

all(colnames(raw_data) %in% rownames(metadata_subset))

grch38

