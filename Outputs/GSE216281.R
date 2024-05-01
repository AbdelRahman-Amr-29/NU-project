setwd("C:/Users/arahm/OneDrive/Desktop/NU/tutorials/GSE216281/Outputs")

# installing packages
#BiocManager::install("biomaRt")
#install.packages("devtools")
#devtools::install_github("stephenturner/annotables")
#BiocManager::install("org.Hs.eg.db")

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
raw_data <- read.table("C:/Users/arahm/OneDrive/Desktop/NU/tutorials/GSE216281/Data/GSE216281_raw_data.txt", header = TRUE, sep = " ")

# Bind the row names with the "gene" column
data <- cbind(ensembl_gene_id = rownames(raw_data), raw_data)

# Remove row names from a raw_data frame
rownames(data) <- NULL

# extract first column
ids <- data[ , 1, drop = FALSE]
head(ids)


# download the first column from the data file contains ids
write.csv(ids, "ensembl_gene_ids.csv", row.names = FALSE)

# read the ensembl_gene_ids file
ensembl_ids <- read.csv("C:/Users/arahm/OneDrive/Desktop/NU/tutorials/GSE216281/Outputs/ensembl_gene_ids.csv", header = TRUE)

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




