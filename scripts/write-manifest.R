# Read in SRA sample IDs and list of raw fastq files to generate a manifest file to import quantsf files

# Required libraries
library(reshape2);library(tidyverse)

# Import all files in currect directory
paths <- as.data.frame(list.files(pattern = "_1.fastq.gz", full.names = TRUE))
colnames(paths)[1]<-"FASTQ"

# Get local path & add to dataframe
path_to_files <- getwd()
paths$PATH <- path_to_files

# Extract sample ID
paths_run <- paths %>% mutate(SAMPLEID = str_replace(FASTQ, "(\\w*?)_(\\w+).fastq.gz","\\1"))
## ^See wildcard options on this line to modify how R script pulls out your sample IDs from fastq files

paths_run$SAMPLEID <- gsub("./","", paths_run$SAMPLEID)
paths_run$FASTQ <- gsub("./","", paths_run$FASTQ)

# Write full path
paths_run$FULL_PATH <- paste(paths_run$PATH, paths_run$FASTQ, sep="/")

# Write manifest:
manifest_orig <- data.frame(paths_run$SAMPLEID, paths_run$FULL_PATH)
colnames(manifest_orig)[1:2]<-c("names","files")

# If R1_001.fastq.gz or _1.fastq.gz ending, fill new column with forward, else reverse
manifest_orig$condition <- c("NS", "NS", "NS", "NS", "NS", "NS", "CS", "CS", "CS", "CS", "CS", "CS")
manifest_orig$files <- str_replace(manifest_orig$files, "_1.fastq.gz","quant.sf")


# Write output as a manifest file
write.table(manifest_orig, "rnaseq_samples.csv",quote=FALSE,col.names=TRUE,row.names=FALSE,sep=",")
