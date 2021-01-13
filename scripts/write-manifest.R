
# Required libraries
library(reshape2);library(tidyverse)

# Import all files in currect directory
paths <- as.data.frame(list.files(pattern = "_quant", full.names = TRUE))
colnames(paths)[1]<-"QUANT"

# Get local path & add to dataframe
path_to_files <- getwd()
paths$PATH <- path_to_files

# Extract sample ID
paths_run <- paths %>% mutate(SAMPLEID = str_replace(QUANT, ""))
## ^See wildcard options on this line to modify how R script pulls out your sample IDs from fastq files

paths_run$SAMPLEID <- gsub("./","", paths_run$SAMPLEID)
paths_run$QUANT <- gsub("./","", paths_run$QUANT)

# Write full path
paths_run$FULL_PATH <- paste(paths_run$PATH, paths_run$QUANT, sep="/")

# Write manifest:
manifest_orig <- data.frame(paths_run$SAMPLEID, paths_run$FULL_PATH)
colnames(manifest_orig)[1:2]<-c("names","files")

# srtings
manifest_orig$names <- str_remove(manifest_orig$names, "_quant")
manifest_orig$files <- str_replace(manifest_orig$files, "_quant","_quant/quant.sf")

# input the treatments here
manifest_orig$condition <- c("NS", "NS", "NS", "NS", "CS", "CS", "CS", "CS")


# Write output as a manifest file
write.table(manifest_orig, "rnaseq_samples.csv",quote=FALSE,col.names=TRUE,row.names=FALSE,sep=",")
