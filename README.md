# RNAseq pipeline for paired-sequencing
last updated 2020-01-07

To do:
- script to generate csv to import quant.sf files for analysis in R

Clone repo
```
git clone https://github.com/nasiegel88/rnaseq.git
cd rnaseq
```

Create conda environment
```
conda env create -n rnaseq -f env/rnaseq-env.yml
```

Download practice data by first installing and activating the ```grabseq``` env
```
conda env create -n grabseq -f env/grabseq.yml
# download data from Duclos et al., 2019
# switch to the raw data directory within rnaseq
cd raw.dir
grabseqs sra PRJNA543474
```
* Note this may take a few hours to download

Or do it all at once with the bash script by running the code below
```
bash raw.dir/download_data.sh
```

Next, change the paths in the configurtion file to map from your computer

```
proj_name: name
contact:
  email: jd@some-address.com
  person: John Doe
raw-data: /$PWD/rnaseq/raw.dir
scratch:  /$PWD/rnaseq/scratch
outputDIR: /$PWD/rnaseq/output
metadata: /$PWD/metadata.tsv
ref: /$PWD/rnaseq/refs

# Fastq file suffix, following the read pair designation
suffix: .fastq.gz
# suffix: _001.fastq.gz

# Read pair designations
r1_suf: 1
r2_suf: 2
#r1_suf: R1
#r2_suf: R2
```

