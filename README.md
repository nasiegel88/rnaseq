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

Next, change the paths in the configurtion file to map from your computer. I am only using TruSeq-PE adapter but it possible to run this workflow by substituting for different adapter sequences.

```
proj_name: name
contact:
  email: jd@some-address.com
  person: John Doe
raw-data: /$PWD/rnaseq/raw.dir
scratch:  /$PWD/rnaseq/scratch
outputDIR: /$PWD/rnaseq/output
metadata: /$PWD/metadata.tsv
ref: /$PWD/rnaseq/refs/

# Organism
species:
  mouse: Mus_musculus.GRCm38.cdna.all.fa.gz
  human: Homo_sapiens.GRCh38.cdna.all.fa.gz

# Adapters
seq:
  name: TruSeq2-PE.fa
  trueseq: curl -L https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-PE.fa

# Downloads
transcriptome:
  human: curl -L ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  mouse: curl -L ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz

# Fastq file suffix, following the read pair designation
suffix: .fastq.gz
# suffix: _001.fastq.gz

# Read pair designations
r1_suf: 1
r2_suf: 2
#r1_suf: R1
#r2_suf: R2
```

Lastly, some changes will need to be made to the snakefile...

```
# Snakemake file
configfile: "config.yaml"
import io 
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#
PROJ = config["proj_name"]
INPUTDIR = config["raw-data"]
SCRATCH = config["scratch"]
REFERENCE = config["ref"]
OUTPUTDIR = config["outputDIR"]
```

If using TruSeq adapters remove dont't change anything here.

```
# Adapters
ADAPTER = config['seq']['name']
SEQUENCE = config['seq']['trueseq']
```

Specify which organsim transcriptome you want to align too. Currently mouse and human are the only transcriptomes listed in the config.yml. Different reference transcriptiome can be found at the [enembl](https://uswest.ensembl.org/info/data/ftp/index.html) database.

```
# Organisim
TRANSCRIPTOME = config['transcriptome']['human']
SPECIES = config['species']['human']
```