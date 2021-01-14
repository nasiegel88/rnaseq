# Bulk RNAseq pipeline for paired-end and single-end reads
last updated 2020-01-07

Author: Noah Siegel

Clone repo
```
git clone https://github.com/nasiegel88/rnaseq.git
cd rnaseq
```

Create conda environment
```
conda env create -n rnaseq -f env/rnaseq-env.yml
```
This environment has everything needed to run snakemake, perform quality control of reads, and quantify reads based on what matches a reference genome.

Download practice data from [Himes at et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24926665/) by first installing and activating the ```grabseq``` env. Note, only do this for paired sequencing data. If you want to test out the snakefiles for single-end sequencing you either need to find different practice data or simply delete the reverse reads from Himes et al., 2014 by using something like ```rm *_2.fastq.gz``` in the raw.dir directory.
```
conda env create -n grabseq -f env/grabseq.yml
conda activate grabseq
# switch to the raw data directory within rnaseq
cd raw.dir
grabseqs sra PRJNA229998
```
* Note this may take a few hours to download

Or do it all at once with the bash script by running the code below...
```
bash raw.dir/download_data.sh
```

Next, change the paths in the configurtion file to map to your computer. I am only using TruSeq-PE adapter but it possible to run this workflow by substituting for different adapter sequences. Currently, the configuration file supports both paired and single end TruSeq adapters

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
  PE: TruSeq2-PE.fa
  trueseq-pe: https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-PE.fa
  SE: TruSeq2-SE.fa
  trueseq-se: https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-SE.fa

# Downloads
transcriptome:
  human: ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
  mouse: ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz


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

Follow variable assigments based on the kind of reads you have. 

Paired end reads:
```
# Adapters
ADAPTER = config['seq']['PE']
SEQUENCE = config['seq']['trueseq-pe']
```

Single end reads:
```
# Adapters
SE_ADAPTER = config['seq']['SE']
SE_SEQUENCE = config['seq']['trueseq-se']
```

Specify which organsim transcriptome you want to align too. Currently mouse and human are the only transcriptomes listed in the config.yml. Different reference transcriptiome can be found at the [enembl](https://uswest.ensembl.org/info/data/ftp/index.html) database.

```
# Organisim
TRANSCRIPTOME = config['transcriptome']['human']
SPECIES = config['species']['human']
```

If working with mouse sample adjusted the above code to...

```
# Organisim
TRANSCRIPTOME = config['transcriptome']['mouse']
SPECIES = config['species']['mouse']
```

Once every thing is set up activate the environment using,
```conda activate rnaseq```
You are ready to run Snakemake.

Paired-end: ```snakemake -s Snakefile-paired.smk -j 2```

Single-end: ```snakemake -s Snakefile-single.smk -j 2```
