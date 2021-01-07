#!/bin/bash

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# activate conda environment
conda activate grabseq

# navigate to the raw.dir directory
cd raw.dir

# download files
grabseqs sra PRJNA543474
