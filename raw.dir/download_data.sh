#!/bin/bash
## Generatate tar.gz file with all files in current directory

# tar cf file.tar.gz *.gz

# Preliminary sequencing run
curl -L https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR9079176/SRR9079176.1 -o file.gz
gunzip *

# wget -c https://osf.io/7uzjb/download -O - | tar -xz