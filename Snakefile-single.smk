# Snakemake file - input raw reads to generate quant files for analysis in R
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

# Adapters
SE_ADAPTER = config['seq']['SE']
SE_SEQUENCE = config['seq']['trueseq-se']

# Organsim
TRANSCRIPTOME = config['transcriptome']['human']
SPECIES = config['species']['human']
SUF = config["suffix"]

# Samples
sample_links = {"file": "https://osf.io/5daup/download"}
SAMPLES=sample_links.keys()

rule all:
    input:
        # create a new filename for every entry in SAMPLES,
        # replacing {name} with each entry.
        expand("rnaseq/quant/{name}_quant/quant.sf", name=SAMPLES),
        "rnaseq/fastqc/multiqc_report.html",

rule download_reads:
    output: "raw.dir/{sample}.fq.gz" 
    params:
        download_link = lambda wildcards: sample_links[wildcards.sample]
    shell:
        """
        curl -L {params.download_link} -o {output}
        tar xf file.tar.gz
        rm file.tar.gz
        """

rule fastqc:
    input: INPUTDIR + "/{sample}" + SUF
        html = SCRATCH + "/fastqc/{sample}_fastqc.html",
        zip = SCRATCH + "/fastqc/{sample}_fastqc.zip"
    params:""
    wrapper:
        "0.35.2/bio/fastqc"

## quality trim reads and assess with fastqc/multiqc
rule download_trimmomatic_adapter_file:
    output: REFERENCE + SE_ADAPTER
    shell:
        """
        curl -L {SE_SEQUENCE} -o {output}
        """

rule trimmmomatic_se:
    input: 
        reads= INPUTDIR + "/{sample}_" + SUF,
        adapters= REFERENCE + SE_ADAPTER,
    output: 
        reads = SCRATCH + "/trimmed/{sample}_" + "_trim.fastq.gz"
    conda: "rnaseq-env.yml"
    shell:
        """
        trimmomatic SE {input.reads} {output.reads} ILLUMINACLIP:{input.adapters}:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25    
        """

rule fastqc_trim:
    input: SCRATCH + "/trimmed/{sample}_" + "_trim.fastq.gz"
    output:
      html = SCRATCH + "/fastqc/{sample}  _trimmed_fastqc.html",
      zip = SCRATCH + "/fastqc/{sample}_trimmed_fastqc.zip"
  params: ""
  log:
    SCRATCH + "/logs/fastqc/{sample}_trimmed.log"
  wrapper:
      "0.35.2/bio/fastqc"

rule multiqc:
    input:
        raw_qc = expand("{scratch}/fastqc/{sample}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET),
        trim_qc = expand("{scratch}/fastqc/{sample}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET)
    output:
        raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html", 
        raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
        trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", 
        trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"
    shell: 
        """
        multiqc -n multiqc.html {input.raw_qc} #run multiqc
        mv multiqc.html {output.raw_multi_html} #rename html
        mv multiqc_data/multiqc_general_stats.txt {output.raw_multi_stats} #move and rename stats
        rm -rf multiqc_data #clean-up
        #repeat for trimmed data
        multiqc -n multiqc.html {input.trim_qc} #run multiqc
        mv multiqc.html {output.trim_multi_html} #rename html
        mv multiqc_data/multiqc_general_stats.txt {output.trim_multi_stats} #move and rename stats
        rm -rf multiqc_data #clean-up
        """ 

rule download_transcriptome:
    output: REFERENCE + SPECIES
    shell:
        """
        curl -L {TRANSCRIPTOME} -o {output}
        """

rule salmon_index:
    input:
        ref = REFERENCE + SPECIES
    output: directory("output/quant/sc_ensembl_index")
    conda: "env/rnaseq-env.yml"
    shell:
        """
        salmon index --index {output} --transcripts {input} # --type quasi
        """

rule salmon_quant:
    input:
        reads = SCRATCH + "/trimmed/{sample}_" + "_trim.fastq.gz",
        index_dir = "output/quant/sc_ensembl_index"
    output: "rnaseq/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "output/" + wildcards.sample + "_quant
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} --seqBias --gcBias --validateMappings
        """