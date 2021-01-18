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
R1_SUF = str(config["r1_suf"])

# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards(INPUTDIR + "/{sample}_{num}" + SUF)
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

rule all:
    input:
        expand("output/{name}_quant/quant.sf", name=SAMPLE_LIST),
        # fastqc output before trimming
        raw_html = expand("{scratch}/fastqc/{sample}_{num}_fastqc.html", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        raw_zip =expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html",
        raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
        # Trimmed data output
        trimmedData = expand("{scratch}/trimmed/{sample}_{num}_trim.fastq.gz", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS), 
        trim_html = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.html", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", #next change to include proj name
        trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"

rule fastqc:
    input: INPUTDIR + "/{sample}_{num}" + SUF
    output:
        html = SCRATCH + "/fastqc/{sample}_{num}_fastqc.html",
        zip = SCRATCH + "/fastqc/{sample}_{num}_fastqc.zip"
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
        reads= INPUTDIR + "/{sample}_" + R1_SUF + SUF,
        adapters= REFERENCE + SE_ADAPTER,
    output: 
        reads = SCRATCH + "/trimmed/{sample}_" + R1_SUF + "_trim.fastq.gz",
        unpaired = SCRATCH + "/trimmed/{sample}_.unpaired.fastq.gz"
    conda: "rnaseq-env.yml"
    shell:
        """
        trimmomatic SE {input.reads} \
        {output.reads} {output.unpaired} \
        ILLUMINACLIP:{input.adapters}:2:0:15 LEADING:2 TRAILING:2 \
        SLIDINGWINDOW:4:2 MINLEN:25    
        """

rule fastqc_trim:
    input: SCRATCH + "/trimmed/{sample}_{num}" + "_trim.fastq.gz"
    output:
      html = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.html",
      zip = SCRATCH + "/fastqc/{sample}_{num}_trimmed_fastqc.zip"
    params: ""
    log:
      SCRATCH + "/logs/fastqc/{sample}_{num}_trimmed.log"
    wrapper:
        "0.35.2/bio/fastqc"

rule multiqc:
    input:
        raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS)
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
        reads = SCRATCH + "/trimmed/{sample}_" + R1_SUF + "_trim.fastq.gz",
        index_dir = "output/quant/sc_ensembl_index"
    output: "output/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "output/" + wildcards.sample + "_quant"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} --seqBias --gcBias --validateMappings
        """