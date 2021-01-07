# Snakemake file - input raw fastq reads to generate ASVs
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
REFERENCE = config["ref"]["human"]
OUTPUTDIR = config["outputDIR"]

SUF = config["suffix"]
R1_SUF = str(config["r1_suf"])
R2_SUF = str(config["r2_suf"])

# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards(INPUTDIR + "/{sample}_{num}" + SUF)
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

# Uncomment to diagnose if snakemake is reading in wildcard variables properly
#print(SAMPLE_SET)
#print(SET_NUMS)
#print(SUF)
#print(R1_SUF)

# Create final manifest file for qiime2
MANIFEST = pd.read_csv(config["manifest"]) #Import manifest
MANIFEST['filename'] = MANIFEST['absolute-filepath'].str.split('/').str[-1] #add new column with only file name
MANIFEST.rename(columns = {'absolute-filepath':'rawreads-filepath'}, inplace = True)
PATH_TRIMMED = "trimmed" # name of directory with trimmed reads
NEWPATH = os.path.join(SCRATCH, PATH_TRIMMED)
MANIFEST['filename'] = MANIFEST['filename'].str.replace(SUF, "_trim.fastq.gz")
MANIFEST['absolute-filepath'] = NEWPATH+ "/" + MANIFEST['filename']    
MANIFEST[['sample-id','absolute-filepath','direction']].set_index('sample-id').to_csv('manifest-trimmed.txt')
MANIFEST = config["manifest"]

rule all:
    input:
        expand("rnaseq/quant/{name}_quant/quant.sf", name=SAMPLES
        # fastqc output before trimming
        raw_html = expand("{scratch}/fastqc/{sample}_{num}_fastqc.html", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        raw_zip = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch = SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        raw_multi_html = SCRATCH + "/fastqc/raw_multiqc.html",
        raw_multi_stats = SCRATCH + "/fastqc/raw_multiqc_general_stats.txt",
        # Trimmed data output
        trimmedData = expand("{scratch}/trimmed/{sample}_{num}_trim.fastq.gz", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS), 
        trim_html = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.html", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        raw_qc = expand("{scratch}/fastqc/{sample}_{num}_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        trim_qc = expand("{scratch}/fastqc/{sample}_{num}_trimmed_fastqc.zip", scratch= SCRATCH, sample=SAMPLE_SET, num=SET_NUMS),
        trim_multi_html = SCRATCH + "/fastqc/trimmed_multiqc.html", #next change to include proj name
        trim_multi_stats = SCRATCH + "/fastqc/trimmed_multiqc_general_stats.txt"




sample_links = {"GSM3773108": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079176/1061_TGACCA_combined_R2_001.fastq.gz",
                "GSM3773109": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079177/1067_ACAGTG_combined_R2_001.fastq.gz",
                "GSM3773111": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079178/1069_GCCAAT_combined_R2_001.fastq.gz",
                "GSM3773112": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079179/1074_ACTTGA_combined_R2_001.fastq.gz",
                "GSM3773114": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079180/1076_GATCAG_combined_R2_001.fastq.gz",
                "GSM3773116": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079181/1079_TAGCTT_combined_R2_001.fastq.gz",
                "GSM3773117": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079182/1017_ATCACG_combined_R2_001.fastq.gz",
                "GSM3773119": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079183/1027_CGATGT_combined_R2_001.fastq.gz",
                "GSM3773121": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079184/1059_TTAGGC_combined_R2_001.fastq.gz",
                "GSM3773122": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079187/1085_CTTGTA_combined_R2_001.fastq.gz",
                "GSM3773124": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079185/1072_CAGATC_combined_R2_001.fastq.gz",
                "GSM3773125": "https://sra-pub-src-2.s3.amazonaws.com/SRR9079186/1082_GGCTAC_combined_R2_001.fastq.gz"}

# the sample names are dictionary keys in sample_links. extract them to a list we can use below
SAMPLES=sample_links.keys()

rule download_reads:
    output: "rnaseq/raw_data/{sample}.fq.gz" 
    params:
        download_link = lambda wildcards: sample_links[wildcards.sample]
    shell:
        """
        curl -L {params.download_link} -o {output}
        """

rule fastqc:
    input:    
        INPUTDIR + "/{sample}_{num}" + SUF
    output:
        html = SCRATCH + "/fastqc/{sample}_{num}_fastqc.html",
        zip = SCRATCH + "/fastqc/{sample}_{num}_fastqc.zip"
    params: ""
    log:
        SCRATCH + "/logs/fastqc/{sample}_{num}.log"
    wrapper:
        "0.35.2/bio/fastqc

## quality trim reads and assess with fastqc/multiqc
rule download_trimmomatic_adapter_file:
    output: "TruSeq2-SE.fa"
    shell:
        """
        curl -L https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-SE.fa -o {output}
        """

rule trimmomatic_pe:
    input:
        r1 = INPUTDIR + "/{sample}_" + R1_SUF + SUF,
        r2 = INPUTDIR + "/{sample}_" + R2_SUF + SUF
    output:
        r1 = SCRATCH + "/trimmed/{sample}_" + R1_SUF + "_trim.fastq.gz",
        r2 = SCRATCH + "/trimmed/{sample}_" + R2_SUF + "_trim.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired = SCRATCH + "/trimmed/{sample}_1.unpaired.fastq.gz",
        r2_unpaired = SCRATCH + "/trimmed/{sample}_2.unpaired.fastq.gz"
    log:
        SCRATCH + "/trimmed/logs/trimmomatic/{sample}.log"
    params:
        trimmer = ["LEADING:2", "TRAILING:2", "SLIDINGWINDOW:4:2", "MINLEN:25"],
        extra = ""
    wrapper:
        "0.35.2/bio/trimmomatic/pe"

rule fastqc_trim:
  input:
    SCRATCH + "/trimmed/{sample}_{num}_trim.fastq.gz"
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

  conda:
   "../envs/multiqc-env.yaml"
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

### download and index the human transcriptome ###
rule download_human_transcriptome:
    output: "rnaseq/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz" 
    shell:
        """
        curl -L ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o {output}
        """

rule salmon_index:
    input:
        ref = REFERENCE + "Homo_sapiens.GRCh38.cdna.all.fa.gz" 
    output: "rnaseq/quant/sc_ensembl_index"
    conda: "env/rnaseq-env.yml"
    shell:
        """
        salmon index --index {output} --transcripts {input} # --type quasi -k 21
        """

### quantify reads with salmon
rule salmon_quantify:
    input:
        reads="rnaseq/quality/{sample}.qc.fq.gz",
        index_dir="rnaseq/quant/sc_ensembl_index"
    output: "rnaseq/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/quant/" + wildcards.sample + "_quant"
    conda: "env/rnaseq-env.yml"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} --seqBias --useVBOpt --validateMappings
        """
        

rule salmon_quant:
    input:
        r1 = SCRATCH + "/trimmed/{sample}_" + R1_SUF + "_trim.fastq.gz",
        r2 = SCRATCH + "/trimmed/{sample}_" + R2_SUF + "_trim.fastq.gz",
        index_dir = "rnaseq/quant/sc_ensembl_index" 
    output: "rnaseq/output/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/output/" + wildcards.sample + "_quant"
    conda: "env/rnaseq-env.yml"
    shell:
        "{SALMON} quant -i {input.index_dir} -l A -p 6 --validateMappings \
         --gcBias --numGibbsSamples 20 -o {params.dir} \
         -1 {input.r1} -2 {input.r2}"
