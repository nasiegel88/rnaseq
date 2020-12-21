# Use glob statement to find all samples in 'raw_data' directory
## Wildcard '{num}' must be equivalent to 'R1' or '1', meaning the read pair designation.
SAMPLE_LIST,NUMS = glob_wildcards("raw.dir/{sample}_{num}.fastq.gz")
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
SET_NUMS = set(NUMS)

# the sample names are dictionary keys in sample_links. extract them to a list we can use below
SAMPLES=sample_links.keys()
SET_NUMS = set(NUMS)

rule all:
    input:
        # create a new filename for every entry in SAMPLES,
        # replacing {name} with each entry.
        expand("rnaseq/quant/{name}_{num}_quant/quant.sf", name=SAMPLES),
        "rnaseq/fastqc/multiqc_report.html",

# download human rna-seq data from Duclose et al, 2019 study
rule download_reads:
    output: 
        r1 = expand("rnaseq/quality/{sample}_{num}.qc.fq.gz", sample=SAMPLES, num=SET_NUMS),
        r2 = expand("rnaseq/quality/{sample}_{num}.qc.fq.gz", sample=SAMPLES, num=SET_NUMS)
    params:
        download_link = lambda wildcards: sample_links[wildcards.sample]
    shell:
        """
        curl -L {params.download_link} -o {output}
        """

rule fastqc_raw:
    input: 
        r1 = expand("rnaseq/quality/{sample}_{num}.qc.fq.gz", sample=SAMPLES, num=SET_NUMS),
        r2 = expand("rnaseq/quality/{sample}_{num}.qc.fq.gz", sample=SAMPLES, num=SET_NUMS)
    output: 
        r1_fastq = "rnaseq/raw_data/fastqc/{sample}_{num}.fastqc.html",
        r2_fastq = "rnaseq/raw_data/fastqc/{sample}_{num}-fastqc.html"
    params:
        outdir="rnaseq/raw_data/fastqc"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        fastqc {input.r1} --outdir {output.r1_fastq}
        fastqc {input.r2} --outdir {output.r2_fastq}
        """

## quality trim reads and assess with fastqc/multiqc
rule download_trimmomatic_adapter_file:
    output: "TruSeq2-SE.fa"
    shell:
        """
        curl -L https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE-2.fa  -o {output}
        """

rule quality_trim:
    input: 
        r1 = expand("rnaseq/quality/{sample}_{num}.qc.fq.gz", sample=SAMPLES, num=SET_NUMS),
        r2 = expand("rnaseq/quality/{sample}_{num}.qc.fq.gz", sample=SAMPLES, num=SET_NUMS),
        adapters="TruSeq3-PE-2.fa"
    output: 
        "rnaseq/quality/{sample}.qc.fq.gz"
        trimmedf = "rnaseq/quality/{sample}_{num}.trimmed.fq.gz",
        trimmedr = "rnaseq/quality/{sample}_{num}.trimmed.fq.gz",
        untrimmedf = "rnaseq/quality/{sample}_{num}.untrimmed.fq.gz",
        untrimmedr = "rnaseq/quality/{sample}_{num}.untrimmed.fq.gz"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        trimmomatic SE {input.reads} {output} ILLUMINACLIP:{input.adapters}:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:20 MINLEN:35
        
        trimmomatic PE -threads 4 {input.r1} {input.r2}  \
              {SRR_1056_1.trimmed.fastq} {SRR_1056_1untrimmed.fastq} \
              {SRR_1056_2.trimmed.fastq} {SRR_1056_2untrimmed.fastq} \
              ILLUMINACLIP:{input.adapters} SLIDINGWINDOW:4:2
        """

rule fastqc_trimmed:
    input: 
        trimmedf = expand("rnaseq/quality/{sample}_{num}.trimmed.fq.gz", sample=SAMPLES, num=SET_NUMS),
        trimmedr = expand("rnaseq/quality/{sample}_{num}.trimmed.fq.gz", sample=SAMPLES, num=SET_NUMS),
    output: 
        "rnaseq/quality/fastqc/{sample}.qc_fastqc.html"
    params:
        outdir="rnaseq/quality/fastqc"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        fastqc {input.trimmedf} --outdir {params.outdir}
        fastqc {input.trimmedr} --outdir {params.outdir}
        """

rule multiqc:
    input: 
        raw=expand("rnaseq/raw_data/fastqc/{sample}_fastqc.html", sample=SAMPLES, num=SET_NUMS),
        trimmed=expand("rnaseq/quality/fastqc/{sample}.qc_fastqc.html", sample=SAMPLES, num=SET_NUMS)
    output: 
        "rnaseq/fastqc/multiqc_report.html"
    params:
        raw_dir="rnaseq/raw_data/fastqc",
        trimmed_dir="rnaseq/raw_data/fastqc",
    conda: 
        "rnaseq-env.yml"
    shell:
        """   
        multiqc {params.raw_dir} {params.trimmed_dir} --filename {output}
        """

rule multiqc:
    input:
        raw_qc = expand("rnaseq/raw_data/fastqc/{sample}_{num}_fastqc.zip", sample=SAMPLES, num=SET_NUMS),
        trim_qc = expand("rnaseq/raw_data/fastqc/{sample}_{num}_trimmed_fastqc.zip", sample=SAMPLES, num=SET_NUMS)  
    output: 
        "rnaseq/fastqc/multiqc_report.html"
    params:
        raw_dir="rnaseq/raw_data/fastqc",
        trimmed_dir="rnaseq/raw_data/fastqc"
    conda: 
        "rnaseq-env.yml"
    shell:
        """   
        multiqc {params.raw_dir} {params.trimmed_dir} --filename {output}
        """

### download and index the human transcriptome ###
rule download_human_transcriptome:
    output: 
        "rnaseq/reference/Mus_musculus.GRCm38.cdna.all.fa.gz" 
    shell:
        """
        curl -L ftp://ftp.ensembl.org/pub/release-101/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz -o {output}
        """

rule salmon_index:
    input: 
        "rnaseq/reference/Mus_musculus.GRCm38.cdna.all.fa.gz" 
    output: 
        "rnaseq/quant/sc_ensembl_index"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        salmon index --index {output} --transcripts {input} # --type quasi -k 31
        """

### quantify reads with salmon
rule salmon_quantify:
    input:
        r1 = "rnaseq/quality/{sample}_1.qc.fq.gz",
        r2 = "rnaseq/quality/{sample}_2.qc.fq.gz",
        index_dir = "rnaseq/quant/sc_ensembl_index"
    output:
        "rnaseq/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/quant/" + wildcards.sample + "_quant"
    conda: 
        "rnaseq-env.yml"
    shell:
        """ 
        salmon quant -i {input.index_dir} -l A -p 6 --validateMappings \
         --gcBias --numGibbsSamples 20 -o {params.outdir} \
         -1 {input.r1} -2 {input.r2}
        """