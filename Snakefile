sample_links = {"ETS_NEO": "https://osf.io/tyce4/download",
                 "FA_NEO": "https://osf.io/whx3s/download"}

# the sample names are dictionary keys in sample_links. extract them to a list we can use below
SAMPLES=sample_links.keys()

rule all:
    input:
        # create a new filename for every entry in SAMPLES,
        # replacing {name} with each entry.
        expand("rnaseq/quant/{name}_quant/quant.sf", name=SAMPLES),
        "rnaseq/fastqc/multiqc_report.html",

# download human rna-seq data from Duclose et al, 2019 study
rule download_reads:
    output: 
        "rnaseq/raw_data/{sample}.fq.gz" 
    params:
        download_link = lambda wildcards: sample_links[wildcards.sample]
    shell:
        """
        curl -L {params.download_link} -o {output}
        """

rule fastqc_raw:
    input: 
        "rnaseq/raw_data/{sample}.fq.gz"
    output: 
        "rnaseq/raw_data/fastqc/{sample}_fastqc.html"
    params:
        outdir="rnaseq/raw_data/fastqc"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        fastqc {input} --outdir {params.outdir}
        """

## quality trim reads and assess with fastqc/multiqc
rule download_trimmomatic_adapter_file:
    output: "TruSeq2-SE.fa"
    shell:
        """
        curl -L https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-SE.fa -o {output}
        """

rule quality_trim:
    input: 
        reads="rnaseq/raw_data/{sample}.fq.gz",
        adapters="TruSeq2-SE.fa",
    output: 
        "rnaseq/quality/{sample}.qc.fq.gz"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        trimmomatic SE {input.reads} {output} ILLUMINACLIP:{input.adapters}:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:20 MINLEN:35    
        """

rule fastqc_trimmed:
    input: 
        "rnaseq/quality/{sample}.qc.fq.gz"
    output: 
        "rnaseq/quality/fastqc/{sample}.qc_fastqc.html"
    params:
        outdir="rnaseq/quality/fastqc"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        fastqc {input} --outdir {params.outdir}
        """

rule multiqc:
    input: 
        raw=expand("rnaseq/raw_data/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        trimmed=expand("rnaseq/quality/fastqc/{sample}.qc_fastqc.html", sample=SAMPLES)
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
        reads="rnaseq/quality/{sample}.qc.fq.gz",
        index_dir="rnaseq/quant/sc_ensembl_index"
    output: "rnaseq/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/quant/" + wildcards.sample + "_quant"
    conda: 
        "rnaseq-env.yml"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} --seqBias --useVBOpt --validateMappings
        """

