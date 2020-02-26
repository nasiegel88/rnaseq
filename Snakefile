sample_links = {"ERR458493": "https://osf.io/5daup/download",
                "ERR458494":"https://osf.io/8rvh5/download",
                 "ERR458495":"https://osf.io/2wvn3/download",
                 "ERR458500":"https://osf.io/xju4a/download",
                 "ERR458501": "https://osf.io/nmqe6/download",
                 "ERR458502": "https://osf.io/qfsze/download"}

# the sample names are dictionary keys in sample_links. extract them to a list we can use below
SAMPLES=sample_links.keys()

rule all:
    input:
        # create a new filename for every entry in SAMPLES,
        # replacing {name} with each entry.
        expand("rnaseq/quant/{name}_quant/quant.sf", name=SAMPLES),
        "rnaseq/fastqc/multiqc_report.html",
        "rnaseq/diffexp/rnaseq_samples.csv"

# download yeast rna-seq data from Schurch et al, 2016 study
rule download_reads:
    output: "rnaseq/raw_data/{sample}.fq.gz" 
    params:
        download_link = lambda wildcards: sample_links[wildcards.sample]
    shell:
        """
        curl -L {params.download_link} -o {output}
        """

rule fastqc_raw:
    input: "rnaseq/raw_data/{sample}.fq.gz"
    output: "rnaseq/raw_data/fastqc/{sample}_fastqc.html"
    params:
        outdir="rnaseq/raw_data/fastqc"
    conda: "rnaseq-env.yml"
    shell:
        """
        fastqc {input} --outdir {params.outdir}
        """

## quality trim reads and assess with fastqc/multiqc
rule download_trimmomatic_adapter_file:
    output: "TruSeq2-SE.fa"
    shell:
        """
        curl -L https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-SE.fa -o {output}
        """

rule quality_trim:
    input: 
        reads="rnaseq/raw_data/{sample}.fq.gz",
        adapters="TruSeq2-SE.fa",
    output: "rnaseq/quality/{sample}.qc.fq.gz"
    conda: "rnaseq-env.yml"
    shell:
        """
        trimmomatic SE {input.reads} {output} ILLUMINACLIP:{input.adapters}:2:0:15 LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2 MINLEN:25    
        """

rule fastqc_trimmed:
    input: "rnaseq/quality/{sample}.qc.fq.gz"
    output: "rnaseq/quality/fastqc/{sample}.qc_fastqc.html"
    params:
        outdir="rnaseq/quality/fastqc"
    conda: "rnaseq-env.yml"
    shell:
        """
        fastqc {input} --outdir {params.outdir}
        """

rule multiqc:
    input: 
        raw=expand("rnaseq/raw_data/fastqc/{sample}_fastqc.html", sample=SAMPLES),
        trimmed=expand("rnaseq/quality/fastqc/{sample}.qc_fastqc.html", sample=SAMPLES)
    output: "rnaseq/fastqc/multiqc_report.html"
    params:
        raw_dir="rnaseq/raw_data/fastqc",
        trimmed_dir="rnaseq/raw_data/fastqc",
    conda: "rnaseq-env.yml"
    shell:
        """
        multiqc {params.raw_dir} {params.trimmed_dir} --filename {output}
        """

### download and index the yeast transcriptome ###
#rule download_yeast_transcriptome:
#    output: "rnaseq/reference/GCA_000146045.2_R64_rna_from_genomic.fna.gz"
#    shell:
        #curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/146/045/GCA_000146045.2_R64/GCA_000146045.2_R64_rna_from_genomic.fna.gz
        #mv GCA_000146045.2_R64_rna_from_genomic.fna.gz {output}
#        """
#        curl -L https://osf.io/yyyyyyy/download -o {output}
#        """

rule salmon_index:
    input:  "rnaseq/reference/GCA_000146045.2_R64_rna_from_genomic.fna.gz"
    output: directory("rnaseq/quant/sc_index")
    conda: "rnaseq-env.yml"
    shell:
        """
        salmon index --index {output} --transcripts {input} # --type quasi
        """

### quantify reads with salmon
rule salmon_quantify:
    input:
        reads="rnaseq/quality/{sample}.qc.fq.gz",
        index_dir="rnaseq/quant/sc_index"
    output: "rnaseq/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/quant/" + wildcards.sample + "_quant"
    conda: "rnaseq-env.yml"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} --seqBias --gcBias --validateMappings
        """

# download sample information for differential expression analysis
rule download_sample_condition_csv: 
    output: "rnaseq/diffexp/rnaseq_samples.csv"
    shell:
        """
        curl -L https://osf.io/cxp2w/download -o {output}
        """


