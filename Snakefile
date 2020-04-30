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

rule all:
    input:
        # create a new filename for every entry in SAMPLES,
        # replacing {name} with each entry.
        expand("rnaseq/quant/{name}_quant/quant.sf", name=SAMPLES),
        "rnaseq/fastqc/multiqc_report.html",

# download human rna-seq data from Duclose et al, 2019 study
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

### download and index the human transcriptome ###
rule download_human_transcriptome:
    output: "rnaseq/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz" 
    shell:
        """
        curl -L ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -o {output}
        """

rule salmon_index:
    input: "rnaseq/reference/Homo_sapiens.GRCh38.cdna.all.fa.gz" 
    output: directory("rnaseq/quant/sc_ensembl_index")
    conda: "rnaseq-env.yml"
    shell:
        """
        salmon index --index {output} --transcripts {input} # --type quasi
        """

### quantify reads with salmon
rule salmon_quantify:
    input:
        reads="rnaseq/quality/{sample}.qc.fq.gz",
        index_dir="rnaseq/quant/sc_ensembl_index"
    output: "rnaseq/quant/{sample}_quant/quant.sf"
    params:
        outdir= lambda wildcards: "rnaseq/quant/" + wildcards.sample + "_quant"
    conda: "rnaseq-env.yml"
    shell:
        """
        salmon quant -i {input.index_dir} --libType A -r {input.reads} -o {params.outdir} --seqBias --useVBOpt --validateMappings
        """

