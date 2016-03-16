from os.path import join

# Import the config file -------------------------------------------------------

configfile: "config.yaml"

# Define some global variable --------------------------------------------------

SPATH = config["script_path"]
FASTQ_DIR = config["fastq_dir"]
SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,[^/]+}.fastq.gz'))
REF = config["ref"]
CONTIGS = config["assembly"]
OUTPUT = config["output"]

# Snakemake --------------------------------------------------------------------

rule all:
    input:
        OUTPUT

rule bwa_sample:
    input:
        fastq = FASTQ_DIR + "/{sample}.fastq.gz",
	    contigs = CONTIGS
    output:
        prefix = "mapped_sample/{sample}",
        bam = "mapped_sample/{sample}.bam",
        bai = "mapped_sample/{sample}.bam.bai"
    run:
        shell(SPATH + "script/launch_bwa.sh -r {input.contigs} -1 {input.fastq} -o {output.prefix}")

rule bwa_ref:
    input:
        fastq = REF,
        contigs = CONTIGS
    output:
        prefix = "mapped_ref/" + REF,
        bam = "mapped_ref/" + REF + ".bam",
        bai = "mapped_ref/" + REF + ".bam.bai"
    run:
        shell(SPATH + "script/launch_bwa.sh -r {input.contigs} -1 {input.fastq} -o {output.prefix}") 

rule mutect:
    input:
        bams = "mapped_sample/{sample}.bam",
        bamr = "mapped_ref/" + REF + ".bam",
        contigs = CONTIGS
    output:
        "mutect/{sample}.stats.txt"
    run:
        shell(SPATH + "script/launch_mutect.sh -a {input.contigs} -r {input.bamr} -i {input.bams} -o {output}")

rule merge_results:
    input:
        stats = expand("mutect/{sample}.stats.txt", sample=SAMPLES),
        contigs = CONTIGS
    output:
        OUTPUT
    params:
        " ".join(config["pipeline_output"]["params"])
    run:
        shell(SPATH + "script/pipeline_output.py {params} -r {input.contigs} -i {input.stats} -o {output}")
