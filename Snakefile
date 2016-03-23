from os.path import join

# Import the config file -------------------------------------------------------

configfile: "config.yaml"

# Define some global variable --------------------------------------------------
get_prefixes = lambda filename: filename.split(".")[0]
REF = get_prefixes(config["ref"])
print(REF)

SPATH = config["script_path"]
FASTQ_DIR = config["fastq_dir"]
SAMPLES, = glob_wildcards(join(FASTQ_DIR, '{sample,[^/]+}.fastq.gz'))
CONTIGS = config["genome"]
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
        bam = "mapped_sample/{sample}.bam",
        bai = "mapped_sample/{sample}.bam.bai"
    params:
        prefix = "mapped_sample/{sample}",
        spath = SPATH,
        src = config["module_src"],
        module = config["mapping_module"]
    shell:
        """
        . {params.src}
        module load {params.module}
        {params.spath}/script/launch_bwa.sh -r {input.contigs} \
        -1 {input.fastq} -o {params.prefix}
        """

rule bwa_ref:
    input:
        fastq = config["ref"],
        contigs = CONTIGS
    output:
        bam = "mapped_ref/" + REF + ".bam",
        bai = "mapped_ref/" + REF + ".bam.bai"
    params:
        prefix = "mapped_ref/" + REF,
        spath = SPATH,
        src = config["module_src"],
        module = config["mapping_module"]
    shell:
        """
        . {params.src}
        module load {params.module}
        {params.spath}/script/launch_bwa.sh -r {input.contigs} \
        -1 {input.fastq} -o {params.prefix}
        """ 

rule mutect:
    input:
        bams = "mapped_sample/{sample}.bam",
        bamr = "mapped_ref/" + REF + ".bam",
        contigs = CONTIGS
    output:
        "mutect/{sample}.stats.txt"
    params:
        spath = SPATH,
        src = config["module_src"],
        module = config["mutect_module"]
    shell:
        """
        . {params.src}
        module load {params.module}
        {params.spath}/script/launch_mutect.sh -a {input.contigs} \
        -r {input.bamr} -i {input.bams} -o {output}
        """

rule merge_results:
    input:
        stats = expand("mutect/{sample}.stats.txt", sample=SAMPLES),
        contigs = CONTIGS
    output:
        OUTPUT
    params:
        spath = SPATH,
        filin = " ".join(config["pipeline_output"]["params"])
    shell:
        """
        {params.spath}/script/pipeline_output.py {params.filin} \
        -r {input.contigs} -i {input.stats} -o {output}
        """
