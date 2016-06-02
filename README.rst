Pool variant calling
#####################

This workflow was developed to analyze SNP in bacteria's pools.

Installation
==============

::

    git clone https://github.com/ddesvillechabrol/pool_variant_calling.git

We download all necessary file with this command. 
The workflow need bwa/0.6.2, samtools/1.2 and muTect/1.1.4 available as module on the cluster.

Example to use it on SGE cluster
==============================================

We will assume you have a working directory with a reference genome (fasta),
reads of a pure sample (fastq) and reads of your pooled sample (fastq) in a
directory.
You need to copy the config.yaml inside your working directory and change the 
script_path to indicate where is the "pool_variant_calling" directory.
Actually, the script is setting for the bic cluster of Institut Pasteur. 
You must change the path of module source (module_src) to use it on classical
SGE cluster.

Inside the config.yaml, you need to change different parameter (ref, fastq_dir,
genome, output).

::

    module load snakemake
    snakemake -p -s path/to/Snakefile --configfile config.yaml --cluster "qsub -q pf4 -cwd -V -b y" --jobs 5

The command runs the workflow on a SGE cluster. But on the bic cluster, you must use this command line:

::

    module load snakemake
    /local/gensoft2/exe/snakemake -p -s path/to/Snakefile --configfile config.yaml --cluster "qsub -q pf4 -cwd -V -b y" --jobs 5
