Pool variant calling
#####################

This workflow was developed to analyze SNP in bacteria's pools.

Installation
==============

::

    git clone git@github.com:ddesvillechabrol/pool_variant_calling.git

We download all necessary file with this command.

Example to use it on SGE cluster
==============================================

After you download the directory from git, you copy this where you want on the cluster.

::

    scp -r pool_variant_calling ddesvill@bic.pasteur.fr:~
    chmod 750 -R pool_variant_calling

We assume you have a working directory with a reference genome (fasta), reads of
a pure sample (fastq) and reads of your pooled sample (fastq) in a directory.
You need to copy the config.yaml inside your working directory and change the 
script_path to indicate where is the "pool_variant_calling" directory.
Actually, the script is setting for the bic cluster of Institut Pasteur. 
You must change the path of module source (module_src).

Inside the config.yaml, you need to change different parameter (ref, fastq_dir,
genome, output).

::

    module load snakemake
    snakemake -p -s path/to/Snakefile --configfile config.yaml --cluster "qsub -q pf4 -cwd -V -b y" --jobs 5

The command runs the workflow on a SGE cluster. But on the bic cluster, you must use this command line:

::

    module load snakemake
    /local/gensoft2/exe/snakemake -p -s path/to/Snakefile --configfile config.yaml --cluster "qsub -q pf4 -cwd -V -b y" --jobs 5
