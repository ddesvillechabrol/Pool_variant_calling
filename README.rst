Pool variant calling
#####################

This workflow was developed to analyze SNP in bacteria's pools.

Installation
==============

::

    git clone git@github.com:ddesvillechabrol/pool_variant_calling.git

We download all necessary file with this command.

Example
========

First of all, you must modify config file for your data. You have an example called config.yaml.
This pipeline is developed for the bic cluster at Institut Pasteur.
If you want to use it on other computer, you will need to change path in bash wrapper.

::

    /local/gensoft2/exe/snakemake/3.5.4/bin/snakemake -p -s path/to/Snakefile --configfile path/to/config.yaml \
    --cluster "qsub -q pf4 -cwd -V -b y" --jobs 5


