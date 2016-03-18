#!/bin/bash

#- ---------------------------------------------------------------------
#-    Copyright (C) 2015 Dimitri Desvillechabrol
#-
#- This program is free software: you can redistribute it and/or modify
#- it under the terms of the GNU General Public License as published by
#- the Free Software Foundation, either version 3 of the License, or
#- (at your option) any later version.
#-
#- This program is distributed in the hope that it will be useful,
#- but WITHOUT ANY WARRANTY; without even the implied warranty of
#- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#- GNU General Public License for more details.
#-
#- You should have received a copy of the GNU General Public License
#- along with this program.  If not, see <http://www.gnu.org/licenses/>.
#- ---------------------------------------------------------------------

## Usage: launch_bwa.sh -r FILE -single FILE -o FILE.bam
## 
## Required option:
##          -h, --help          Show help options.
##          -r, --reference     Fasta file for reference.
##          -o, --output        Output prefix.
##          -1, --file1         Input fastq file.
##
## Optional input option:
##          -2, --file2         Paired input fastq file.
## 
## You must use this versions of software:
##          - bwa          v0.6.2
##          - samtools     v1.2

mapping_pe(){
    reference=$1
    read1=$2
    read2=$3
    output=$4 # bam of the output
    map_dir=$output"_tmp"
    bwa index $reference;
    bwa aln -t 8 $reference $read1 > $map_dir"_1.sai"
    bwa aln -t 8 $reference $read2 > $map_dir"_2.sai"
    bwa sampe -r '@RG\tID:ref\tSM:ref\tPL:ILLUMINA' $reference \
                 $map_dir"_1.sai" $map_dir"_2.sai" $read1 $read2 > $map_dir.sam
    samtools faidx $reference;
    samtools import $reference.fai $map_dir.sam $map_dir.bam;
    samtools sort $map_dir.bam $output;
    samtools index $output.bam;
}

mapping_se(){
    reference=$1
    read1=$2
    output=$3
    map_dir=$output"_tmp"
    bwa index $reference;
    bwa aln -t 8 $reference $read1 > $map_dir"_1.sai"
    bwa samse -r '@RG\tID:ref\tSM:ref\tPL:ILLUMINA' $reference \
                 $map_dir"_1.sai" $read1 > $map_dir.sam
    samtools faidx $reference;
    samtools import $reference.fai $map_dir.sam $map_dir.bam;
    samtools sort $map_dir.bam $output;
    samtools index $output.bam;
}

make_dir() {
    if [[ ! -s $1 ]]
    then
        echo -e "[INFO] $1 directory is create"
        mkdir $1
    fi
}

extract_name(){
    filename=${1##*/}
    name=${filename%%.*}
    echo $name
}

HELP=`grep "^## " "${BASH_SOURCE[0]}" | cut -c 4-`
LICENCE=`grep "^#- " "${BASH_SOURCE[0]}" | cut -c 4-`

while [[ $# > 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo "${HELP}"
    exit 1
    ;;
    -l|--licence)
    echo "${LICENCE}"
    exit 1
    ;;
    -r|--reference)
    REFERENCE="$2"
    shift
    ;;
    -o|--output)
    OUTPUT="$2"
    shift
    ;;
    -1|--file1)
    INPUT1="$2"
    shift
    ;;
    -2|--file2)
    PAIRED=YES
    INPUT2="$2"
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

if [[ -z ${REFERENCE} ]]; then
    echo "${HELP}"
    echo -e '\n /!\ ERROR /!\ \n'
    echo "You must put a reference file (fasta)"
    exit 1
fi


if [[ -z ${INPUT1} ]]; then
    echo "${HELP}"
    echo -e '\n /!\ ERROR /!\ \n'
    echo "You must put reads files (fastq.gz)"
    exit 1
fi

make_dir ${OUTPUT%/*}

if [[ -n $PAIRED ]]; then
    mapping_pe ${REFERENCE} ${INPUT1} ${INPUT2} ${OUTPUT}
else    
    mapping_se ${REFERENCE} ${INPUT1} ${OUTPUT}
fi

