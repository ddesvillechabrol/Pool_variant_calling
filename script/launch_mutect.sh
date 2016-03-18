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

## Usage: launch_mutect.sh -a FILE -r FILE -i FILE -o DIR
##  
##          -h, --help          Show help options.
##          -l, --licence       Print licence info.
##          -a, --assembly      File with contigs.
##          -r, --reference     Reference's bam file.
##          -i, --input         Pool's bam file.
##          -o, --output        Output directory.
## 
## You must use this versions of software:
##          -muTect v1.1.4
##

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

HELP=$(grep "^## " "${BASH_SOURCE[0]}" | cut -c 4-)
LICENCE=$(grep "^#- " "${BASH_SOURCE[0]}" | cut -c 4-)

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
    -a|--assembly)
    ASSEMBLY="$2"
    shift
    ;;
    -r|--reference)
    REFERENCE="$2"
    shift
    ;;
    -o|--output)
    OUTPUT="$2"
    shift
    ;;
    -i|--input)
    INPUT="$2"
    shift
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done


if [[ -z ${ASSEMBLY} ]]
then
    echo "${HELP}"
    echo -e '\n /!\ ERROR /!\ \n'
    echo "You must put a genome reference file (fasta)."
    exit 1
fi

if [[ -z ${REFERENCE} ]]
then
    echo "${HELP}"
    echo -e '\n /!\ ERROR /!\ \n'
    echo "You must put a reference BAM file."
    exit 1
fi

if [[ -z ${INPUT} ]]
then
    echo "${HELP}"
    echo -e '\n /!\ ERROR /!\ \n'
    echo "You must put a input BAM file."
    exit 1
fi

if [[ -z ${OUTPUT} ]]
then
    echo "${HELP}"
    echo -e '\n /!\ ERROR /!\ \n'
    echo "You must put an output file.(.stats.txt)"
    exit 1
fi

make_dir ${OUTPUT%/*}

muTect --analysis_type MuTect --reference_sequence ${ASSEMBLY}  \
    --input_file:normal ${REFERENCE} --input_file:tumor ${INPUT} --out ${OUTPUT} 
