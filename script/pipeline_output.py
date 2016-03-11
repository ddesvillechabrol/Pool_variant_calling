#!/usr/bin/env python2.7
# _*_ coding:utf-8 _*_
#
# ---------------------------------------------------------------------
#    Copyright (C) 2015 Dimitri Desvillechabrol
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ---------------------------------------------------------------------
#


""" Tool to generate usefull ouput file with MuTect output.

This script take MuTect stats files and assembly file to generate
usefull tabulate file. If you use annotation file (.gbf or .gbk) and a
complete genome, script add usefull information.

"""

import sys
import re
import argparse
import string

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


VARIANT = set()


class Contig:
    """Contains annotation information for a contig.

    It takes a SeqRecord object to generate a Contig object which
    have all information from prokka or rast (.gbf or .gbk).

    """
    def __init__(self, seq_record):
        self.contig = seq_record.name  # Name of the contig
        # Retrieve all genomic informations inside SeqRecord object
        self.region = []
        for feature in seq_record.features:
            if feature.type.find("source") < 0:  # Check if feature is present
                new_region = Region(feature, seq_record.seq)
                if new_region.remove_gene():
                    self.region.append(new_region)


class Region:
    """Contains annotation of one region of a contig.

    It takes a feature from SeqRecord object and sequence to
    generate a region object.

    """
    def __init__(self, feature, seq):
        self.start = feature.location.start  # Start of region in contigs
        self.end = feature.location.end  # End of region in contigs
        self.size = self.end - self.start  # Size of region
        self.strand = feature.location.strand  # Strand of region
        self.gene = feature.type  # Type of region (CDS or other)
        # Check if the region is a repeat region.
        try:
            self.product = feature.qualifiers["product"][0]  # Catch annotation
        except KeyError:
            self.product = "n/a"
        self.seq = seq[self.start:self.end]  # Retrieve sequence in contigs

    def remove_gene(self):
        if self.gene.startswith("gene"):
            return False
        return True

    def do_mutation(self, pos_base, alt_base):
        """ Return the sequence with a mutation.


        It mutates the sequence (self.seq). It takes the position
        and the alternate base as input.

        """
        mutable_seq = self.seq.tomutable()  # Sequence can be mutated
        mutable_seq[pos_base] = alt_base  # Mutates base
        return mutable_seq.toseq()

    def get_aa(self, pos, alt_base):
        """Return amino acid position, reference amino acid and
        alternate amino acid.

        It takes the position and the alternate base as input.

        """
        # Check if the genes is on the current strand or on the other.
        if self.strand > 0:
            pos_nuc = pos - 1 - self.start  # Set position in region
            pos_aa = int(pos_nuc/3)  # Position of the amino acid
            ref_aa = self.seq.translate()[pos_aa]  # Catch reference amino acid
            alt_seq = self.do_mutation(pos_nuc, alt_base)  # Do mutation
            alt_aa = alt_seq.translate()[pos_aa]  # Catch alternate amino acid
            return [str(pos_aa), ref_aa, alt_aa]
        else:
            pos_nuc = pos - 1 - self.start  # Set position in region
            pos_aa = int((self.end - pos - 1) / 3)  # Position of the AA
            ref_aa = self.seq.reverse_complement().translate()[pos_aa]
            alt_seq = self.do_mutation(pos_nuc, alt_base)  # Do mutation
            alt_aa = alt_seq.reverse_complement().translate()[pos_aa]
            return [str(pos_aa), ref_aa, alt_aa]


class Mutect:
    """Contains variants detected by MuTect.

    It takes a MuTect file to generate a MuTect object.

    """
    def __init__(self, fl):
        # Take the name of file without path and extension
        resultat = re.findall("/?([^/]+)", fl)[-1]
        resultat = re.search("([^\.]+)\.", resultat)
        try:
            self.name = resultat.group(1)
        except ValueError:
            self.name = fl
        # Take variant selected and rejected by MuTect.
        self.variant_keep, self.variant_reject = self.extract_variant(fl)

    def extract_variant(self, fl):
        """ Return dictionaries of variant keep and of variant reject.

        It takes a MuTect file as input.

        """
        keep_list = {}
        reject_list = {}
        # Regex to retrieve conserved variants by MuTect
        regex_keep = re.compile("KEEP")
        with open(fl, "r") as fp:
            for line in fp:
                # Do not read header of files
                if line.startswith("#") or line.startswith("contig"):
                    continue
                new_variant = Variant(line)
                if new_variant.end_filter():
                    continue
                if regex_keep.search(line):
                    keep_list["{}_{}".format(new_variant.contig,
                                             new_variant.pos)] = new_variant
                    VARIANT.add("{}_{}".format(new_variant.contig,
                                               new_variant.pos))
                else:
                    reject_list["{}_{}".format(new_variant.contig,
                                               new_variant.pos)] = new_variant
        return keep_list, reject_list

    def get_name(self):
        """Return name of MuTect file.

        """
        return self.name

    def get_variant(self, variant, recheck):
        """Return variant's name.

        Check in rejected variant if a variant is present with a
        frequency higher than 2%.

        """
        if variant in self.variant_keep:
            return self.variant_keep[variant]
        if recheck:
            if variant in self.variant_reject:
                if self.variant_reject[variant].get_freq() >= 0.02:
                    return self.variant_reject[variant]
        return None


class Variant:
    """Contains informations of one variant.

    It takes a line from MuTect file and keep information.
    """
    def __init__(self, line):
        row = line.split()
        contig_info = row[0].split("_")
        self.contig = int(contig_info[1])  # Contig number
        self.size = int(contig_info[3])  # Contig size
        self.pos = int(row[1])  # Position in contig
        self.ref_base = row[3]  # Reference base
        self.alt_base = row[4]  # Alternate base
        self.freq = float(row[17])  # Alternate base frequency

    def get_freq(self):
        """Return frequency of alternate allele.

        """
        return self.freq

    def get_contig(self):
        """Return number of contig

        """
        return self.contig

    def end_filter(self):
        """Return boolean.

        Check if a variant is from the end of a contig. False for no
        and True for yes.

        """
        if self.pos < 100 or self.size - self.pos < 100:
            return True
        return False


def arguments():
    """Return Namespace object.

    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", help="Complete genome (fasta)",
                        default="check_string_for_empty")
    parser.add_argument("-a", "--annotation",
                        help="Prokka genbank file (.gbf)",
                        default="check_string_for_empty")
    parser.add_argument("--recheck", help="Recheck the MuTect file", 
                        action="store_true")
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument("-r", "--reference",
                               help="Reference file (fasta)", required=True)
    requiredNamed.add_argument("-i", "--input", nargs="+", help="MuTect file",
                               required=True)
    requiredNamed.add_argument("-o", "--output", help="Output file name",
                               required=True)
    return parser.parse_args()


def multi_gbk_parse(gbk_file):
    """Return dictionarie of annotation.

    It takes a genbank file (.gbk or .gbf) as input.

    """
    annotation = {}
    for seq_record in SeqIO.parse(gbk_file, "genbank"):
        name = seq_record.name.split("_")[1]  # Takes contig number
        annotation[name] = Contig(seq_record)  # Takes all contigs information
    return annotation


def assembly_parse(fasta_file):
    """Return dictionarie with different contigs.

    """
    assembly = {}
    for contigs in SeqIO.parse(fasta_file, "fasta"):
        name = contigs.id.split("_")[1]  # Takes contig number
        assembly[name] = contigs.seq  # Takes contigs sequence
    return assembly


def read_mutect_files(mutect_files):
    """Return list of MuTect object.

    """
    mutect_list = []
    for mutect in mutect_files:
        mutect_list.append(Mutect(mutect))
    return mutect_list


def check_optional_argument(argument):
    """Return boolean.

    Check if an argument is present.

    """
    if argument == "check_string_for_empty":  # Check the default value
        return False
    return True


def get_header(annotation, ref_genome, mutect):
    """Return string of header.

    """
    header = "contig\tposition\tref_res\talt_res\tsequence"
    if ref_genome:
        header = "{}\t{}".format(header, "sense")
    for pool in mutect:
        header = "{}\t{}".format(header, pool.get_name())
    if annotation:
        annot = "type\tsize\taa_pos\taa_ref\taa_mut\tannotation"
        header = "{}\t{}".format(header, annot)
    header = "{}\n".format(header)
    return header


def seek_motif(seq, genome_ref, alt):
    """Return if the sequence is wild type or not.

    """
    table = string.maketrans('CGAT', 'GCTA')
    seq_alt = "{}{}{}".format(seq[0:10], alt, seq[11:])  # Do mutation
    # Set the second strand of the short sequence
    seq_comp = seq.translate(table)
    seq_alt_comp = seq.translate(table)
    # Check all possibility for the short sequence against the sequence
    for genome in genome_ref:
        ref = genome.seq
        if ref.find(seq) >= 0:
            return "ref"
        if ref.find(seq[::-1]) >= 0:
            return "ref"
        if ref.find(seq_comp) >= 0:
            return "ref"
        if ref.find(seq_comp[::-1]) >= 0:
            return "ref"
        if ref.find(seq_alt) >= 0:
            return "alt"
        if ref.find(seq_alt[::-1]) >= 0:
            return "alt"
        if ref.find(seq_alt_comp) >= 0:
            return "alt"
        if ref.find(seq_alt_comp[::-1]) >= 0:
            return "alt"
    return "abs"


def get_variant_info(ref_genome, assembly, mutect, variant, recheck):
    """Return string to write in the output.

    It takes reference genome, assembly, list of MuTect and a variant
    and makes a line with position of variant, what variant is and
    frequencies inside each pool.

    """
    freq_list = []
    for pool in mutect:
        pool_variant = pool.get_variant(variant, recheck)
        # Retrieve all frequencies in each pool
        if pool_variant:
            freq_list.append("{0:.2f}".format(pool_variant.get_freq()))
            info = pool_variant
        else:
            freq_list.append("0.00")
    pos = int(info.pos) - 1
    seq = assembly[variant.split("_")[0]][pos-10:pos+11]
    seq_lower = "".join(c.lower() if i != 10 else c  for i, c in enumerate(seq))
    # Line with information from MuTect
    row = "{}\t{}\t{}\t{}\t{}".format(info.contig, info.pos,
                                      info.ref_base, info.alt_base, seq_lower)
    if ref_genome:  # Check optional option
        # Return the sense of the mutation
        sense = seek_motif(seq.tostring(), ref_genome, info.alt_base)
        return "{}\t{}\t{}".format(row, sense, "\t".join(freq_list))
    return "{}\t{}".format(row, "\t".join(freq_list))


def get_annotation(annotation, variant_line):
    """Return annotation field for the output.

    """
    variant_info = variant_line.split()  # Takes information of variant
    try:
        contig_info = annotation[variant_info[0]]
    except KeyError:
        print(" /!\\ ERROR /!\\ ")
        print("Check if you use the good annotation file.")
        sys.exit()
    # Check each region inside the contig which contain the variant
    for region_info in contig_info.region:
        if region_info.end > int(variant_info[1]):
            if region_info.start < int(variant_info[1]):
                if region_info.gene.startswith("CDS"):
                    # Retrieve amino acid information
                    aa_info = region_info.get_aa(int(variant_info[1]),
                                                 variant_info[3])
                    return "{}\t{}\t{}\t{}".format(region_info.gene,
                                                   region_info.size/3,
                                                   "\t".join(aa_info),
                                                   region_info.product)
                # If is not in CDS region
                return "{}\tn/a\tn/a\tn/a\tn/a\t{}".format(region_info.gene,
                                                           region_info.product)
    return "Intergenic\tn/a\tn/a\tn/a\tn/a\tn/a"  # If there are no annotation


def read_genome(fl_genome):
    """Return list of seq_record.

    """
    genome = []
    for seq_record in SeqIO.parse(fl_genome, "fasta"):
        genome.append(seq_record)
    return genome


def write_file(args, assembly, mutect):
    """Write the output file.

    """
    ref_genome = None  # Complete genome of reference
    annotation = None  # Annotation of the reference
    variant_annotation = ""
    # Check if we have access of optional arguments
    if check_optional_argument(args.annotation):
        annotation = multi_gbk_parse(args.annotation)  # Read gbk file
    if check_optional_argument(args.genome):
        ref_genome = read_genome(args.genome)  # Read fasta file
    with open(args.output, "w") as fp:
        fp.write(get_header(annotation, ref_genome, mutect))  # Write header
        for variant in VARIANT:
            line = get_variant_info(ref_genome, assembly, mutect, variant,
                                    args.recheck)
            if annotation:
                variant_annotation = get_annotation(annotation, line)
            fp.write("{}\t{}\n".format(line, variant_annotation))


if __name__ == "__main__":
    # Parse arguments with argparse
    args = arguments()
    # Read assembly file
    assembly = assembly_parse(args.reference)
    # Read MuTect files
    mutect = read_mutect_files(args.input)
    # Write file
    write_file(args, assembly, mutect)
