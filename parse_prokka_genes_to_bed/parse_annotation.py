#!/usr/bin/env python

"""
Documentation for script
"""

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_file", dest="input_file", type=str, default=None, required=True, action="store",
                    help="input file (filtered gff file)")
parser.add_argument("-c", "--contig_length", dest="contig_length_file", type=str, default=None, required=True,
                    action="store", help="contig length file")
parser.add_argument("-g", "--gene_name", dest="gene_name", type=str, default=None, required=True, action="store",
                    help="gene to search for")
parser.add_argument("-o", "--output_file", dest="output_file", type=str,default=None, required=True, action="store",
                    help="output file (bed file)")
args = parser.parse_args()

phntmp_filename = open(args.input_file, "r")
contiglength_filename = open(args.contig_length_file)
output_filename = open(args.output_file, "w")

gene = args.gene_name
search_string = "gene=" + gene
d = {}

def write_output_to_file(contig_name, start, end, gene_name):
    return "\t".join(str(x) for x in (contig_name, start, end, gene_name)) + "\n"

# create dictionary with contig names and contig lengths
for entry in contiglength_filename:
    if entry.startswith("##sequence-region"):
        entry_split = entry.split()
        if entry_split[1] not in d:
            d[entry_split[1]] = int(entry_split[-1])

for line in phntmp_filename:
    line_split = line.rstrip().split("\t")
    # add entry to output file if search_string (=gene) is found
    if search_string in line_split[3]:
        contig_name = line_split[0]
        contig_length = d[contig_name]
        start = line_split[1]
        end = line_split[2]
        gene_name = re.findall(search_string, line_split[3])[0]
        gene_name = gene_name.split("gene=")[1]


        if contig_length < 20000:
            output_filename.write(write_output_to_file(contig_name, start, end, gene_name))
        else:
            index = 0
            while contig_length > 20000:
                contig_length = contig_length - 10000
                contig_name_new = contig_name + "." + str(index)
                if index != 0:
                    start = str(int(start) - 10000)
                    end = str(int(end) - 10000)
                if int(start) > 0 and int(start) < 10000:
                    output_filename.write(write_output_to_file(contig_name_new, start, end, gene_name))
                index += 1
            contig_name_new = contig_name + "." + str(index)
            start = str(int(start) - 10000)
            end = str(int(end) - 10000)
            if int(start) > 0 and int(start) < 10000:
                output_filename.write(write_output_to_file(contig_name_new, start, end, gene_name))
