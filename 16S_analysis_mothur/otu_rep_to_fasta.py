#!/bin/python

"""
script takes the output from mothurs get.oturep() command and turns it into a
standard fasta file which can be used in R as DNAstrings object and
therefore included in phyloseqs refseq() slot
"""

import argparse
import re

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", dest = "input_file", type = str,
					default = None, required = True, action = "store",
                    help = "input file (get.oturep() fasta file from mothur)")
parser.add_argument("-o", "--output_file", dest = "output_file", type = str,
					default = None, required = True, action = "store",
                    help="output file")
args = parser.parse_args()

input_filename = open(args.input_file, "r")
output_filename = open(args.output_file, "w")
seq_gap_symbol = [".", "-"]

def clean_fasta_alignment(line, symbols):
	"""
	the for loop deletes the gap symbols in each fasta line.

	Args:
        line (str)
        symbols (str): The characters to be removed from the fasta sequence

	Return:
		sequence lines without gap chars
	"""
    for symbol in symbols:
        line = line.rstrip().replace(symbol, "")
    return line

# this is the actual formatting process
with input_filename as fasta:
    for line in fasta:
        if line.startswith(">"):
            fasta_header = re.match(r"(.+?)\|", line.split()[1]).group(1)
            output_filename.write(">" + fasta_header + "\n")
        else:
            seq = clean_fasta_alignment(line, seq_gap_symbol)
            output_filename.write(seq + "\n")
