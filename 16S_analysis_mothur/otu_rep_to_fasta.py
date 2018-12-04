#!/bin/python

"""
script takes the output from mothurs get.oturep() command and turns it into a
standard fasta file which can be used in R as DNAstrings object and
therefore included in phyloseqs refseq() slot
"""

import argparse
import re

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

def main():	
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input_file", dest = "input_file", type = str,
						default = None, required = True, action = "store",
	                    help = "input file (get.oturep() fasta file from mothur)")
	parser.add_argument("-o", "--output_file", dest = "output_file", type = str,
						default = None, required = True, action = "store",
	                    help = "output file")
	args = parser.parse_args()

	seq_gap_symbol = [".", "-"]

	with open(args.input_file) as open_fasta_line, open(args.output_file, "w") as open_output_file:
		for line in open_fasta_line:
			if line.startswith(">"):
			# line.split selects the part of the header to look for regex
				fasta_header = re.match(r"(.+?)\|", line.split()[1]).group(1)
				open_output_file.write(">" + fasta_header + "\n")
			else:
				seq = clean_fasta_alignment(line, seq_gap_symbol)
				open_output_file.write(seq + "\n")

if __name__ == "__main__":
	# execute only if run as a script
	main()
