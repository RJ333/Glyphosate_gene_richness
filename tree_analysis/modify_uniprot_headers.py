#!/usr/bin/env python

# I guess this script modifies uniprot headers? :-)

# use python2

import re
import sys

# file = "/tmp/genes_degradation_phn/sequences/Brevundimonas_phnC.fasta"
file = sys.argv[1]

with open(file) as f:
    for line in f:
        line_strip = line.rstrip()
        rx_object = re.search(r"^(>[a-z]+[A-Z0-9\|_]+) .+( OS=.+) OX=", line_strip)
        if rx_object:
            new_header = rx_object.group(1) + rx_object.group(2)
            print(new_header.replace(" ", "@"))
        else:
            print(line_strip)
