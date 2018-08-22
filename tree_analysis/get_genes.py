#!/usr/bin/env python

import os
import re
import urllib.parse
import urllib.request
from ete3 import NCBITaxa

NCBI = NCBITaxa()

def get_desired_ranks(taxid):
    if taxid == -1:
        return {"superkingdom": -1, "phylum": -1, "class": -1, "order": -1, "family": -1,
                "genus": -1}
    lineage = NCBI.get_lineage(taxid)
    lineage2ranks = NCBI.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    for taxrank in ["superkingdom", "phylum", "class", "order", "family", "genus"]:
        if taxrank not in ranks2lineage:
            ranks2lineage[taxrank] = -1
    return ranks2lineage

def create_tax_dict(abspath_names_dmp):
    ncbi_tax_dict = {}
    ncbi_tax_dict[-1] = -1
    with open(abspath_names_dmp) as names_dmp_open:
        for line in names_dmp_open:
            curr_line = re.split(r"\t*\|\t*", line.rstrip())
            if curr_line[3] == "scientific name":
                ncbi_tax_dict[int(curr_line[0])] = curr_line[1]

    return ncbi_tax_dict

tax_dict = create_tax_dict("/tmp/genes_degradation_phn/names.dmp")

tax_list = [31988, 69278, 150203, 41275, 75,
            106589, 323413, 1231242, 745410, 274591,
            85, 135575, 665874, 245186, 149698,
            68287, 359407, 33882, 500577, 256616,
            556257, 286, 238783, 445219, 379,
            1060, 13687, 165697, 1434046, 296051,
            526215]
output_folder = os.path.join(os.getcwd(), "sequences")

phosphonate_degrading = ["phnC", "phnD", "phnE", "phnF", "phnG",
                         "phnH", "phnI", "phnJ", "phnK", "phnL",
                         "phnM", "phnN", "phnO", "phnP"]
sarcosine_oxidizing = [] # still to be resolved
phosphorous_sensing_starvation = ["phoB", "phoR"]
pi_starvation_inducible_psi_genes = ["psiF"]

gene_list = phosphonate_degrading
# gene_list = phosphonate_degrading + sarcosine_oxidizing + phosphorous_sensing_starvation + pi_starvation_inducible_psi_genes

for taxid in tax_list:
    for current_gene in gene_list:
        filename = os.path.join(output_folder, tax_dict[taxid] + "_" + str(current_gene) + ".fasta")

        print(tax_dict[taxid] + "_" + current_gene)
        taxon_queries = ['taxonomy:"%s"' % tid for tid in [taxid]]
        taxon_query = " OR ".join(taxon_queries)
        reviewed = False
        rev = " reviewed:%s" % reviewed if reviewed else ''
        gene = " gene_exact:%s" %current_gene

        url = "https://www.uniprot.org/uniprot/"
        query = "%s%s%s" % (taxon_query, rev, gene)
        params = {"query": query, "force": "yes", "format": "fasta"}
        data = urllib.parse.urlencode(params).encode("utf-8")
        msg = urllib.request.urlretrieve(url=url, filename=filename, data=data)[1]
        headers = {j[0]: j[1].strip() for j in [i.split(':', 1)
                                                for i in str(msg).strip().splitlines()]}

        if "Content-Length" in headers and headers["Content-Length"] == "0":
            os.remove(filename)

