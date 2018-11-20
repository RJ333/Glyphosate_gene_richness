# Idea 

(already described in https://gitlab.io-warnemuende.de/Workflows/calculate-functional-trees/blob/master/README.md )

In this workflow we perform multipe sequence alignments on protein fasta sequence. The data stems from metagenomic samples annotated with prokka. We also include reference genes.

## Setup directories, get prokka data

This script gets the annotations from all samples and transforms the prokka annotated sequences to a single line `01_general_setup_for_tree_analysis.sh`

## Build a tree with a specific gene

To build a tree you can use `02_build_tree_for_specific_genes.sh`, you need to know the name of the gene you want to perform MSA on

## python scripts by Johannes

`get_genes.py`: I assume this is the step to fetch a list of reference gens via the Uniprot API

`modify_uniprot_headers.py`: This is probably to perform some changes and substitutions within the Uniprot fasta headers

