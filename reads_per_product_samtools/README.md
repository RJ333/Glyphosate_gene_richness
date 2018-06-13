# README.md

These first script copies the prokka annotation, modifies the product annotations within and uses the modified products to generate a list of unique products across all samples with a given threshold
"Modified" means that the product descriptions are adjusted to decrease redundancy (e.g. for protein_A, Protein_A, Protein-A and protein a) and substitute problematic characters ("/", "_", ...). 


## How to run the script

```
conda activate Renv
time Rscript 02_modify_prokka_products.r -p="/data/Rene/glyph/prokka/prokka_all_modified.tsv" -t=1 -o="/data/Rene/glyph/prokka/"
```

Based on the unique product list a bed file
 is created to perform an intersect between the mapped reads and the products (with samtools). As this is a very computational extensive task, GNU parallel was implemented for full use of cores.
 If anything goes wrong regarding the name of the results files (e.g. you forgot to substitute / or _) the `rename.sh` script helps you. 
 Important: this script moves the files to ./renamed. If you don't use it, set the input folder correct in the later scripts 

When all the files are present and optionally renamed, the last script combines all the data to a single table, containing genes per contigs per sample with the reads per genes. This is the input for R

## How to run the script
start a tmux session as the script will take many days and you wanna logoff securely

All scripts have a $BASE_DIR variable, which needs to be set within the script

adjust the settings within `get_reads_per_prokka_product_GNU_parallel.sh`, namely
```
parallel -j $GNU_CORES sam_parallel ::: A1 A2 A3 A4 A5 A6 A7 B8 B9 B10 :::: $OUTPUT_DIR/unified_unique_prokka_products_greater_${THRESHOLD}.tsv

# how often needs a product be present per sample
THRESHOLD=1
# how many cores should GNU Parallel use for samtools?
GNU_CORES=44

# directories of original data
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP/
ORIGINAL_PROKKA_DIR=/data/jwerner/glyphosate/IMP/${sample}/output_IMP/Analysis/annotation/ 
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools

# working directory
OUTPUT_DIR=/data/Rene/glyph/prokka  

# results directories
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp
```

## Questions/bugs?

Please create an issue.