# README.md

These two scripts `parse_prokka_output.sh` and `parse_annotation.py` parse the split prokka output and prepare a bed file to perform an intersect between the mapped reads and the genes of interest (with samtools). 

## How to run the script

```
./parse_prokka_output.sh -i <INPUT_DIR> -o <OUTPUT_DIR> [-g <GENE_NAME>] [-b <BED_FILE>]
```

## How to use the bed file afterwards

Use samtools, e.g. like this

```
cd <MAPPING_DIR>
for i in *.bam
do
    echo $i
    samtools view -@ 26 -L <BED_FILE> $i | cut -f 3 | sort | uniq | wc -l
done
```

## Questions/bugs?

Please create an issue.