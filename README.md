**Abstract**

The widespread herbicide glyphosate has been detected in aquatic coastal zones of the southern Baltic Sea with yet unknown consequences for this brackish ecosystem. We investigated the impact of glyphosate on the succession of bacterial community assemblages and functions in brackish and phosphorus-limited microcosms. These were incubated for 69 days prior adding of glyphosate to obtain stable community dynamics; eventually, the microcosms were incubated for further 72 days. The system was sampled up to twice a week for the whole incubation period and analyzed concerning glyphosate degradation and bacterial succession. Succession of bacterial assemblages and their functions was determined by cell count analyses, next generation 16S rRNA (gene) amplicon as well as shotgun metagenomic sequencing, respectively. As result, glyphosate was degraded without detectable amounts of aminomethylphosphonic acid. Glyphosate addition revealed shifts in the bacterial community compositions, as well as increases in cell counts, the abundance of specific clusters on genus level, microbial diversity, and cluster richness. Metagenomic analyses revealed a shift in the stoichiometric composition of the phn operon, resulting in increased numbers of the phnM gene. The phn operon shift was detected till the end of the experiment whereas all the other described responses disappeared at latest after 29 days when glyphosate concentrations were below 1 mg∙L-1. The phn operon probably encoded for cleaving the C-P bond in glyphosate, yielding sarcosine. Its shift towards clusters harboring more phnM copies was probably an effect of glyphosate utilization as a P resource by a P-limited bacterial community. Thus, glyphosate impacted bacterial assemblages in a brackish system on different levels, but especially their functional gene compositions for a prolonged period of time. Further analyses should deepen the understanding if whether accordant functional fingerprints could be used as an environmental glyphosate indicator.

-----

**Scripts**

***16S amplicon analysis***

The program `mothur` and the R package `dada2` (together with `cutadapt`) were used to process amplicon reads. 

****set up dada2****

based upon https://f1000research.com/articles/5-1492/v2

Reads for dada2 have to be without primers, we use `cutadapt`. dada2 infers sequencing errors, which usually differ between sequencing runs and habitats. Therefore, I put all data from different runs into separate directories. 

The following directories were set up:

```bash
# these dirs include the gzipped raw reads
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_dna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_cdna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/biofilm

# for reads without primers
/data/projects/glyphosate/reads/reads_16S_cutadapt/water_dna
...

# for reads filtered and trimmed by dada2
/data/projects/glyphosate/reads/dada2_processed/water_dna
...

# we also need database for taxonomy in a dada2 appropriate format
mkdir -p /data/db
wget -O /data/db/silva_nr_v132_train_set.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1'
wget -O /data/db/silva_species_assignment_v132.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1'
```
The primers are removed using the (currently) hard-coded version of a cutadapt script `01_cutadapt_tested_static.sh`, which does it for each directory separately. The script has to be called from within the raw reads directory and output and primer sequences have to be adjusted accordingly. See the manual: https://cutadapt.readthedocs.io/en/stable/guide.html

```bash
#!/bin/bash
for xy in *R*.gz
do
  cutadapt -j 20 -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -o /data/projects/glyphosate/reads/reads_16S_cutadapt/water_dna/cut_${xy} $xy
done
```

We are then ready to install dada2 in a conda environment

```
# install Bioconductor and dada2 in conda, 
conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8 cutadapt=1.18
conda activate dada2
R
```



***Metagenomic analysis***

```parse_prokka_genes_to_bed``` contains scripts that were used to extract a bed.file from prokka annotations and then perform an intersect via samtools to get mapping information on specific genes. Originally developed for my co assembled metagenomes in Sweden, a branch contains the adapted version for the IMP single assemblies

```process_coverage_information``` contains a set of scripts to process the per base coverage of specific genes. It also includes scripts to determine gene length and contig length from the prokka annotation file. This information is used to apply length and coverage thresholds, the remaining data is imported to R and the number of genes, bases and contigs can be calculated there.

```read_abundance``` contains richness information and the script to perform a Mann-Whitney test in R on the reads mapping to a certain gene, testing treatment against control.

```samtools_depth``` gives the per base coverage of a certain gene by using samtools
