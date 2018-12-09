****set up dada2****

based upon https://f1000research.com/articles/5-1492/v2. If you want to follow this tutorial, install dada2 as below and look into `tutorial_dada2.r`

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
There is also a version including argument passing, which is not tested yet (`01_cutadapt_untested.sh`).

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

I performed parallel runs of dada2, as each sequencing run data needs interactive and individual caretaking: the plotting of the sequence quality will show you where to trim your reads. The outcome of the scripts `02_dada2...r` are `.RData`-files, which contain a sequence table (Seqtab). At the moment, all workspaces need to be loaded to merge the Seqtabs from the parallel runs into one object (`03_dada2_glyph_water_mergetest.r`). Those sequences will then be assigned taxonomically and turned into an `phyloseq`-object.

****phyloseq analysis****

about to come, this is supposed to be the same for `dada2` and `mothur`-processed data