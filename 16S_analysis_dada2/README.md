# Preparing and analysing 16S data using dada2 and the R phyloseq-package 

dada2 is an R package, which uses machine learning to determine sequencing errors in amplicon reads and correct them, thereby reducing falsely increased species richness. As each sequencing run and each sample type has its own characteristics, dada2 should be run separately on those data subsets. 

dada2 as well provides full quality trimming and annotation support (except primer removal) and has a direct connection to the phyloseq package. Very helpful is the visualization of amplicon read quality (as known by fastqc), which can be used independently of the rest of the workflow.

The complete workflow was adapted from this source:  https://f1000research.com/articles/5-1492/v2

## install dada2 and setup directories

A possible way to set up directories, download required data bases and link your data to the work directory is described in `00_setup_conda_R_for_dada2.sh`. Please note that the bioconductor installation of phyloseq, dada2 and dependencies takes quite long in R as well.

## cutadapt

To use dada2, we need to remove primersequences from the amplicon read libraries. We will use `cutadapt` to do so in script `01_cutadapt_tested_static.sh`. I also tried to develop a version with argument parsing, but I haven't tested this yet: `01_cutadapt_untested.sh`

## dada2 (V 1.8 in conda environment)

With adapterfree reads we can start with the dada2 workflow. If you want to learn or test dada2, you can follow the `xx_tutorial_dada2.r` script.

I splitted my data set into the 3 different sequence runs and ran an adjusted script on each of those. For each data set, many trimming parameters were tested. To work in parallel, I performed each data subset analysis on another server.

`02_dada2_glyph_water_dna.r`

`02_dada2_glyph_water_cdna.r`

`02_dada2_glyph_biofilm_dna_cdna.r`

Work in progress:

* the workflow suggests generating a tree before merging the phyloseq object. I would recommend to first get the data (OTU table, taxonomy table, meta table, OTU seqs) into phyloseq, remove singletons and then start calculating to reduce computational effort.

* I'm trying to combine the separated outputs from the sequencing data. This worked fine (`03_dada2_glyph_water_mergetest.r`), but it might be possible in an more efficient way, where I don't save and load whole workspaces, but only the specific sequence table?

## phyloseq

Once the separated data is combined and with its representative sequences completely in phyloseq, we can interchange scripts for mothur and dada2 processed data. At the moment, mothur analysis is more progressed in phyloseq:

The following scripts are based on the phyloseq-objects or the phyloseq-melt-objects to target specific questions:

* get the sequence data into phyloseq `missing`

* generate relative abundance, remove singletons, melt phyloseq object for plotting: copy from `16S_analysis_mothur/02_import_mothur_to_phyloseq.r`

* calculate phylogenetic tree within phyloseq: copy from `16S_analysis_mothur/02_import_mothur_to_phyloseq.r`

* evaluate the distribution of OTUs per genus, their abundance and diversity: copy from `16S_analysis_mothur/02_overview_on_abundant_OTUs_and_genera.r`

* bar-plot the microbial composition on order level (water column, DNA and RNA) over all time points `16S_analysis_mothur/03_community_overview_plot.r`

* plot the abundance of a specific OTU over time for water column and biofilm, treatment and control, DNA and RNA `16S_analysis_mothur/04_*_OTU_abundance_plots.r`

* ordination plots `missing`

* DESeq2-Test for differentially abundant OTUs before and after glyphosate addition `missing`

* Heatmap showing phylogenetic relation and reaction to glyphosate `missing`






