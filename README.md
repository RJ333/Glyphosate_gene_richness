##Abstract

The widespread herbicide glyphosate has been detected in aquatic coastal zones of the southern Baltic Sea with yet unknown consequences for this brackish ecosystem. We investigated the impact of glyphosate on the succession of bacterial community assemblages and functions in brackish and phosphorus-limited microcosms. These were incubated for 69 days prior adding of glyphosate to obtain stable community dynamics; eventually, the microcosms were incubated for further 72 days. The system was sampled up to twice a week for the whole incubation period and analyzed concerning glyphosate degradation and bacterial succession. Succession of bacterial assemblages and their functions was determined by cell count analyses, next generation 16S rRNA (gene) amplicon as well as shotgun metagenomic sequencing, respectively. As result, glyphosate was degraded without detectable amounts of aminomethylphosphonic acid. Glyphosate addition revealed shifts in the bacterial community compositions, as well as increases in cell counts, the abundance of specific clusters on genus level, microbial diversity, and cluster richness. Metagenomic analyses revealed a shift in the stoichiometric composition of the phn operon, resulting in increased numbers of the phnM gene. The phn operon shift was detected till the end of the experiment whereas all the other described responses disappeared at latest after 29 days when glyphosate concentrations were below 1 mg∙L-1. The phn operon probably encoded for cleaving the C-P bond in glyphosate, yielding sarcosine. Its shift towards clusters harboring more phnM copies was probably an effect of glyphosate utilization as a P resource by a P-limited bacterial community. Thus, glyphosate impacted bacterial assemblages in a brackish system on different levels, but especially their functional gene compositions for a prolonged period of time. Further analyses should deepen the understanding if whether accordant functional fingerprints could be used as an environmental glyphosate indicator.

-----

**Scripts**

###16S amplicon analysis

The program `mothur` (1.39.5) and the R package `dada2` 1.8 (together with `cutadapt` 1.8.3) were used to process amplicon reads.

****set up mothur****

set up the directories: some of the structure was already generated for the `dada2`-analysis. We will use the same raw reads and add a directory which will become the working directory

```bash
# the dirs with raw reads
/data/projects/glyphosate/reads/raw_reads_16S/*
# create working dir
mkdir -p /data/projects/glyphosate/reads/mothur_processed

# add databases
wget -O /data/db/Silva.seed_v132.tgz 'https://mothur.org/w/images/7/71/Silva.seed_v132.tgz'
wget -O /data/db/Silva.nr_v132.tgz 'https://mothur.org/w/images/3/32/Silva.nr_v132.tgz'
cd /data/db 
tar -xvzf /data/db/Silva.seed_v132.tgz
tar -xvzf /data/db/Silva.nr_v132.tgz
# link databases to working dir
ln -s /data/db/silva.seed_v132.align /data/projects/glyphosate/reads/mothur_processed

# create soft links to reads from different dirs to mothur input dir
ln -s /data/projects/glyphosate/reads/raw_reads_16S/*/*.gz /data/projects/glyphosate/reads/raw_reads_16S/

# you can remove the controls
cd /data/projects/glyphosate/reads/raw_reads_16S/
rm *pos*.gz
rm *neg*.gz

# install mothur into conda environment (versions from 1.4* are slow)
conda env create -f conda_mothur.conf 

# you will also need a file called `oligo.txt` in your working directory which specifies the primers. Below is an example for the primer pair 341f - 805r
cat oligo.txt
forward CCTACGGGNGGCWGCAG
reverse GACTACHVGGGTATCTAATCC

# additionally you will have to find out where your primers bind. This allows you to generate a database for only this region, which reduces the computational effort

# run mothur with
conda activate mothur_1395
mothur
```
****mothur workflow****

The mothur workflow is described in  `01_mothur_workflow_1395.h`. The goal of this script is to perform all steps to generate an OTU table, a taxonomic annotation and the amplicon sequences representing each OTUs. Relative abundance, singleton removal and tree building will later be performed within `phyloseq`.

I suggest to separately store the output from the `summary.seqs()`-command, which gives a good overview how your data set changes (`summary_seqs_1395_collected.txt`). The output is also contained in the log file. 

****adjust stability.files****

The first step in the workflow is to generate a file listing all the sample read pairs. The group names in this `stability.files` are wrong (they contain the full path) and it is created in the wrong directory, you have to adjust them using awk. "group" means sample or library here. . Use this two-liner, which renames `stability.files` and moves it to the input folder (where the reads are):
The workflow contains the outcommented version of the required code, it has to be executed outside of mothur

****characteristics of mothur (1.39.5)****

mothur is written in C++. You need to set the input and output in the same command. mothur clears the values that are not set when calling `set.dir()`. mothur also looks in different directories, when it can't find a file in the input dir, but there is a pasting bug with "/" so it won't work until it is in the inputdir.

mothur does not understand line breaks! There is also a command line mode and a batch mode, but I have not used those. 

****memory demands****

The command `cluster.split()` is likely to demand a huge amount of RAM, also based on the number of cores you used before. It crashed for me several times, even on `phy-2` with 1 TB of RAM. The workflow includes commands to transfer the relevant data between e.g. the cloud and phy-2.

```bash
# copy the relevant files from your work dir to phy-2 (via bio-48)
ssh bio-48
cd /.../<mothur_work_dir_on_bio48>
scp user_name@IP:/path/to/latest_fasta .
scp user_name@IP:/path/to/latest_count_table .
scp user_name@IP:/path/to/latest_taxonomy .

# run mothur 1.39.5 with screen on phy-2, set input and output dir to working dir
ssh phy-2
screen
/dss6/bio49/bin/mothur/mothur
set.dir(input = /dss6/bio49/.../<mothur_work_dir_on_bio48>, output = /dss6/bio49/.../<mothur_work_dir_on_bio48>)
# specify processors and cutoff value
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 25)

# copy files back, except for dist-file, which is huge and not needed
# on bio-48:
cd /.../<mothur_work_dir_on_bio48>
scp stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list user_name@IP:/path/to/mothur_work_dir
scp stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.sensspec user_name@IP:/path/to/mothur_work_dir
```
****generate OTU sequence fasta file****

Within this workflow, we will generate a file containing aligned sequences of all OTUs, based (here) on the most abundant sequence per OTU. Using the python script `otu_rep_to_fasta.py`, we can turn this output into proper fasta format, which can be integrated into `phyloseq`. Run it as 

```bash
cd <mothur_working_dir>
python /path/to/scripts/otu_rep_to_fasta.py \
  -i <your_file.rep.fasta> \
  -o <your_proper_fasta.fasta>
```

****Metagenomic analysis****

```parse_prokka_genes_to_bed``` contains scripts that were used to extract a bed.file from prokka annotations and then perform an intersect via samtools to get mapping information on specific genes. Originally developed for my co assembled metagenomes in Sweden, a branch contains the adapted version for the IMP single assemblies

```process_coverage_information``` contains a set of scripts to process the per base coverage of specific genes. It also includes scripts to determine gene length and contig length from the prokka annotation file. This information is used to apply length and coverage thresholds, the remaining data is imported to R and the number of genes, bases and contigs can be calculated there.

```read_abundance``` contains richness information and the script to perform a Mann-Whitney test in R on the reads mapping to a certain gene, testing treatment against control.

```samtools_depth``` gives the per base coverage of a certain gene by using samtools
