#!/usr/bin/bash




# go to reference data and modify the fasta headers
cd /data/Rene/glyph/trees_degradation/uniprot_sequences/sequences_original

for file in *.fasta; 
do 
	echo $file; 
    python2 /data/Rene/git/Glyphosate_gene_richness/tree_analysis/modify_uniprot_headers.py $file >\
      ../sequences_modified/${file}; 
done

# move to msa folder and and copy modified uniprot sequences there
cd /data/Rene/glyph/trees_degradation/msa/phnC
cp ../../uniprot_sequences/sequences_modified/*phnC* .

# take all sequences from reference data
cat *.fasta >> all.phnC.reference.fasta

# perform cd-hit for reference data to reduce amount of sequences for tree
cd-hit -i all.phnC.reference.fasta -o phnC.ref.cdhit75.faa -c 0.75 -d 0

# removing duplicates (own data)
cd-hit-dup -i PhnC_glyphosate_metagenomes.faa -o phnC.own.cdhit100.faa

# merge two one file with all sequences for alignment and tree
cat phnC.own.cdhit100.faa > phnC.all.faa
cat phnC.ref.cdhit75.faa >> phnC.all.faa

# multipe sequence alignment
mafft --auto phnC.all.faa > phnC.all.msa.faa

#convert to phylip
/data/Rene/glyph/trees_degradation/tree/fasta_to_phylip.py \
phnC.all.msa.faa phnC.all.msa.faa.phy

# remove unwanted characters
python2 /data/Rene/glyph/trees_degradation/tree/remove_special_characters.py phnC.all.msa.faa.phy \
  > phnC.all.msa.faa.clean.phy
  
cat phnC.all.msa.faa.clean.phy > phnC.all.cdhit75.msa.faa.clean.phy
  
  
# build tree with RAXML  
raxmlHPC-PTHREADS-SSE3 \
  -T 28 \
  -p 12345 \
  -s phnC.all.msa.cdhit.clean.phy \
  -m PROTCATAUTO \
  -n phnC \
  -x 12345 \
  -f a \
  -N 100 \
  -w /data/glyphosate/trees/phnC/raxml
