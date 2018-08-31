## TO DO
# add variables
# include sequence length in distribution plot
# check number of sequences for prokka and ref > 0



cd /data/Rene/glyph/trees_degradation/msa
mkdir phnJ
cd phnJ
cp  /data/Rene/glyph/trees_degradation/uniprot/sequences_modified/*phnJ* .


# take all sequences from reference data
cat *.fasta >> all.phnJ.reference.fasta


# take all phnJ sequences from own data
egrep -h -A 1 "^>" /data/Rene/glyph/tree/prokka_modified/* |\
  egrep -i -A 1 "phosphate@c-p@lyase" | sed '/^--$/d' > own.phnJ.faa

# perform cd-hit for reference data to reduce amount of sequences for tree
cd-hit -i all.phnJ.reference.fasta -o phnJ.ref.cdhit75.faa -c 1 -d 0

# removing duplicates (own data)
cd-hit-dup -i own.phnJ.faa -o phnJ.own.cdhit100.faa

# merge two one file with all sequences for alignment and tree
cat phnJ.own.cdhit100.faa > phnJ.all.faa
cat phnJ.ref.cdhit75.faa >> phnJ.all.faa

# multipe sequence alignment
mafft --auto phnJ.all.faa > phnJ.all.msa.faa

#convert to phylip
/data/Rene/glyph/trees_degradation/tree/fasta_to_phylip.py \
phnJ.all.msa.faa phnJ.all.msa.faa.phy

# remove unwanted characters
python2 /data/Rene/glyph/trees_degradation/tree/remove_special_characters.py phnJ.all.msa.faa.phy \
  > phnJ.all.msa.faa.clean.phy
  
cat phnJ.all.msa.faa.clean.phy > phnJ.all.cdhit100.msa.faa.clean.phy
  
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
  -w /data/glyphosate/trees/phnJ/raxml