cd /data/Rene/glyph/trees_degradation/msa
mkdir phnH
cd phnH
cp  /data/Rene/glyph/trees_degradation/uniprot/sequences_modified/*phnH* .


# take all sequences from reference data
cat *.fasta >> all.phnH.reference.fasta


# take all phnH sequences from own data
egrep -h -A 1 "^>" /data/Rene/glyph/tree/prokka_modified/* |\
  egrep -A 1 "PhnH" | sed '/^--$/d' > own.phnH.faa

# perform cd-hit for reference data to reduce amount of sequences for tree
cd-hit -i all.phnH.reference.fasta -o phnH.ref.cdhit75.faa -c 1 -d 0

# removing duplicates (own data)
cd-hit-dup -i own.phnH.faa -o phnH.own.cdhit100.faa

# merge two one file with all sequences for alignment and tree
cat phnH.own.cdhit100.faa > phnH.all.faa
cat phnH.ref.cdhit75.faa >> phnH.all.faa

# multipe sequence alignment
mafft --auto phnH.all.faa > phnH.all.msa.faa

#convert to phylip
/data/Rene/glyph/trees_degradation/tree/fasta_to_phylip.py \
phnH.all.msa.faa phnH.all.msa.faa.phy

# remove unwanted characters
python2 /data/Rene/glyph/trees_degradation/tree/remove_special_characters.py phnH.all.msa.faa.phy \
  > phnH.all.msa.faa.clean.phy
  
cat phnH.all.msa.faa.clean.phy > phnH.all.cdhit95.msa.faa.clean.phy
  
  