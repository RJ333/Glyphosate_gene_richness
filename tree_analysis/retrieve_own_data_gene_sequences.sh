# most abundant or responding genera by 16S data

# should we use the respective gene from all available reference genomes for the tree?
Aminobacter
Aquamicrobium
Blastomonas
Brevundimonas
Caulobacter
Cupriavidus
Dokdonella
Ferrovibrio
Gallaecimonas
Hoeflea
Hyphomonas
Idiomarina
Limnohabitans
Loktanella
Massilia
Mesorhizobium
Methylotenera
Microbacterium
Nesiotobacter
Parvibaculum
Pseudolabrys
Pseudomonas
Pseudorhodobacter
Reyranella
Rhizobium
Rhodobacter
Sphingomonas
Sphingopyxis
Sphingorhabdus
Terrimonas
Thalassobaculum
uncultured_Caulobacteraceae
uncultured_Rhodobacteraceae
uncultured_Rhodospirillaceae



# phosphonate degrading

gene:	product2:

phnC	phosphate@import@atp@binding@protein@phnc
phnD	phosphate@import@protein@phnd@precursor
phnE	phosphate@import@permease@protein@phne
phnF	putative@transcriptional@regulator@phnf
phnG	alpha@d@ribose@1@methylphosphonate@5@triphosphate@synthase@subunit@phng
phnH	alpha@d@ribose@1@methylphosphonate@5@triphosphate@synthase@subunit@phnh
phnI	alpha@d@ribose@1@methylphosphonate@5@triphosphate@synthase@subunit@phni
phnJ	alpha@d@ribose@1@methylphosphonate@5@phosphate@c@p@lyase
phnK	putative@phosphonates@utilization@atp@binding@protein@phnk
phnL	alpha@d@ribose@1@methylphosphonate@5@triphosphate@synthase@subunit@phnl
phnM	alpha@d@ribose@1@methylphosphonate@5@triphosphate@diphosphatase
phnN	ribose@1%2c5@bisphosphate@phosphokinase@phnn
		aminoalkylphosphonic@acid@n@acetyltransferase
phnP	phosphoribosyl@1%2c2@cyclic@phosphodiesterase


# sarcosine oxidizing

product2: 

monomeric@sarcosine@oxidase
sarcosine@oxidase@subunit@beta
sarcosine@oxidase%2c@gamma@subunit@family
sarcosine@oxidase%2c@delta@subunit@family


# phosphorus sensing/starvation

gene:	product2: 

phoB	phosphate@regulon@transcriptional@regulatory@protein@phob
phoR	phosphate@regulon@sensor@protein@phor

# Pi starvation-inducible (psi) genes,
 
gene:	product2: 

psiF 	phosphate@starvation@inducible@protein@psif@precursor

#################################### pieces of code?


#################### this code can be used to adress amino acid sequences!

# prokka already provides the corresponding faa-sequences per annotation
# the file headers are not appropriate yet: need to be adjusted accordingly to above
# ## replace white lines and special characters in product names with "@"
# prokka_select$product2 <- gsub("_|-|'| |/|:|\\(|\\.|\\)|\\|", "@", prokka_select$product2)
# tolower

# additionally, the header should involve the sample names to distinguish contigs across samples


# put fasta output with header into one line
awk '{if (!/>/) {printf "%s",$0;next} else {printf "%s%s%s","\n",$0,"\t"}}' prokka.faa > prokka_1line.faa

# grep line and the following
grep -A 1 "phosphonate" prokka_1line.faa












#################### the code below can be used if you want to extract DNA substrings from the contigs

# read in the combined prokka data from all samples

print("reading prokka data ...")
prokka_tables <- fread("/data/Rene/glyph/prokka/prokka_all_modified.tsv",
  col.names = c("sample", "contig_id", "gene_start", "gene_end", "annotation"),
  colClasses = c("factor", "factor", "integer", "integer", "character"))

# splitting annotation-column
prokka_neat <- prokka_tables %>% 
  cSplit("annotation", sep = ";") %>%
  gather(Key, value, -c(sample, contig_id, gene_start, gene_end)) %>% 
  separate(value, c("Col", "Value"), sep = "=") %>% 
  select(-Key) %>% 
  filter(!(is.na(Col) & is.na(Value))) %>% 
  spread(Col, Value)
  
# combine information from columns "note" and "product"
print("combining columns product and note")
prokka_neat$product2 <- ifelse(prokka_neat$product == "hypothetical protein" & !is.na(prokka_neat$note), 
  prokka_neat$note, prokka_neat$product)

# calculate gene length
print("calculating gene length")
prokka_neat$gene_length <- prokka_neat$gene_end - prokka_neat$gene_start

# remove now unneeded columns "note" and "product"
# gene start and end is needed for bed file generation!
print("selecting important columns")
prokka_select <- prokka_neat[, c("sample", "contig_id", "gene_start", "gene_end", "gene_length",
								 "eC_number", "gene", "product2", "locus_tag")]

## unify product names , all to lower
print("unifying product annotations")
prokka_select$product2 <- tolower(prokka_select$product2)

## replace white lines and special characters in product names with "@"
prokka_select$product2 <- gsub("_|-|'| |/|:|\\(|\\.|\\)|\\|", "@", prokka_select$product2)

# unify gene numbering e.g. genE_01, genE_02 etc --> genE
print("unifying gene names")
prokka_select$gene <- gsub("_.*$","", prokka_select$gene)

# generate the pos.txt file for the awk script




# stackoverflow: https://stackoverflow.com/questions/45153483/how-to-use-info-on-substring-position-from-one-file-to-extract-substring-from-an

cut.awk

# We are reading two files: pos.txt and strings.txt
# NR is equal to FNR as long as we are reading the
# first file.
NR==FNR{
    pos[">"$1]=$2 # Store the startpoint in an array pos (indexed by $1)
    len[">"$1]=$4 # Store the length in an array len (indexed by $1)
    next # skip the block below for pos.txt
}

# This runs on every line of strings.txt
$1 in pos {
    # Extract a substring of $2 based on the position and length
    # stored above
    key=$1
    mod=substr($2,pos[key],len[key])
    $2=mod
    print # Print the modified line
}

Call it like this:

awk -f cut.awk pos.txt strings.txt

# One important thing to mention. substr() assumes strings to start at index 1 - in opposite to most programming languages where strings start at index 0. If the positions in pos.txt are 0 based, the substr() must become:

mod=substr($2,pos[key]+1,len[key])

# I recommend to test it with simplified, meaningful versions of:

pos.txt

foo  2  5  3    phnW  
bar  4  5  1    phnW
test 1  5  4    phnW

and strings.txt

>foo 123456  
>bar 123456
>non 123456

Output:

>foo 234
>bar 4