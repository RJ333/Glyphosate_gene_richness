```bash
#################### this code can be used to adress amino acid sequences!

# put fasta output with header into one line
awk '{if (!/>/) {printf "%s",$0;next} else {printf "%s%s%s","\n",$0,"\n"}}' prokka.faa > prokka_1line.faa

# grep line and the following
grep -A 1 "phosphonate" prokka_1line.faa
```











#################### the code below can be used if you want to extract DNA substrings from the contigs

```r
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
```

```bash
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
```