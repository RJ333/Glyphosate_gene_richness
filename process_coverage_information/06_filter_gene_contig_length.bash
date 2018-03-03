###filter out genes and contigs which do not match the quality requirements
#not tested yet

genename="phnM"
min_gene_length="400"
min_contig_length="1000"

echo "gene length threshold is set to $min_gene_length"
echo "contig length threshold is set to $min_contig_length";
#maybe get rid of for loop, not really necessary here?
for length_file in */${genename}_contig_and_gene_length.txt
do
	echo $length_file;
	awk -v a="$min_gene_length" b="$min_contig_length" '$2 >= a && $4 >= b {print $0}' $length_file > ${length_file}_trimmed.txt
done
