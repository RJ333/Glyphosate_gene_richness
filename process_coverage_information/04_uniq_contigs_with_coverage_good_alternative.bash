awk '
BEGIN {
    gene = "phnM"
    threshold = "5"
}
FNR==1 { prt() }
{
    genes[$5]
    contigs[$1]
}
END { prt() }
function prt() {
    if (fname != "") {
        printf "%s was found %d times on %d contigs with minimum coverage of %d in %s\n",
            gene, length(genes), length(contigs), threshold, fname
        delete genes
        delete contigs
    }
    fname = FILENAME
}
' */*.cov