#counting occurrences of $gene and $contigs in columns 5 and 1 respectively in all cov files found
#printing information including previously applied parameters
#all output except for the last input gene will be printed by  FNR==1 { prt() }, last input file uses END { prt() }
#adding step to calculate duplicates? not tested yet

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
        printf "%s was found %d times on %d contigs with minimum coverage of %d in %s, inferring %d duplicates\n",
            gene, length(genes), length(contigs), threshold, fname, length(genes) - length(contigs)
        delete genes
        delete contigs
    }
    fname = FILENAME
}
' */*.cov

###thanks to Ed Morton on Stack Overflow
#https://stackoverflow.com/questions/49067678/print-output-of-user-defined-function-in-awk-gives-unexpected-token-error