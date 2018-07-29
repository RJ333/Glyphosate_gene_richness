# how to find the exact products for single copy genes

# dna gyrase a und b
head(subset(prokka_all, grepl("gyr", gene)))
table(droplevels(subset(prokka_all, grepl("dna", product2) & 
  grepl("gyrase", product2))$product2))
table(droplevels(subset(prokka_all, grepl("dna", product2) & 
  grepl("gyrase", product2))$gene))
subset(prokka_all, grepl("@@dna@gyrase", product2))

# phosphoglycerate@kinase

head(subset(prokka_all, grepl("phosphoglycerate@kinase", product2)))
table(droplevels(subset(prokka_all, grepl("phosphoglycerate@kinase", product2))$product2))
table(droplevels(subset(prokka_all, grepl("phosphoglycerate@kinase", product2))$gene))
subset(prokka_all, grepl("2@phosphoglycerate@kinase", product2))

# trna synthetases class I ( in verschiedenen kombis, nix gefunden)

head(subset(prokka_all, grepl("trna@synthetases", product2) & grepl("class", product2)))

# mraw methylase ( in verschiedenen kombis, nix für mra gefunden)
head(subset(prokka_all, grepl("mraW", gene)))
subset(prokka_all, grepl("mra", gene))$product2


table(droplevels(subset(prokka_all, grepl("methylase", product2) & 
  grepl("mra", product2))$product2))
table(droplevels(subset(prokka_all, grepl("dna", product2) & 
  grepl("gyrase", product2))$gene))
subset(prokka_all, grepl("methylase", product2))

# tig trigger factor
subset(prokka_all, grepl("trigger@factor", product2))
table(droplevels(subset(prokka_all, grepl("tig", gene))$product2))
table(droplevels(subset(prokka_all, grepl("trigger", product2))$product2))

# pyrG ctp@synthase
subset(prokka_all, grepl("ctp@synthase", product2))
table(droplevels(subset(prokka_all, grepl("pyrG", gene))$product2))
table(droplevels(subset(prokka_all, grepl("ctp@synthase", product2))$product2))


#------------------------tRNA synthase ligase synthetase (same meaning)

# don't forget to check with grep for e.g. "glycine" and "glycyl"!!!
# as usual hits are either leucine@@trna@ligase or leucyl@trna@synthetase



# glyS glycine trna ligase (glyS, glyQS, glyQ)
subset(prokka_all, grepl("glycine", product2))
table(droplevels(subset(prokka_all, grepl("glyS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("glycine", product2) & grepl("trna", product2))$product2))

# glyS isoleucine trna ligase
table(droplevels(subset(prokka_all, grepl("ileS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("isoleucine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("isoleucine", product2) & grepl("trna", product2))$gene))


# leuS leucine trna ligase
table(droplevels(subset(prokka_all, grepl("leuS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("leucine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("leucine", product2) & grepl("trna", product2))$gene))

# serS serine trna ligase
table(droplevels(subset(prokka_all, grepl("serS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("serine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("serine", product2) & grepl("trna", product2))$gene))

# thrS threonine trna ligase
table(droplevels(subset(prokka_all, grepl("thrS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("threonine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("threonine", product2) & grepl("trna", product2))$gene))

# valS valine trna ligase
table(droplevels(subset(prokka_all, grepl("valS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("valine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("valine", product2) & grepl("trna", product2))$gene))

# cysS cysteine trna ligase
table(droplevels(subset(prokka_all, grepl("cysS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("cysteine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("cysteine", product2) & grepl("trna", product2))$gene))

# tyrS tyrosine trna ligase
table(droplevels(subset(prokka_all, grepl("tyrS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("tyrosine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("tyrosine", product2) & grepl("trna", product2))$gene))

# hisS histidine trna ligase
table(droplevels(subset(prokka_all, grepl("hisS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("histidine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("histidine", product2) & grepl("trna", product2))$gene))

# pheS pheT phenylalanine trna ligase
table(droplevels(subset(prokka_all, grepl("pheS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("pheT", gene))$product2))
table(droplevels(subset(prokka_all, grepl("phenylalanine", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("phenylalanine", product2) & grepl("trna", product2))$gene))

# aspS aspartate trna ligase
table(droplevels(subset(prokka_all, grepl("aspS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("aspartate", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("aspartate", product2) & grepl("trna", product2))$gene))

# proS proline trna ligase
table(droplevels(subset(prokka_all, grepl("proS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("proline", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("proline", product2) & grepl("trna", product2))$gene))

# alaS alanine trna ligase
table(droplevels(subset(prokka_all, grepl("alaS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("alan", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("alanine", product2) & grepl("trna", product2))$gene))

#--------------------------------------------------------------------

# peptide chain release factor
table(droplevels(subset(prokka_all, grepl("prfA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("peptide@chain", product2) & grepl("factor", product2))$product2))
table(droplevels(subset(prokka_all, grepl("peptide@chain", product2) & grepl("factor", product2))$gene))
subset(prokka_all, grepl("peptide@chain", product2) & grepl("factor", product2))$gene

# transcription termination factor
table(droplevels(subset(prokka_all, grepl("nusG", gene))$product2))
table(droplevels(subset(prokka_all, grepl("transcription@termination", product2) & grepl("factor", product2))$product2))
table(droplevels(subset(prokka_all, grepl("transcription@termination", product2) & grepl("factor", product2))$gene))
subset(prokka_all, grepl("transcription@termination", product2) & grepl("factor", product2))$gene

# transcription termination factor
table(droplevels(subset(prokka_all, grepl("nusA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("nusa", product2)))$product2))
table(droplevels(subset(prokka_all, grepl("transcription@termination", product2) & grepl("factor", product2))$gene))
subset(prokka_all, grepl("transcription@termination", product2) & grepl("factor", product2))$gene

# secA E Y preprotein translocase
table(droplevels(subset(prokka_all, grepl("secA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("translocase", product2) & grepl("preprotein", product2))$product2))
table(droplevels(subset(prokka_all, grepl("translocase", product2) & grepl("preprotein", product2))$gene))
subset(prokka_all, grepl("translocase", product2) & grepl("preprotein", product2))$gene

# dnaA replication chromosomal
table(droplevels(subset(prokka_all, grepl("dnaA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("chromosomal", product2) & grepl("replication", product2))$product2))
table(droplevels(subset(prokka_all, grepl("chromosomal", product2) & grepl("replication", product2))$gene))
subset(prokka_all, grepl("chromosomal", product2) & grepl("replication", product2))$gene

# trmU trna methylamin
table(droplevels(subset(prokka_all, grepl("trmU", gene))$product2))
table(droplevels(subset(prokka_all, grepl("methylamin", product2) & grepl("trna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("methylamin", product2) & grepl("trna", product2))$gene))
subset(prokka_all, grepl("methylamin", product2) & grepl("trna", product2))$gene

# smpB trna binding@protein
table(droplevels(subset(prokka_all, grepl("smpB", gene))$product2))
table(droplevels(subset(prokka_all, grepl("binding@protein", product2) & grepl("ssr", product2))$product2))
table(droplevels(subset(prokka_all, grepl("binding@protein", product2) & grepl("ssr", product2))$gene))
subset(prokka_all, grepl("binding@protein", product2) & grepl("ssr", product2))$gene

# grpE trna grpe
table(droplevels(subset(prokka_all, grepl("grpE", gene))$product2))
table(droplevels(subset(prokka_all, grepl("grpe", product2))$product2))
table(droplevels(subset(prokka_all, grepl("grpe", product2))$gene))
subset(prokka_all, grepl("grpe", product2))$gene

# ychF (vermutlich falsch, nur ribosome@binding@atpase@ychf gefunden)
table(droplevels(subset(prokka_all, grepl("ychF", gene))$product2))
table(droplevels(subset(prokka_all, grepl("binding@protein", product2) & grepl("GTP", product2))$product2))
table(droplevels(subset(prokka_all, grepl("binding@protein", product2) & grepl("GTP", product2))$gene))
subset(prokka_all, grepl("binding@protein", product2) & grepl("GTP", product2))$gene

# infABC translation initiation factor
table(droplevels(subset(prokka_all, grepl("inf", gene))$product2))
table(droplevels(subset(prokka_all, grepl("translation@initiation", product2) & grepl("factor", product2))$product2))
table(droplevels(subset(prokka_all, grepl("translation@initiation", product2) & grepl("factor", product2))$gene))
subset(prokka_all, grepl("translation@initiation", product2) & grepl("factor", product2))$gene

# frr ribosome recycling factor
table(droplevels(subset(prokka_all, grepl("frr", gene))$product2))
table(droplevels(subset(prokka_all, grepl("recycling", product2) & grepl("factor", product2))$product2))
table(droplevels(subset(prokka_all, grepl("recycling", product2) & grepl("factor", product2))$gene))
subset(prokka_all, grepl("recycling", product2) & grepl("factor", product2))$gene

# dnlj DNA ligase NAD dependent (unsichere Angaben, ligA, ligB und ykoU)
table(droplevels(subset(prokka_all, grepl("dnlj", gene))$product2))
table(droplevels(subset(prokka_all, grepl("dna", product2) & grepl("ligase", product2))$product2))
table(droplevels(subset(prokka_all, grepl("dna", product2) & grepl("ligase", product2))$gene))
subset(prokka_all, grepl("dna", product2) & grepl("ligase", product2))$gene

# uvrB excinuclease abc nuclease, uvrabc@system
table(droplevels(subset(prokka_all, grepl("uvrB", gene))$product2))
table(droplevels(subset(prokka_all, grepl("abc", product2) & grepl("uvr", product2))$product2))
table(droplevels(subset(prokka_all, grepl("abc", product2) & grepl("nuclease", product2))$gene))
subset(prokka_all, grepl("abc", product2) & grepl("nuclease", product2))$gene

# dnaN dna@polymerase@iii subunit
table(droplevels(subset(prokka_all, grepl("dnaN", gene))$product2))
table(droplevels(subset(prokka_all, grepl("beta", product2) & grepl("dna@polymerase@iii", product2))$product2))
table(droplevels(subset(prokka_all, grepl("beta", product2) & grepl("dna@polymerase@iii", product2))$gene))
subset(prokka_all, grepl("dna@polymerase@iii", product2) & grepl("subunit", product2))$gene

# tsf translation@elongation
table(droplevels(subset(prokka_all, grepl("tsf", gene))$product2))
table(droplevels(subset(prokka_all, grepl("factor@ts", product2) & grepl("elongation", product2))$product2))
table(droplevels(subset(prokka_all, grepl("factor@ts", product2) & grepl("elongation", product2))$gene))
subset(prokka_all, grepl("elongation", product2) & grepl("factor@ts", product2))$gene

# fmt methionyl formyl
table(droplevels(subset(prokka_all, grepl("fmt", gene))$product2))
table(droplevels(subset(prokka_all, grepl("formyl", product2) & grepl("methionyl", product2))$product2))
table(droplevels(subset(prokka_all, grepl("formyl", product2) & grepl("methionyl", product2))$gene))
subset(prokka_all, grepl("methionyl", product2) & grepl("formyl", product2))$gene

# infABC translation initiation factor
table(droplevels(subset(prokka_all, grepl("infC", gene))$product2))
table(droplevels(subset(prokka_all, grepl("translation@initiation", product2) & grepl("factor@if@3", product2))$product2))
table(droplevels(subset(prokka_all, grepl("translation@initiation", product2) & grepl("factor@if@3", product2))$gene))
subset(prokka_all, grepl("translation@initiation", product2) & grepl("factor@if@3", product2))$gene

# dnaG dna primase
table(droplevels(subset(prokka_all, grepl("dnaG", gene))$product2))
table(droplevels(subset(prokka_all, grepl("primase", product2) & grepl("dna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("primase", product2) & grepl("dna", product2))$gene))
subset(prokka_all, grepl("dna", product2) & grepl("primase", product2))$gene

# dnaX dna
table(droplevels(subset(prokka_all, grepl("dnaX", gene))$product2))
table(droplevels(subset(prokka_all, grepl("tau", product2) & grepl("dna@polymerase", product2))$product2))
table(droplevels(subset(prokka_all, grepl("chaperone", product2) & grepl("protein", product2))$gene))
subset(prokka_all, grepl("protein", product2) & grepl("chaperone", product2))$gene

# lepA gtp binding elongation factor 4
table(droplevels(subset(prokka_all, grepl("lepA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("binding", product2) & grepl("gtp", product2))$product2))
table(droplevels(subset(prokka_all, grepl("binding", product2) & grepl("gtp", product2))$gene))
subset(prokka_all, grepl("gtp", product2) & grepl("binding", product2))$gene

# recA protein 
table(droplevels(subset(prokka_all, grepl("rec", gene))$product2))
table(droplevels(subset(prokka_all, grepl("reca", product2) & grepl("protein", product2))$product2))
table(droplevels(subset(prokka_all, grepl("reca", product2))$product2))
subset(prokka_all, grepl("protein", product2) & grepl("reca", product2))$gene

# rpo dna polymerase elongation factor 4
table(droplevels(subset(prokka_all, grepl("rpo[ABC]", gene))$product2))
table(droplevels(subset(prokka_all, grepl("polymerase", product2) & grepl("dna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("polymerase", product2) & grepl("beta", product2))$product2))
table(droplevels(subset(prokka_all, grepl("polymerase", product2) & grepl("dna", product2))$gene))
subset(prokka_all, grepl("dna", product2) & grepl("polymerase", product2))$gene

# ribonuclease iii 3
table(droplevels(subset(prokka_all, grepl("ribonuclease", product2))$product2))
table(droplevels(subset(prokka_all, grepl("ribonuclease", product2) & grepl("3", product2))$gene))
subset(prokka_all, grepl("iii", product2) & grepl("ribonuclease", product2))$gene

# guanylate kinase
table(droplevels(subset(prokka_all, grepl("gmk", gene))$product2))
table(droplevels(subset(prokka_all, grepl("guanylate", product2) & grepl("kinase", product2))$product2))
table(droplevels(subset(prokka_all, grepl("guanylate", product2) & grepl("kinase", product2))$gene))
subset(prokka_all, grepl("kinase", product2) & grepl("guanylate", product2))$gene

# lysidine trna synthetase
table(droplevels(subset(prokka_all, grepl("tilS", gene))$product2))
table(droplevels(subset(prokka_all, grepl("lysidine", product2) & grepl("rna", product2))$product2))
table(droplevels(subset(prokka_all, grepl("lysidine", product2) & grepl("rna", product2))$gene))
subset(prokka_all, grepl("rna", product2) & grepl("lysidine", product2))$gene

# gtpase obg cgtA
table(droplevels(subset(prokka_all, grepl("cgtA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("gtpase", product2) & grepl("obg", product2))$product2))
table(droplevels(subset(prokka_all, grepl("gtpase", product2) & grepl("obg", product2))$gene))
subset(prokka_all, grepl("obg", product2) & grepl("gtpase", product2))$gene

# gtpase ribosome engA
table(droplevels(subset(prokka_all, grepl("engA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("gtpase", product2) & grepl("ribosome", product2))$product2))
table(droplevels(subset(prokka_all, grepl("gtpase", product2) & grepl("ribosome", product2))$gene))
subset(prokka_all, grepl("ribosome", product2) & grepl("gtpase", product2))$gene

# gtp@binding era era
table(droplevels(subset(prokka_all, grepl("era", gene))$product2))
table(droplevels(subset(prokka_all, grepl("gtp", product2) & grepl("era", product2))$product2))
table(droplevels(subset(prokka_all, grepl("gtp", product2) & grepl("era", product2))$gene))
subset(prokka_all, grepl("era", product2) & grepl("gtp", product2))$gene

# coa kinase dephospho coaE
table(droplevels(subset(prokka_all, grepl("coaE", gene))$product2))
table(droplevels(subset(prokka_all, grepl("coa", product2) & grepl("dephospho", product2))$product2))
table(droplevels(subset(prokka_all, grepl("coa", product2) & grepl("dephospho", product2))$gene))
subset(prokka_all, grepl("dephospho", product2) & grepl("coa@kinase", product2))$gene

# coa kinase metalloprotein coaE
table(droplevels(subset(prokka_all, grepl("ybey", gene))$product2))
table(droplevels(subset(prokka_all, grepl("ybe", product2) & grepl("metalloprotein", product2))$product2))
table(droplevels(subset(prokka_all, grepl("ybe", product2) & grepl("endoribonuclease", product2))$product2))
subset(prokka_all, grepl("endoribonuclease", product2) & grepl("ybe", product2))$gene

# coa kinase antiporter nhaD
table(droplevels(subset(prokka_all, grepl("nhaD", gene))$product2))
table(droplevels(subset(prokka_all, grepl("nhad", product2) & grepl("antiporter", product2))$product2))
table(droplevels(subset(prokka_all, grepl("nhad", product2) & grepl("antiporter", product2))$gene))
subset(prokka_all, grepl("antiporter", product2) & grepl("", product2))$gene

# coa kinase antiporter ftsY
table(droplevels(subset(prokka_all, grepl("ftsY", gene))$product2))
table(droplevels(subset(prokka_all, grepl("signal@recognition@particle", product2) & grepl("recept", product2))$product2))
table(droplevels(subset(prokka_all, grepl("signal", product2) & grepl("recognition", product2))$gene))
subset(prokka_all, grepl("recognition", product2) & grepl("", product2))$gene

# coa kinase antiporter ftsY
table(droplevels(subset(prokka_all, grepl("ffh", gene))$product2))
table(droplevels(subset(prokka_all, grepl("signal@recognition", product2) & grepl("particle@protein", product2))$product2))
table(droplevels(subset(prokka_all, grepl("signal", product2) & grepl("recognition", product2))$gene))
subset(prokka_all, grepl("recognition", product2) & grepl("", product2))$gene

# rbfA: ribosome-binding factor A
table(droplevels(subset(prokka_all, grepl("rbfA", gene))$product2))
table(droplevels(subset(prokka_all, grepl("ribosome", product2) & grepl("binding", product2))$product2))
table(droplevels(subset(prokka_all, grepl("ribosome", product2) & grepl("binding", product2))$gene))
subset(prokka_all, grepl("binding", product2) & grepl("ribosome", product2))$gene


table(droplevels(subset(prokka_all, grepl("trna", product2) & grepl("phenyl", product2)& grepl("beta", product2))$product2))
sum(table(droplevels(subset(prokka_all, grepl("trna", product2) & grepl("phenyl", product2)& grepl("beta", product2))$product2)))
