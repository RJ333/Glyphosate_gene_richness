# plots based on prokka_all (contains reads per product per taxonomy) 
# and product_reads2 (contains richness per product)

# data is prepared in 03_calculate_and_plot_gene_richness.R
#and 01_funct_tax_data_merging.bash

library(ggplot2)
library(dplyr)
# plot taxonomic richness

ggplot(tax_rich_omics, aes(x = days, colour = treatment))+
  geom_line(aes(y = richness), linetype = 1) +
  geom_line(aes(y = richness_rel), linetype = 2) +
  facet_wrap(~ treatment, nrow = 2)
  
  
  
# plot gene richness
sarc <- subset(gene_richness2, grepl("arcosin", product2))
sarc <- subset(sarc, grepl("oxidase", product2))
sarc <- subset(sarc, grepl("delta", product2))

# visualization per product
ggplot(sarc, aes (x = new_day, colour = product2, fill = contig_id))+
  geom_bar(stat = "identity", position = position_dodge(), aes(y = product_rpm))+
  #geom_line(aes(y = gene_richness_relative), size = 1.5)+
  #geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel, colour = "black"), linetype = 2, size = 2, alpha = 0.7) +
  theme(legend.position="none")+
  facet_wrap(~treatment, nrow = 2)

# cumulated product's visualization  
ggplot(sarc, aes(colour = product2))+
  geom_bar(stat = "identity", aes(x = new_day - 0.75, y = product_rpm), fill = "black", width = 1)+
  geom_bar(stat = "identity", aes(x = new_day + 0.75, y = gene_richness_relative), fill = "blue", width = 1)+
  geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel, colour = "black"), linetype = 2) + 
  facet_wrap(~treatment, nrow = 2) 


phn_op <- subset(prokka_all, treatment == "glyph" & grepl("phn[C-N]", gene))
sox_all <- subset(prokka_all, treatment == "glyph" & grepl("sarcosine@oxidase", product2))


test <- sox_all[with(sox_all, order(sox_all$product2, sox_all$sample, -sox_all$average_coverage)),]

test <- table(prokka_all$average_coverage) 
# with this we can check what average coverage is neccessary to
# generate a contig and check, how close our interesting genes are to that
# lowest value  for phn_op is 2.296, for sox_all is 2.55, 
# but we have at least ~7500 contigs with lower average coverage
sum(head(test, 4600))

phn_op_rich <- subset(gene_richness2, treatment == "glyph" & grepl("phn[C-N]", gene))
sox_all_rich <- subset(gene_richness2, treatment == "glyph" & grepl("sarcosine@oxidase", product2))



#Anteil von phn Genen und sox Genen ohne taxonomische Zuordnung?
table(is.na(phn_op$genus))
FALSE  TRUE 
  126   958 
  
table(is.na(sox_all$genus))
FALSE  TRUE 
   71   415 


#------------------------------------------------------------------ assessment of single copy and housekeeping genes

   
# what genes do we have from Wu et al marker genes? 
# additionally Mads Albertsen et al, "recovery of genomes", 2013, supplementary Table S4
levels(droplevels(subset(gene_richness2, grepl("rpl", gene))$gene))
levels(droplevels(subset(gene_richness2, grepl("rps", gene))$gene))   

table(droplevels(subset(gene_richness2, grepl("rps", gene))$gene))
table(droplevels(subset(gene_richness2, grepl("rpl", gene))$gene))

table(droplevels(subset(prokka_all, grepl("rpl", gene))$gene))
table(droplevels(subset(prokka_all, grepl("rpm", gene))$gene))
table(droplevels(subset(prokka_all, grepl("rps", gene))$gene))


table(droplevels(subset(prokka_all, grepl("era", gene))$product2))
table(droplevels(subset(prokka_all, grepl("trna", product2) & grepl("class", product2))$product2))

table(droplevels(subset(prokka_all, grepl("trna", product2) & grepl("ligase", product2))$product2))

#translation table product2 gene
gene_product2 <- unique(gene_richness2[, c(1, 2)])


# get richness values per sample for selected large ribosomal proteins

pat_rpl <- c("50s@ribosomal@protein@l1",
			"50s@ribosomal@protein@l2",
			"50s@ribosomal@protein@l3",
			"50s@ribosomal@protein@l4",
			"50s@ribosomal@protein@l5",
			"50s@ribosomal@protein@l6",
			"50s@ribosomal@protein@l7@l12",
			"50s@ribosomal@protein@l10",
			"50s@ribosomal@protein@l11",
			"50s@ribosomal@protein@l13",
			"50s@ribosomal@protein@l14",
			"50s@ribosomal@protein@l15",
			"50s@ribosomal@protein@l16",
			"50s@ribosomal@protein@l17",
			"50s@ribosomal@protein@l18",
			"50s@ribosomal@protein@l19",
			"50s@ribosomal@protein@l20",
			"50s@ribosomal@protein@l21",
			"50s@ribosomal@protein@l22",
			"50s@ribosomal@protein@l23",
			"50s@ribosomal@protein@l24",
			"50s@ribosomal@protein@l25",
			"50s@ribosomal@protein@l27",
			"50s@ribosomal@protein@l28",
			"50s@ribosomal@protein@l29",
			"50s@ribosomal@protein@l32",
			"50s@ribosomal@protein@l34",
			"50s@ribosomal@protein@l35")

query_strings_rpl <- paste0("^", pat_rpl, "$")

rpl_A1 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A1")$product2)))
rpl_A2 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A2")$product2)))
rpl_A3 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A3")$product2)))
rpl_A4 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A4")$product2)))
rpl_A5 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A5")$product2)))
rpl_A6 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A6")$product2)))
rpl_A7 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "A7")$product2)))
rpl_B8 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "B8")$product2)))
rpl_B9 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "B9")$product2)))
rpl_B10 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rpl, collapse = "|"), product2) & sample == "B10")$product2)))

rpl <- Reduce(function(x,y) merge(x = x, y = y, by = "Var1"), 
  list(rpl_A1, rpl_A2, rpl_A3, rpl_A4, rpl_A5, rpl_A6, rpl_A7, rpl_B8, rpl_B9, rpl_B10))
names(rpl) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
rpl$group <- c(rep("rpl", nrow(rpl)))

# get richness values per sample for selected small ribosomal proteins

pat_rps <- c("30s@ribosomal@protein@s2",
			"30s@ribosomal@protein@s3",
			"30s@ribosomal@protein@s4",
			"30s@ribosomal@protein@s5",
			"30s@ribosomal@protein@s6",
			"30s@ribosomal@protein@s7",
			"30s@ribosomal@protein@s8",
			"30s@ribosomal@protein@s9",
			"30s@ribosomal@protein@s10",
			"30s@ribosomal@protein@s11",
			"30s@ribosomal@protein@s12",
			"30s@ribosomal@protein@s13",
			"30s@ribosomal@protein@s15",
			"30s@ribosomal@protein@s16",
			"30s@ribosomal@protein@s17",
			"30s@ribosomal@protein@s18",
			"30s@ribosomal@protein@s19",
			"30s@ribosomal@protein@s20")

query_strings_rps <- paste0("^", pat_rps, "$")

rps_A1 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A1")$product2)))
rps_A2 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A2")$product2)))
rps_A3 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A3")$product2)))
rps_A4 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A4")$product2)))
rps_A5 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A5")$product2)))
rps_A6 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A6")$product2)))
rps_A7 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "A7")$product2)))
rps_B8 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "B8")$product2)))
rps_B9 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "B9")$product2)))
rps_B10 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_rps, collapse = "|"), product2) & sample == "B10")$product2)))

rps <- Reduce(function(x,y) merge(x = x, y = y, by = "Var1"), 
  list(rps_A1, rps_A2, rps_A3, rps_A4, rps_A5, rps_A6, rps_A7, rps_B8, rps_B9, rps_B10))
names(rps) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
rps$group <- c(rep("rps", nrow(rps)))


# get richness values per sample for selected tRNA synthases/ligases

pat_trna_ligase <- c("glycine@@trna@ligase@alpha@subunit",
						"glycine@@trna@ligase@beta@subunit",
						"glycine@@trna@ligase",
						"glycyl@trna@synthetase@subunit@alpha",
						"glycyl@trna@synthetase@subunit@beta",
						"isoleucine@@trna@ligase",
						"isoleucyl@trna@synthetase",
						"leucine@@trna@ligase",
						"leucyl@trna@synthetase",
						"leucine@@trna@ligase@subunit@alpha",
						"leus@bact@@leucine@@trna@ligase",
						"serine@@trna@ligase",
						"seryl@trna@synthetase",
						"threonine@@trna@ligase",
						"threonine@@trna@ligase@1",
						"threonine@@trna@ligase@2",
						"threonyl@trna@synthetase",
						"valine@@trna@ligase",
						"vals@@valine@@trna@ligase",
						"valyl@trna@synthetase",
						"cysteine@@trna@ligase",
						"cysteinyl@trna@synthetase",
						"tyrosine@@trna@ligase",
						"tyrosine@@trna@ligase@1",
						"tyrosyl@trna@synthetase",
						"histidine@@trna@ligase",
						"hiss@@histidine@@trna@ligase",
						"histidyl@trna@synthetase",
						"phenylalanine@@trna@ligase@alpha@subunit",
						"phes@@phenylalanine@@trna@ligase%2c@alpha@subunit",
						"phenylalanine@@trna@ligase@beta@subunit",
						"phet@bact@@phenylalanine@@trna@ligase%2c@beta@subunit",
						"phenylalanyl@trna@synthetase@subunit@beta",
						"phenylalanyl@trna@synthetase@subunit@alpha",
						"arginine@@trna@ligase",
						"arginyl@trna@synthetase",
						"methionyl@trna@synthetase",
						"methionine@@trna@ligase",
						"alanine@@trna@ligase",
						"alas@@alanine@@trna@ligase",
						"alanyl@trna@synthetase",
						"trna@ile@@lysidine@synthase",
						"lysidine@tils@n@@trna@ile@@lysidine@synthetase",
						"trna@ile@@lysidine@synthetase",
						"aspartate@@trna@ligase",
						"aspartyl@trna@synthetase",
						"proline@@trna@ligase",
						"prolyl@trna@synthetase",
						"trna@pseudouridine@synthase@b")

						
query_strings_trna_ligase <- paste0("^", pat_trna_ligase, "$")
						
trna_ligase_A1 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A1")$product2)))
trna_ligase_A2 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A2")$product2)))
trna_ligase_A3 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A3")$product2)))
trna_ligase_A4 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A4")$product2)))
trna_ligase_A5 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A5")$product2)))
trna_ligase_A6 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A6")$product2)))
trna_ligase_A7 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "A7")$product2)))
trna_ligase_B8 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "B8")$product2)))
trna_ligase_B9 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "B9")$product2)))
trna_ligase_B10 <- as.data.frame(table(droplevels(subset(prokka_all, grepl(paste(query_strings_trna_ligase, collapse = "|"), product2) & sample == "B10")$product2)))

trna_ligase <- Reduce(function(x,y) merge(x = x, y = y, by = "Var1", all = TRUE), 
  list(trna_ligase_A1, trna_ligase_A2, trna_ligase_A3, trna_ligase_A4, trna_ligase_A5, trna_ligase_A6, trna_ligase_A7, trna_ligase_B8, trna_ligase_B9, trna_ligase_B10))
names(trna_ligase) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
trna_ligase[is.na(trna_ligase)] <- 0

						
group_vector_trna_ligase <- as.data.frame(cbind(pat_trna_ligase, c("trna_glyc1",
							"trna_glyc2",
							"trna_glyc",
							"trna_glyc1",
							"trna_glyc2",
							"trna_ile",
							"trna_ile",
							"trna_leu",
							"trna_leu",
							"trna_leu1",
							"trna_leu",
							"trna_ser",
							"trna_ser",
							"trna_thr",
							"trna_thr1",
							"trna_thr2",
							"trna_thr",
							"trna_val",
							"trna_val",
							"trna_val",
							"trna_cys",
							"trna_cys",
							"trna_tyr",
							"trna_tyr1",
							"trna_tyr",
							"trna_his",
							"trna_his",
							"trna_his",
							"trna_phe1",
							"trna_phe1",
							"trna_phe2",
							"trna_phe2",
							"trna_phe2",
							"trna_phe1",
							"trna_arg",
							"trna_arg",
							"trna_met",
							"trna_met",
							"trna_ala",
							"trna_ala",
							"trna_ala",
							"trna_til",
							"trna_til",
							"trna_til",
							"trna_asp",
							"trna_asp",
							"trna_pro",
							"trna_pro",
							"trna_psu")))

trna_ligase2 <- merge(trna_ligase, group_vector_trna_ligase, by.x = "product2", by.y = "pat_trna_ligase")			
names(trna_ligase2) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10", "group")
trna_ligase3 <- trna_ligase2[,-1]

trna_ligase4 <- aggregate(. ~ V2, data = trna_ligase3, sum)
trna_ligase4$group <- trna_ligase4$V2
names(trna_ligase4) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10", "group")

phn_op_A1 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A1")$product2)))
phn_op_A2 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A2")$product2)))
phn_op_A3 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A3")$product2)))
phn_op_A4 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A4")$product2)))
phn_op_A5 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A5")$product2)))
phn_op_A6 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A6")$product2)))
phn_op_A7 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "A7")$product2)))
phn_op_B8 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "B8")$product2)))
phn_op_B9 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "B9")$product2)))
phn_op_B10 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("phn[CDEGHIJKLMN]", gene) & sample == "B10")$product2)))

phn_op <- Reduce(function(x,y) merge(x = x, y = y, by = "Var1", all = TRUE), 
  list(phn_op_A1, phn_op_A2, phn_op_A3, phn_op_A4, phn_op_A5, phn_op_A6, phn_op_A7, phn_op_B8, phn_op_B9, phn_op_B10))
names(phn_op) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
phn_op[is.na(phn_op)] <- 0
phn_op$group <- c(rep("phn_op", nrow(phn_op)))

single_copy_richness <- rbind(rps, rpl) #phn_op, trna_ligase4)

# weitere Gruppen erstellen für mittelwert (referenz) gegen mittelwert(probe) ?
# wie willkürlich ist die auswahl für Referenz und Probe?
# weitere referenzen hinzufügen, checken, was die standardanzahl ist, rest entfernen, begründung dafür?
# was sagt die Abundanz?
# sox gene hinzufügen
# alle phn gene gefunden?
# welche referenzgene sind nutzlos?


single_copy_richness2 <- subset(single_copy_richness, A1 > 10 & A1 < 25)
single_copy_richness2 <- single_copy_richness

hist(single_copy_richness2$A1, breaks = 10)

library(reshape2)
single_copy_richness2_molten <- melt(single_copy_richness2, id = c("product2", "group"))
meta_single_copy_richness2_molten <- merge(single_copy_richness2_molten, meta_data_omics_small, by.x = "variable", by.y ="sample")
meta_single_copy_richness2_molten$group <- as.factor(meta_single_copy_richness2_molten$group)
library(data.table)
bla <- meta_single_copy_richness2_molten
setDT(bla)
bla[,value_relative := value/value[new_day == 0]*100, by = .(product2, group, treatment)]

gd <- bla %>% 
        group_by(variable, treatment, new_day, group) %>% 
        summarise(value = mean(value))
#subset(bla, group == "rps" | group == "rpl" | group == "phn_op") %>% 

#displaying species and rpl rps richness with two y axes
ggplot(data = subset(bla, treatment == "glyph"), aes(x = new_day, y = value*4, colour = group, group = group)) +
  #geom_line(aes(group = product2), alpha = 0.4, size = 1)+
  geom_point(data = subset(gd, treatment == "glyph"), size = 3, alpha = 1, colour = "black") +
  geom_line(data = subset(gd, treatment == "glyph"), size = 0.8, alpha = 1) +
  #stat_summary(fun.y="mean", geom="line") +
  geom_point(data = subset(tax_rich, treatment == "glyph"), size = 3, aes(x = days, y = richness)) +
  geom_line(data = subset(tax_rich, treatment == "glyph"), size = 0.8, aes(x = days, y = richness)) +
  scale_y_continuous(name = "Species richness based on 16S rRNA gene amplicons", 
						sec.axis = sec_axis(~./4, name = "single copy gene richness")) +
  scale_colour_manual(values = c("group_tax" = "black", "rpl" = "darkgreen", "rps" = "lightgreen"),
						name = "Gene groups",
						breaks = c("group_tax", "rpl", "rps"),
						labels = c("16S rRNA gene\n (averaged, amplicons)", 
									"50S ribosomal proteins\n (averaged, metagenome)", 
									"30S ribosomal proteins\n (averaged, metagenome)")) +
	labs(x = "Days")
	
ggsave("richness_plot_y_axes.png", width = 12, height = 10)	
#displaying species and rpl rps richness on shared y axis	
ggplot(data = subset(bla, treatment == "glyph"), aes(x = new_day, y = value, colour = group, group = group)) +
  #geom_line(aes(group = product2), alpha = 0.4, size = 1)+
  geom_point(data = subset(gd, treatment == "glyph"), size = 3, alpha = 1, colour = "black") +
  geom_line(data = subset(gd, treatment == "glyph"), size = 0.8, alpha = 1) +
  #stat_summary(fun.y="mean", geom="line") +
  geom_point(data = subset(tax_rich, treatment == "glyph"), size = 3, aes(x = days, y = richness)) +
  geom_line(data = subset(tax_rich, treatment == "glyph"), size = 0.8, aes(x = days, y = richness)) +
  scale_y_continuous(name = "Species or gene richness") +
  scale_colour_manual(values = c("group_tax" = "black", "rpl" = "darkgreen", "rps" = "lightgreen"),
						name = "Gene groups",
						breaks = c("group_tax", "rpl", "rps"),
						labels = c("16S rRNA gene\n (averaged, amplicons)", 
									"50S ribosomal proteins\n (averaged, metagenome)", 
									"30S ribosomal proteins\n (averaged, metagenome)")) +
	labs(x = "Days")	
ggsave("richness_plot_shared.png", width = 12, height = 10)
  #facet_wrap(~treatment, nrow = 2)

  
# plotting sequencing depth of metagenomes
ggplot(data = subset(meta_omics_small, treatment == "glyph"), aes(x = new_day, y = total_reads, colour = treatment)) +
  geom_point(size = 3) +
  geom_line(size = 0.8) +
  scale_y_continuous(name = "Sequencing depth per sample") +
  scale_colour_manual(values = c("glyph" = "black"),
						name = "Sequencing depth",
						breaks = "glyph",
						labels = c("glyph" = "Metagenomic\nsamples")) +
  labs(x = "Days")	
 ggsave("sequencing_depth.png", width = 12, height = 10)
rps_genus <- as.data.frame(table(droplevels(subset(prokka_all, grepl("30s@ribosomal@protein@s", product2) & treatment == "glyph")$genus)))
rps_family <- as.data.frame(table(droplevels(subset(prokka_all, grepl("30s@ribosomal@protein@s", product2) & treatment == "glyph")$family)))
  
rpl_genus <- as.data.frame(table(droplevels(subset(prokka_all, grepl("50s@ribosomal@protein@l", product2) & treatment == "glyph")$genus)))
rpl_family <- as.data.frame(table(droplevels(subset(prokka_all, grepl("50s@ribosomal@protein@l", product2) & treatment == "glyph")$family)))

write.table(sep = "\t", rps_genus, file = "rps_genus.tsv")
write.table(sep = "\t", rps_family, file = "rps_family.tsv")
write.table(sep = "\t", rpl_family, file = "rpl_family.tsv")
write.table(sep = "\t", rpl_genus, file = "rpl_genus.tsv")

# phn G H I J L M
 table(droplevels(subset(prokka_all, grepl("phosphate", product2) & grepl("phosphon", product2))$product2))

 table(droplevels(subset(prokka_all, grepl("phosphate", product2) & grepl("transport", product2))$product2)
 
#-----------------------------------------------------------------------



#fragen?

#ist richness empfindlicher als abundanz? wie teste ich das im vergleich
#zu "normalen" Genen?

# relative abundance and relative richness are not correlated?? (phnC, phnD, phnH
# more correlation for phnE,
# there is an excel table on that findings "richness_abundance comparison"
# differences in reaction time? because of detection limit per contig?
# back to starting niveau for abundance? or richness? elevated richness over time?

# To do
#plots aller relevanten OTUS

# taxonomy + gene:

# gallaeci besitzt soxBCD gene auf vier contigs und phnCDE 1 contig H noch ein contig
# methylo shikimate dehydrogenase? aroE, 1.1.1.25
#Pseudomonas hat phn[C-N] und alle sox gene
#Hoeflea soxBCD, alle phn außer JK
#taxonomie sollte nicht überbewerttet werden





as.data.frame(table(droplevels(phn_op$genus)))

as.data.frame(table(droplevels(sox_op$genus)))
 as.data.frame(table(droplevels(subset(phn_op, genus == "Pseudomonas")$gene)))


phn <- subset(gene_richness2, grepl("phnN", gene))
#phn <- subset(gene_richness2, grepl("phn[C-N]", gene))

# visualization per product



# cumulated product's visualization  
ggplot(sarc, aes (fill = product2, colour = product2))+
  geom_bar(stat = "identity", aes(x = new_day - 0.75, y = product_rpm), fill = "black", width = 1)+
  geom_bar(stat = "identity", aes(x = new_day + 0.75, y = gene_richness_relative), fill = "blue", width = 1)+
  facet_wrap(~treatment, nrow = 2) 
  
  
###### plots for quick investigation: relative richness containing tax richness
  
ggplot(phn_operon_sox_richness_subset, aes(x = new_day, y = gene_richness2)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(sox_richness_subset, aes(x = new_day, y = gene_richness2_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(housekeeping_richness, aes(x = new_day, y = gene_richness2_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel), size = 1.3)+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(new_genes, aes(x = new_day, y = gene_richness2)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
  
# plot species abundance
Paracoccus <- subset(prokka_all, genus == "Paracoccus")
ggplot(Paracoccus, aes(x = new_day, y = kaiju_rpm, fill = treatment)) +
geom_bar(stat = "identity") +
facet_wrap(~ treatment)


# plot gene abundance with taxonomic annotation
phnJ <- subset(prokka_all, gene == "phnJ")

ggplot(phnJ, aes(x = new_day, y = product_rpm, fill = genus)) +
geom_bar(stat = "identity") +
facet_wrap(~ treatment)


# plot multiple genes or products
phn <- subset(prokka_all, grepl("phn", gene))

ggplot(phn, aes(x = new_day, y = product_rpm, fill = order, colour = gene)) +
#geom_bar(stat = "identity") +
geom_bar(stat = "identity", position = position_dodge(), size = 2) +
facet_wrap(~ treatment)

# use multiple subsets to refine results
sarcos <- subset(prokka_all, grepl("sarcosine", product2))
sarcos <- subset(sarcos, grepl("oxidase", product2))

ggplot(sarcos, aes(x = new_day, y = product_rpm, fill = product2, colour = genus)) +
geom_bar(stat = "identity") +
#geom_bar(stat = "identity", position = position_dodge(), size = 1.5) +
facet_wrap(~ treatment, nrow = 2)

# search for specific gene and specific organism
sarcos_methylo <- subset(sarcos, grepl("Methylo", genus))

ggplot(sarcos_methylo, aes(x = new_day, y = product_rpm, fill = product2, colour = genus)) +
geom_bar(stat = "identity") +
#geom_bar(stat = "identity", position = position_dodge(), size = 1.5) +
facet_wrap(~ treatment, nrow = 2)