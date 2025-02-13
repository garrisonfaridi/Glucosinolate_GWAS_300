# Purpose: Tallies up the number of significant SNPs and the top p-value by candidate gene
#          for each trait (column) in a the merged p-value table created by script_combineGwasWideFormat.R.
#          The output includes tables of the number of significant SNPs and top p-values per
#          glucosinolate candidate gene, and also a table summarizing per study (Brachi, Katz, WU, TOU-A)
#          rather than per trait.

# ALIPHATIC
# specify maf cutoff for the analysis
maf.cutoff = 0.05

##### DEFINE SNP LIST, CANDIDATE GENES #####

# load SNPs, maf > 0.05, TOU-A
snps.toua        = subset(read.delim("output_300/300_gluc1_maf0.03_miss0.1GF.assoc.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)

# get full list of SNPs
snps.all         = snps.toua[,c("chr","ps")]

# rs from vcf sometimes has chr first (tou), sometimes ps first (regmap/1001 imputed), so standardize
snps.all$rs      = paste(snps.all$chr, snps.all$ps, sep = "_")
snps.all         = unique(snps.all)

# # save SNPs list so it can be loaded directly later
# write.csv(snps.all, "full_SNP_list_maf05_wTouFrachonMiss10.csv")
# # to load:
# snps.all = read.csv("full_SNP_list_maf05_wTouFrachonMiss10.csv")

# load candidate gene list
candidate.genes  = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt", sep = "\t")
candidate.genes  = subset(candidate.genes, aliphatic == "yes" | indolic == "yes")
nrow(candidate.genes) # 45 candidate genes

# load table of gene models (left and right boundaries)
gene.models      = read.delim("TAIR/TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt", sep = "\t")

# merge gene model info with candidate gene list
candidate.genes  = merge(candidate.genes, gene.models, by = "gene_id")
nrow(candidate.genes) # 45 genes still retained with gene model coordinates now included

##### ALIPHATICS - prep #####

# Optional: subset candidate genes to the significant loci for aliphatic GSL traits
# GS-OX1, GS-OH, BCAT3, AOP2, MAM1

#candidate.genes.ali = subset(candidate.genes, gene_id %in% c("AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT5G23010"))

candidate.genes.ali = subset(candidate.genes, aliphatic == "yes")

candidate.snps.per.gene = list()

# get list of candidate SNPs per gene
for (i in 1:nrow(candidate.genes.ali)){
  
  snps = subset(snps.all, chr == candidate.genes.ali[i,"chr"] &
                  ps > candidate.genes.ali[i,"ps_left"] &
                  ps < candidate.genes.ali[i,"ps_right"],
                select = rs
  )
  
  candidate.snps.per.gene[[ as.character(candidate.genes.ali[i,"gene_id"]) ]] = snps
  
}

##### ALIPHATICS - count sig SNPs #####

# specify which sets of GWAS to parse for top p-values and number of significant SNPs
blups.tou.m10  = subset(read.csv("output_300/snp_table_all_aliphatic_pvals.csv"), af > maf.cutoff)

print(dim(blups.tou.m10))  # Check dimensions of blups.tou.m10
gwas.to.parse = list("blups.tou.m10" = blups.tou.m10)  # Only one dataset

names(gwas.to.parse) = c("blups.tou.m10")

# create data frame to hold top significant SNP pval
top.pvals.sig.best = data.frame(matrix(nrow = nrow(candidate.genes.ali), ncol = length(gwas.to.parse)+1))
colnames(top.pvals.sig.best) = c("gene_id", names(gwas.to.parse))
top.pvals.sig.best$gene_id = candidate.genes.ali$gene_id

print("Script is running up to this point")



for (i in 1:length(names(gwas.to.parse))){
  
  # print progress update
  print(paste0("Parsing: ", names(gwas.to.parse)[[i]], " ..."))
  
  # specify focal gwas
  focal.gwas = gwas.to.parse[[i]]  # Use direct indexing
  
  print(dim(focal.gwas))
  # create data frame to hold top SNP p-values
  # dimensions: rows = number of loci, columns = number of traits 
  top.pvals = data.frame(matrix(nrow = nrow(candidate.genes.ali), ncol = ncol(focal.gwas) - 3))
  colnames(top.pvals) = c("gene_id", colnames(focal.gwas[5:ncol(focal.gwas)]))
  top.pvals$gene_id = candidate.genes.ali$gene_id
  
  # create identical data frame, but for only considering significant SNPs
  top.pvals.sig = top.pvals
  
  # create identical data frame to hold number of sig. SNPs
  num.sig.snps = top.pvals
  
  for (j in top.pvals$gene_id) {
    
    #for (k in 4:ncol(focal.gwas)){
    for (k in colnames(focal.gwas)[5:ncol(focal.gwas)]){
      
      # subset to candidate SNPs for the trait
      sig.snps = focal.gwas[focal.gwas[,k] < maf.cutoff / nrow(focal.gwas) / (ncol(focal.gwas)-5)
                            & focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      candidate.snps = focal.gwas[focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      
      # top p value
      top.pvals[top.pvals$gene_id == j, k] = min(candidate.snps)
      
      # top p value (sig. SNPs only)
      top.pvals.sig[top.pvals$gene_id == j, k] = min(sig.snps)
      
      
      # number of significant SNPs
      num.sig.snps[num.sig.snps$gene_id == j, k] = length(sig.snps)
      
    }
    
  }
  
  top.pvals.sig.best[,names(gwas.to.parse)[i]] = top.pvals.sig[,"p_min"]
  
  # individually name and save these as desired  
  write.csv(top.pvals, paste0("output_300/", "aliphatic.top.pvals.300", names(gwas.to.parse)[i], ".csv"))
  write.csv(top.pvals.sig, paste0("output_300/", "aliphatic.top.pvals.sig.300", names(gwas.to.parse)[i], ".csv"))
  write.csv(num.sig.snps, paste0("output_300/", "aliphatic.num.sig.snps.300", names(gwas.to.parse)[i], ".csv"))
  
  # print data frame of top SNP p-values to file
  
}

write.csv(top.pvals.sig.best, "top_pvals_sig_best_aliphatic_allMiss10_300.csv")

# Purpose: Tallies up the number of significant SNPs and the top p-value by candidate gene
#          for each trait (column) in a the merged p-value table created by script_combineGwasWideFormat.R.
#          The output includes tables of the number of significant SNPs and top p-values per
#          glucosinolate candidate gene, and also a table summarizing per study (Brachi, Katz, WU, TOU-A)
#          rather than per trait.

# INDOLIC
# specify maf cutoff for the analysis
maf.cutoff = 0.05

##### DEFINE SNP LIST, CANDIDATE GENES #####

# load SNPs, maf > 0.05, TOU-A
snps.toua        = subset(read.delim("output_300/300_gluc1_maf0.03_miss0.1GF.assoc.txt", sep="\t")[,c("chr","ps","af")], af > 0.05)

# get full list of SNPs
snps.all         = snps.toua[,c("chr","ps")]

# rs from vcf sometimes has chr first (tou), sometimes ps first (regmap/1001 imputed), so standardize
snps.all$rs      = paste(snps.all$chr, snps.all$ps, sep = "_")
snps.all         = unique(snps.all)

# # save SNPs list so it can be loaded directly later
# write.csv(snps.all, "full_SNP_list_maf05_wTouFrachonMiss10.csv")
# # to load:
# snps.all = read.csv("full_SNP_list_maf05_wTouFrachonMiss10.csv")

# load candidate gene list
candidate.genes  = read.delim("pathway_BiosyntheticGenes_HarunTable2.txt", sep = "\t")
candidate.genes  = subset(candidate.genes, aliphatic == "yes" | indolic == "yes")
nrow(candidate.genes) # 45 candidate genes

# load table of gene models (left and right boundaries)
gene.models      = read.delim("TAIR/TAIR10_GFF3_geneBoundariesTableADG_extendedWindows.txt", sep = "\t")

# merge gene model info with candidate gene list
candidate.genes  = merge(candidate.genes, gene.models, by = "gene_id")
nrow(candidate.genes) # 45 genes still retained with gene model coordinates now included

##### ALIPHATICS - prep #####

# Optional: subset candidate genes to the significant loci for aliphatic GSL traits
# GS-OX1, GS-OH, BCAT3, AOP2, MAM1

#candidate.genes.ali = subset(candidate.genes, gene_id %in% c("AT1G65860","AT2G25450","AT3G49680","AT4G03060","AT5G23010"))

candidate.genes.ind = subset(candidate.genes, indolic == "yes")

candidate.snps.per.gene = list()

# get list of candidate SNPs per gene
for (i in 1:nrow(candidate.genes.ind)){
  
  snps = subset(snps.all, chr == candidate.genes.ind[i,"chr"] &
                  ps > candidate.genes.ind[i,"ps_left"] &
                  ps < candidate.genes.ind[i,"ps_right"],
                select = rs
  )
  
  candidate.snps.per.gene[[ as.character(candidate.genes.ind[i,"gene_id"]) ]] = snps
  
}

##### ALIPHATICS - count sig SNPs #####

# specify which sets of GWAS to parse for top p-values and number of significant SNPs
blups.tou.m10  = subset(read.csv("output_300/snp_table_all_indolic_pvals.csv"), af > maf.cutoff)

print(dim(blups.tou.m10))  # Check dimensions of blups.tou.m10
gwas.to.parse = list("blups.tou.m10" = blups.tou.m10)  # Only one dataset

names(gwas.to.parse) = c("blups.tou.m10")

# create data frame to hold top significant SNP pval
top.pvals.sig.best = data.frame(matrix(nrow = nrow(candidate.genes.ind), ncol = length(gwas.to.parse)+1))
colnames(top.pvals.sig.best) = c("gene_id", names(gwas.to.parse))
top.pvals.sig.best$gene_id = candidate.genes.ind$gene_id

print("Script is running up to this point")



for (i in 1:length(names(gwas.to.parse))){
  
  # print progress update
  print(paste0("Parsing: ", names(gwas.to.parse)[[i]], " ..."))
  
  # specify focal gwas
  focal.gwas = gwas.to.parse[[i]]  # Use direct indexing
  
  print(dim(focal.gwas))
  # create data frame to hold top SNP p-values
  # dimensions: rows = number of loci, columns = number of traits 
  top.pvals = data.frame(matrix(nrow = nrow(candidate.genes.ind), ncol = ncol(focal.gwas) - 3))
  colnames(top.pvals) = c("gene_id", colnames(focal.gwas[5:ncol(focal.gwas)]))
  top.pvals$gene_id = candidate.genes.ind$gene_id
  
  # create identical data frame, but for only considering significant SNPs
  top.pvals.sig = top.pvals
  
  # create identical data frame to hold number of sig. SNPs
  num.sig.snps = top.pvals
  
  for (j in top.pvals$gene_id) {
    
    #for (k in 4:ncol(focal.gwas)){
    for (k in colnames(focal.gwas)[5:ncol(focal.gwas)]){
      
      # subset to candidate SNPs for the trait
      sig.snps = focal.gwas[focal.gwas[,k] < maf.cutoff / nrow(focal.gwas) / (ncol(focal.gwas)-5)
                            & focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      candidate.snps = focal.gwas[focal.gwas$rs %in% candidate.snps.per.gene[[j]]$rs, k]
      
      # top p value
      top.pvals[top.pvals$gene_id == j, k] = min(candidate.snps)
      
      # top p value (sig. SNPs only)
      top.pvals.sig[top.pvals$gene_id == j, k] = min(sig.snps)
      
      
      # number of significant SNPs
      num.sig.snps[num.sig.snps$gene_id == j, k] = length(sig.snps)
      
    }
    
  }
  
  top.pvals.sig.best[,names(gwas.to.parse)[i]] = top.pvals.sig[,"p_min"]
  
  # individually name and save these as desired  
  write.csv(top.pvals, paste0("output_300/", "indolic.top.pvals.300", names(gwas.to.parse)[i], ".csv"))
  write.csv(top.pvals.sig, paste0("output_300/", "indolic.top.pvals.sig.300", names(gwas.to.parse)[i], ".csv"))
  write.csv(num.sig.snps, paste0("output_300/", "indolic.num.sig.snps.300", names(gwas.to.parse)[i], ".csv"))
  
  # print data frame of top SNP p-values to file
  
}

write.csv(top.pvals.sig.best, "output_300/top_pvals_sig_best_indolic_allMiss10_300.csv")

