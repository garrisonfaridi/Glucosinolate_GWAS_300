# Manhattan Plotting

# Load candidate SNP files
candidate.snps = list()
candidate.snps[["aliphatic"]] = read.table("candidate_snps_aliphatic_30kb_majorLociExtended_maf05_wTouFrachonGenoMiss10.txt", h = F) 

# Load table with GWAS p-values combined across all traits for Tou-A 
tou.ali.blup.bestp    = read.csv("output_300/snp_table_all_aliphatic_pvals.csv")
#tou.ali.ratio.bestp   = read.csv("/ratio_aliphatic_assoc/snp_table_all_pvals.csv")

# Run function above to make Manhattan Plots
# (user: manually modify gwas.all and gwas.name to use names from directly above)
gwas.all  = tou.ali.blup.bestp
gwas.name = "tou.ali.bestp300"
gsl.type  = "aliphatic"
manhattan_single(rs               = gwas.all$rs, # position of each SNP
                 chr              = gwas.all$chr, # chromosome of each SNP
                 ps               = gwas.all$ps, # unique ID of each SNP
                 pval             = gwas.all$p_min, # p value for each SNP
                 max.y            = 25, # specify max -log10(p) value
                 outfile          = paste0("output_300/",gwas.name,".png"), # specify location for 
                 snps2highlight   = subset(gwas.all, rs %in% candidate.snps[[gsl.type]]$V1)[,"rs"], # specify candidate SNPs
                 sig.thresh       = .05 / nrow(gwas.all), # single GWAS significance threshold
                 sig.thresh.multi = .05 / nrow(gwas.all)/ (ncol(gwas.all) - 4 - 1), # multi-GWAS combined significance threshold
                 point.scale      = 1.3, # size of points
                 gap.size         = 2e6 # size of gap (number of positions) between each chromosome on x-axis of plot
)

# Indolic GWAS

# Load candidate SNP files
candidate.snps = list()
candidate.snps[["indolic"]]   = read.table("candidate_snps_indolic_30kb_majorLociExtended_maf05_wTouFrachonGenoMiss10.txt", h = F)

# Load table with GWAS p-values combined across all traits for Tou-A 
tou.ind.blup.bestp    = read.csv("output_300/snp_table_all_indolic_pvals.csv")
#tou.ind.ratio.bestp   = read.csv("/ratio_indolic_assoc/snp_table_all_pvals.csv")
#tou.mvlmm             = subset(read.delim("mvlmmA_maf03.assoc.fixed.txt"), af > 0.05)

# Run function above to make Manhattan Plots
# (user: manually modify gwas.all and gwas.name to use names from directly above)
gwas.all  = tou.ind.blup.bestp
gwas.name = "tou.ind.bestp300"
gsl.type  = "indolic"
manhattan_single(rs               = gwas.all$rs,
                 chr              = gwas.all$chr,
                 ps               = gwas.all$ps,
                 pval             = gwas.all$p_min,
                 max.y            = 15,
                 outfile          = paste0("output_300/",gwas.name,".png"),
                 snps2highlight   = subset(gwas.all, rs %in% candidate.snps[[gsl.type]]$V1)[,"rs"],
                 sig.thresh       = .05 / nrow(gwas.all),
                 sig.thresh.multi = .05 / nrow(gwas.all)/ (ncol(gwas.all) - 4 - 1),
                 point.scale      = 1.3,
                 gap.size         = 2e6
)

#mvlmm indolic manhattan
#mvlmm Manhattan

# Load candidate SNP files
candidate.snps = list()
candidate.snps[["indolic"]]   = read.table("candidate_snps_indolic_30kb_majorLociExtended_maf05_wTouFrachonGenoMiss10.txt", h = F)

# Load table with GWAS p-values combined across all traits for Tou-A 
#tou.ind.mvlmm   = read.csv("output_300/mvlmm_maf03_300_GF.assoc.txt", sep = "\t")
#tou.ind.ratio.bestp   = read.csv("/ratio_indolic_assoc/snp_table_all_pvals.csv")
tou.mvlmm             = subset(read.delim("output_300/mvlmm_maf03_300_GF.assoc.txt"), af > 0.05)
tou.mvlmm$rs <- gsub(":", "_", tou.mvlmm$rs)

# Run function above to make Manhattan Plots
# (user: manually modify gwas.all and gwas.name to use names from directly above)
gwas.all  = tou.mvlmm
gwas.name = "tou.ind.300.mvlmm"
gsl.type  = "indolic"
manhattan_single(rs               = gwas.all$rs,
                 chr              = gwas.all$chr,
                 ps               = gwas.all$ps,
                 pval             = gwas.all$p_wald,
                 max.y            = 15,
                 outfile          = paste0("output_300/",gwas.name,".png"),
                 snps2highlight   = subset(gwas.all, rs %in% candidate.snps[[gsl.type]]$V1)[,"rs"],
                 sig.thresh       = .05 / nrow(gwas.all),
                 sig.thresh.multi = .05 / nrow(gwas.all),
                 point.scale      = 1.3,
                 gap.size         = 2e6
)