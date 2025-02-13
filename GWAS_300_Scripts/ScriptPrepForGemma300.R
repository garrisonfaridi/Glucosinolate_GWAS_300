library(tidyverse)
#PrepForGemma
# Load BLUPs and Accessions
dat_300 = read.csv("300_Inputs/TOUA-glucs_compiledBlUPs_ADG-2024-11-12.tsv", header = T, sep = "\t")
acc_300 = read.csv("2022_Inputs/TOU_accessions.csv", header = T, sep = "\t")

acc_300_cleaned <- acc_300 %>%
  separate(SAMPLE.ID_short.ecotype_id.pop_code, into = c("SAMPLE", "ID", "ecotype_id", "pop_code"), sep = ",") 


# merge with list of genotypes, ordered as in the .fam file
# (from the set of PLINK binary files used as genotypes in GWAS with GEMMA)
fam.file.IDs.300 = read.csv("300_Inputs/GEMMA_Inputs/genotypes_305.fam", sep = " ", header = F)[,1]

dat_300 = merge(acc_300_cleaned, dat_300, by = "ID", all.x = T)

dat_300 = dplyr::left_join(data.frame(ecotype_id= as.character(fam.file.IDs.300)), dat_300, by = "ecotype_id")

# glucosinolates to retain
dat_300 = dat_300[,c("ecotype_id","gsl.7mSh.blup","gsl.8MTO.blup","gsl.3mSOp.blup","gsl.4mSOb.blup",
                     "gsl.5mSOp.blup","gsl.6mSOh.blup","gsl.7mSOh.blup","gsl.8mSOo.blup","gsl.Pren.blup",
                     "gsl.Buen.blup","gsl.Peen.blup","gsl.S2hBuen.blup","gsl.2hPeen.blup","gsl.IM.blup",
                     "gsl.1moIM.blup","gsl.1hIM.blup","gsl.4moIM.blup")]

colnames(dat_300)[1] = "ID"

write.table(dat_300[,c(2:ncol(dat_300))], "dat_tot_ForGemmaGF_300.csv",
            sep = ",", row.names = F, col.names = F, quote = F, eol = "\n")
