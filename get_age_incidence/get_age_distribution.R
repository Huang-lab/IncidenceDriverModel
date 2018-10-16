##### get_age_distribution.R #####
# Kuan-lin Huang 2018
# find non-Cancer pathogenic variant in the ExAC cohort
# targeting specific ethnicity

# # set work dir for testing in dev environ
# bdir = "/Users/khuang/Google Drive/ResearchProjects/IncidenceDriverModel/analysis"
# setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")


##### analyze #####

# colnames(clin)
clin_germ = merge(clin,pathVarP[,c(1,which(!(colnames(pathVarP) %in% colnames(clin))))],by="bcr_patient_barcode",all.x=T)

clin_germ$age_bracket = cut(clin_germ$age_at_initial_pathologic_diagnosis, seq(0,100,5), include.lowest=T, right = F)
clin_germ = clin_germ[!duplicated(clin_germ$bcr_patient_barcode),]
clin_germ$condition = "WT"
clin_germ$condition[clin_germ$HUGO_Symbol=="BRCA1"] = "BRCA1"
clin_germ$condition[clin_germ$HUGO_Symbol=="BRCA2"] = "BRCA2"
clin_germ$condition[clin_germ$HUGO_Symbol=="ATM"] = "ATM"

clin_germ = clin_germ[order(clin_germ$type),]

write.table(clin_germ[,c("type","age_at_initial_pathologic_diagnosis","age_bracket","condition")], quote=F, sep="\t", 
            file = "out/TCGA_all_onset_with_germline_genes.tsv", row.names = F)

clin_germ_table = data.frame(table(clin_germ$type,clin_germ$condition,clin_germ$age_bracket))
colnames(clin_germ_table) = c("cancer","condition","age","count")
clin_germ_table = clin_germ_table[order(clin_germ_table$cancer),]

write.table(clin_germ_table, quote=F, sep="\t", 
            file = "out/TCGA_all_onset_with_germline_genes_summarized.tsv", row.names = F)

