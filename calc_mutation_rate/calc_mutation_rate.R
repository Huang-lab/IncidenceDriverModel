##### calc_mutation_rate.R #####
# Kuan-lin Huang 2018

# # set work dir for testing in dev environ
# bdir = "/Users/khuang/Google Drive/ResearchProjects/IncidenceDriverModel/analysis/calc_mutation_rate"
# setwd(bdir)
source("../global_aes_out.R")
source("../dependency_files.R")


##### analyze #####

colnames(clin)
clin_germ = merge(clin,pathVarP[,c(1,which(!(colnames(pathVarP) %in% colnames(clin))))],by="bcr_patient_barcode",all.x=T)

clin_germ$age_bracket = cut(clin_germ$age_at_initial_pathologic_diagnosis, seq(0,100,5), include.lowest=T, right = F)
clin_germ = clin_germ[!duplicated(clin_germ$bcr_patient_barcode),]
clin_germ$condition = "WT"
clin_germ$condition[clin_germ$HUGO_Symbol=="BRCA1"] = "BRCA1"
clin_germ$condition[clin_germ$HUGO_Symbol=="BRCA2"] = "BRCA2"
clin_germ$condition[clin_germ$HUGO_Symbol=="ATM"] = "ATM"

clin_germ = clin_germ[order(clin_germ$type),]
clin_germ_table = data.frame(table(clin_germ$type,clin_germ$condition,clin_germ$age_bracket))
colnames(clin_germ_table) = c("cancer","condition","age","count")
clin_germ_table = clin_germ_table[order(clin_germ_table$cancer),]
clin_germ$cancer_germ = paste(clin_germ$type,clin_germ$HUGO_Symbol,sep="_")
clin_germ$cancer_germ[clin_germ$cancer_germ %in% names(table(clin_germ$cancer_germ)[table(clin_germ$cancer_germ)<5])] = NA

# by cancer type
somatic_count_f = "~/Box Sync/Huang_lab/manuscripts/germlineSomatic/analysis/mutation_rate/out/TCGA_sample_mutation_driver_count.tsv"
somatic_count = read.table(header=T, quote = "", sep="\t", file = somatic_count_f, stringsAsFactors=FALSE)
somatic_count_cancer = merge(clin[,c("bcr_patient_barcode","type")],somatic_count,by="bcr_patient_barcode") 

cat("conduct somatic mutation count analysis in samples with both germline and somatic information:",nrow(somatic_count_cancer),"samples.\n")
somatic_count_by_cancer = aggregate(. ~ type, somatic_count_cancer[,-1], function(x) c(median = median(x), sd = sd(x)))
write.table(somatic_count_by_cancer, quote=F, sep="\t", 
            file = "out/TCGA_mutation_rate_by_cancer_type.tsv", row.names = F)

# germline
somatic_count_germline = merge(clin_germ[,c("bcr_patient_barcode","cancer_germ","age_at_initial_pathologic_diagnosis")],somatic_count_cancer,by="bcr_patient_barcode") 
somatic_count_germline = somatic_count_germline[!is.na(somatic_count_germline$cancer_germ),]
somatic_count_germline$somatic_mutation_count_per_year = somatic_count_germline$somatic_mutation_count/somatic_count_germline$age_at_initial_pathologic_diagnosis
  
#somatic_count_by_germline_cancer = aggregate(. ~ cancer_germ, somatic_count_germline[,-c(1,3,4)], function(x) c(median = median(x), sd = sd(x)))
somatic_count_by_germline_cancer = aggregate(. ~ cancer_germ, somatic_count_germline[,-c(1,3,4)], function(x) c(median = median(x), sd = sd(x)))
write.table(somatic_count_by_germline_cancer, quote=F, sep="\t", 
            file = "out/TCGA_mutation_rate_by_cancer_type_germline_status.tsv", row.names = F)

# subtype analysis
tn = paste("../../data/SubtypeAssignments.xlsx")
subtype = data.frame(readxl::read_xlsx(tn))
colnames(subtype) = c("bcr_patient_barcode","SAMPLE_BARCODE","type","subtype")
subtype$subtype[!is.na(subtype$subtype) & subtype$subtype=="Not_Applicable"] = NA

somatic_count_cancer_subtype = merge(subtype, somatic_count_cancer, by=c("bcr_patient_barcode","type"))
cat("conduct somatic mutation count analysis in samples with subtype information:",nrow(somatic_count_cancer_subtype),"samples.\n")
somatic_count_cancer_subtype$cancer_subtype = paste(somatic_count_cancer_subtype$type,somatic_count_cancer_subtype$subtype,sep="_")
somatic_count_by_subtype_cancer = aggregate(. ~ cancer_subtype, somatic_count_cancer_subtype[,-c(1:4)], function(x) c(median = median(x), sd = sd(x)))
write.table(somatic_count_by_subtype_cancer, quote=F, sep="\t", 
            file = "out/TCGA_mutation_rate_by_cancer_type_subtype.tsv", row.names = F)

### plotting ###
somatic_count_cancer_m = melt(somatic_count_cancer,id.vars = c("bcr_patient_barcode", "type" ))
p = ggplot(data=somatic_count_cancer_m,aes(x = type, y=log2(value)))
p = p + facet_grid(variable~., scale="free_y")
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
#p = p + geom_dotplot(aes(fill= type),binaxis="y", stackdir = "center",colour=NA,binwidth=100,dotsize=0.1)
p = p + geom_violin(aes(fill = type), alpha=0.3,color=NA)
p = p + geom_jitter(aes(color= type),height=0, alpha = 0.1, size=1,stroke=0) 
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "log2(Count of mutations)")
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "none")
p = p + expand_limits(y=0) 
p
fn = 'out/mutation_driver_rate_by_cancer_type.pdf'
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)

#
somatic_count_germline_m = melt(somatic_count_germline[,-c(3,4)],id.vars = c("bcr_patient_barcode", "cancer_germ" ))
p = ggplot(data=somatic_count_germline_m,aes(x =cancer_germ, y=log2(value)))
p = p + facet_grid(variable~., scale="free_y")
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
#p = p + geom_dotplot(aes(fill= type),binaxis="y", stackdir = "center",colour=NA,binwidth=100,dotsize=0.1)
p = p + geom_violin(aes(fill = cancer_germ), alpha=0.3,color=NA)
p = p + geom_jitter(aes(color= cancer_germ),height=0, alpha = 0.1, size=1,stroke=0) 
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "log2(count of mutations)") 
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "none")
p = p + expand_limits(y=0) 
p
fn = 'out/mutation_driver_rate_by_cancer_type_germline.pdf'
ggsave(file=fn, w = 10, useDingbats=FALSE,limitsize=FALSE)

#
somatic_count_cancer_subtype_m = melt(somatic_count_cancer_subtype[,-c(2:4)],id.vars = c("bcr_patient_barcode", "cancer_subtype" ))
p = ggplot(data=somatic_count_cancer_subtype_m,aes(x =cancer_subtype, y=log2(value)))
p = p + facet_grid(variable~., scale="free_y")
p = p + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                     geom = "crossbar", width = 0.8)
#p = p + geom_dotplot(aes(fill= type),binaxis="y", stackdir = "center",colour=NA,binwidth=100,dotsize=0.1)
p = p + geom_violin(aes(fill = cancer_subtype), alpha=0.3,color=NA)
p = p + geom_jitter(aes(color= cancer_subtype),height=0, alpha = 0.1, size=1,stroke=0) 
p = p  + theme_bw() #+ theme(legend.position="none")
p = p + labs(x = "Cancer", y = "log2(count of mutations)") 
p = p + theme(axis.title = element_text(size=18), axis.text.y = element_text(size=14), axis.text.x = element_text(colour="black", size=14, angle = 90, vjust = 0.5))
p = p + theme(legend.position = "none")
p = p + expand_limits(y=0) 
p
fn = 'out/mutation_driver_rate_by_cancer_type_subtype.pdf'
ggsave(file=fn, w = 12, useDingbats=FALSE,limitsize=FALSE)
