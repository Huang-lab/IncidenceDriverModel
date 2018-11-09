##### plot_obs....R #####
# Kuan-lin Huang 2018

# # set work dir for testing in dev environ
bdir = "/Users/khuang/Box Sync/Huang_lab/manuscripts/IncidenceDriverModel/analysis/compare_rates"
setwd(bdir)
source("../global_aes_out.R")

subtype_n = paste("../../data/SubtypeAssignments.xlsx")
subtype = data.frame(readxl::read_xlsx(subtype_n))

clin_n = paste("../../data/PanCan_ClinicalData_V4_wAIM.txt")
clin = read.table(header=T, sep="\t", file=clin_n, fill=T)

tn = paste("../calc_mutation_rate/out/TCGA_mutation_rate_by_cancer_type.xlsx")
tt = data.frame(readxl::read_xlsx(tn))

# erlang_n = paste("../../data/Belikov2017_Erlang_estimates.xlsx")
# erlang = data.frame(readxl::read_xlsx(erlang_n))
# 
# tcga_n = paste("../../data/TCGA_mutation_rate_by_cancer_type.tsv")
# tcga = read.table(header=TRUE, sep="\t", file=tcga_n, fill=T)

dNdS_n = paste("../../data/dNdS_Fig4_exome_wide_number_drivers.xlsx")
dNdS = data.frame(readxl::read_xlsx(dNdS_n))
colnames(dNdS)[2:4] = paste("dNdS",colnames(dNdS)[2:4],sep="_")

select_n = paste("../../data/pancan23_select_driver_count_by_sample.tsv")
select = read.table(header=F, sep="\t", file=select_n, fill=T)
colnames(select) = c("sampleID","select_driver")
select$bcr_patient_barcode = gsub("-01","",gsub("\\.","-",select$sampleID))
select_clin = merge(select,clin,by="bcr_patient_barcode")

select_clin_aggregated = data.frame(aggregate(. ~ type, select_clin[,c("type","select_driver")], function(x) c(median = median(x), sd = sd(x))))
select_clin_aggregated = data.frame(aggregate(. ~ type, select_clin[,c("type","select_driver")], function(x) c(median = median(x))))
#colnames(select_clin_aggregated) = c("select_driver","select_driver_sd")
# perhaps should line up sample to select subset

tt_merge = merge(tt,dNdS,by="type")
tt_all_merge = merge(tt_merge,select_clin_aggregated,by="type")
colnames(tt_all_merge) = gsub("\\.","_",colnames(tt_all_merge))

tt_all_merge_medians = tt_all_merge[,c(grep("median",colnames(tt_all_merge)),which(colnames(tt_all_merge) %in% c("type","Erlang_estimation_driver","dNdS_MLE","select_driver")))]
tt_all_merge_medians = tt_all_merge_medians[,-c(grep("per_year",colnames(tt_all_merge_medians)))]

# change column names for plotting
colnames(tt_all_merge_medians) = gsub("_median","",colnames(tt_all_merge_medians))
colnames(tt_all_merge_medians) = gsub("driver_count","panSoftware_driver",colnames(tt_all_merge_medians))
colnames(tt_all_merge_medians)[colnames(tt_all_merge_medians)=="age_at_initial_pathologic_diagnosis"] = "onset_age"
colnames(tt_all_merge_medians)[colnames(tt_all_merge_medians)=="panSoftware_driver"] = "Bailey_panSoftware"
colnames(tt_all_merge_medians)[colnames(tt_all_merge_medians)=="loose_panSoftware_driver"] = "Bailey_panSoftware_loose"
colnames(tt_all_merge_medians)[colnames(tt_all_merge_medians)=="dNdS_MLE"] = "Martincorena_dNdS"
colnames(tt_all_merge_medians)[colnames(tt_all_merge_medians)=="select_driver"] = "Ciriello_SELECT"

library("GGally")
ggpairs(data=tt_all_merge_medians[,-c(which(colnames(tt_all_merge_medians) =="type"))], #columns = 1:ncol(data), 
        axisLabels = "show")  + theme_bw()
fn = 'out/all_theoretical_vs_observed_driver_pairwise_comparisons.pdf'
ggsave(fn,useDingbat=F, h = 10, w = 10)

tt_all_merge_medians_m = melt(tt_all_merge_medians,id.vars = "type")
tt_all_merge_medians_m$type = factor(tt_all_merge_medians_m$type, levels=unique(tt_all_merge_medians_m$type)[order(tt_all_merge_medians_m[tt_all_merge_medians_m$variable=="somatic_mutation_count",]$value)])

p = ggplot(data=tt_all_merge_medians_m,aes(x=type,y=value, color = type))
p = p + facet_grid(variable~.,scale="free")
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=somatic_mutation_count.median-somatic_mutation_count.sd, ymax=somatic_mutation_count.median+somatic_mutation_count.sd), width=.2)
#,position=position_dodge(0.05))
#p = p + geom_label_repel(aes(label=type), alpha=0.8)
p = p  + theme_bw() + theme_nogrid() + 
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p = p + getPCACancerColor()
p + expand_limits(y=0)#+ labs(x="Signature",y = "Cancer") 

fn = 'out/all_theoretical_vs_observed_driver_points.pdf'
ggsave(fn,h=8,useDingbat=F)

### previous plotting ###

p = ggplot(data=tt,aes(x=Erlang_estimation,y=somatic_mutation_count_per_year.median, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=somatic_mutation_count.median-somatic_mutation_count.sd, ymax=somatic_mutation_count.median+somatic_mutation_count.sd), width=.2)
                      #,position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 

fn = 'out/somatic_mutation_rate_vs_erlang_somatic_prediction.pdf'
ggsave(fn,useDingbat=F)

p = ggplot(data=tt,aes(x=Erlang_estimation_driver,y=loose_driver_count_per_year.median, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=loose_driver_count_per_year.median-loose_driver_count_per_year.sd, ymax=loose_driver_count_per_year.loose_driver_count_per_year.sd), width=.2),position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8) + geom_abline(slope=1, alpha=0.2)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 

fn = 'out/loose_driver_rate_vs_erlang_prediction.pdf'
ggsave(fn,useDingbat=F)


p = ggplot(data=tt,aes(x=Erlang_estimation_driver,y=loose_driver_count.median, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=loose_driver_count_per_year.median-loose_driver_count_per_year.sd, ymax=loose_driver_count_per_year.loose_driver_count_per_year.sd), width=.2),position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8) + geom_abline(slope=1, alpha=0.2)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 
fn = 'out/loose_driver_vs_erlang_prediction.pdf'
ggsave(fn,useDingbat=F)

p = ggplot(data=tt,aes(x=Erlang_estimation_driver,y=somatic_mutation_count_per_year.median, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=loose_driver_count_per_year.median-loose_driver_count_per_year.sd, ymax=loose_driver_count_per_year.loose_driver_count_per_year.sd), width=.2),position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 
fn = 'out/somatic_mutation_rate_vs_erlang_prediction.pdf'
ggsave(fn,useDingbat=F)

p = ggplot(data=tt,aes(x=Erlang_estimation_driver,y=somatic_mutation_count_per_year.median, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=loose_driver_count_per_year.median-loose_driver_count_per_year.sd, ymax=loose_driver_count_per_year.loose_driver_count_per_year.sd), width=.2),position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 
fn = 'out/somatic_mutation_rate_vs_erlang_prediction.pdf'
ggsave(fn,useDingbat=F)

p = ggplot(data=tt,aes(x=Erlang_estimation_driver,y=somatic_mutation_count.median, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=loose_driver_count_per_year.median-loose_driver_count_per_year.sd, ymax=loose_driver_count_per_year.loose_driver_count_per_year.sd), width=.2),position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 
fn = 'out/somatic_mutation_vs_erlang_prediction.pdf'
ggsave(fn,useDingbat=F)

p = ggplot(data=tt,aes(x=Erlang_estimation_driver,y=dNdS_estimation_driver, color = class))
p = p + geom_point()
#p = p + geom_errorbar(aes(ymin=loose_driver_count_per_year.median-loose_driver_count_per_year.sd, ymax=loose_driver_count_per_year.loose_driver_count_per_year.sd), width=.2),position=position_dodge(0.05))
p = p + geom_label_repel(aes(label=type), alpha=0.8) + geom_abline(slope=1, alpha=0.2)
p = p  + theme_bw() + theme_nogrid() +
  theme(axis.text.x = element_text(colour="black", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(colour="black", size=12),axis.ticks = element_blank())#element_text(colour="black", size=14))
p + expand_limits(y=0,x=0)#+ labs(x="Signature",y = "Cancer") 
fn = 'out/dnds_vs_erlang_prediction.pdf'
ggsave(fn,useDingbat=F)
