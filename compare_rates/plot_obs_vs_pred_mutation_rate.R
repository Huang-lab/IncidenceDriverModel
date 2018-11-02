##### plot_obs....R #####
# Kuan-lin Huang 2018

# # set work dir for testing in dev environ
# bdir = "/Users/khuang/Google Drive/ResearchProjects/IncidenceDriverModel/analysis/calc_mutation_rate"
# setwd(bdir)
source("../global_aes_out.R")


tn = paste("out/TCGA_mutation_rate_by_cancer_type.xlsx")
tt = data.frame(readxl::read_xlsx(tn))

### plotting ###

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
