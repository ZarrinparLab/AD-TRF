setwd("~/scratch/AD_EM/metabolomics/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggvenn)
library(VennDiagram) 
library(grDevices)
library(gplots)
library(ggpubr)
library(BioVenn) 
#################################################
#clean log table

mtb<-fread( "data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_log.csv")%>%
  dplyr::rename(FeatureID=row.ID)

md<-fread("data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces")%>%
  mutate(ATTRIBUTE_sex=ifelse(ATTRIBUTE_sex=="F","female",ATTRIBUTE_sex),
         genofedZT=paste(ATTRIBUTE_ZT_collection_time, ATTRIBUTE_geno_n_feedingcond, sep = "_"))%>%
  mutate(ATTRIBUTE_timepoint = as.numeric(sub("^ZT", "", ATTRIBUTE_ZT_collection_time)))

#################################################

#just tg-ALF
mdTg<-md%>%filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond=="Tg-ALF")%>%
  mutate(ATTRIBUTE_ZT_collection_time=factor(ATTRIBUTE_ZT_collection_time,
                                             levels=c("ZT0","ZT4","ZT8","ZT12","ZT16","ZT20")))%>%
  arrange(ATTRIBUTE_ZT_collection_time)

mtbTg<-mtb%>%dplyr::select(FeatureID, all_of(mdTg$sampleid))
write.table(mtbTg,"data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_fecal_tgALF_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just ntg-ALF
mdNTg<-md%>%filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond=="NTg-ALF")%>%
  mutate(ATTRIBUTE_ZT_collection_time=factor(ATTRIBUTE_ZT_collection_time,
                                             levels=c("ZT0","ZT4","ZT8","ZT12","ZT16","ZT20")))%>%
  arrange(ATTRIBUTE_ZT_collection_time)

mtbNTg<-mtb%>%dplyr::select(FeatureID, all_of(mdNTg$sampleid))
write.table(mtbNTg,"data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_fecal_ntgALF_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just tg-TRF
mdTRF<-md%>%filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond=="Tg-TRF")%>%
  mutate(ATTRIBUTE_ZT_collection_time=factor(ATTRIBUTE_ZT_collection_time,
                                             levels=c("ZT0","ZT4","ZT8","ZT12","ZT16","ZT20")))%>%
  arrange(ATTRIBUTE_ZT_collection_time)

mtbTRF<-mtb%>%dplyr::select(FeatureID, all_of(mdTRF$sampleid))
write.table(mtbTRF,"data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_fecal_tgTRF_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#################################################

#run metacycle

#on tg-alf
meta2d(infile="data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_fecal_tgALF_log.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="metacycle/TgALF_metacycle",timepoints=rep(seq(0, 20, by=4),times=c(9,9,5,9,9,9)),minper=20,maxper=24)

#on ntg-alf
meta2d(infile="data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_fecal_ntgALF_log.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="metacycle/NTgALF_metacycle",timepoints=rep(seq(0, 20, by=4),times=c(8,7,7,8,8,8)),minper=20,maxper=24)

#on tg-trf
meta2d(infile="data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_fecal_tgTRF_log.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="metacycle/TgTRF_metacycle",timepoints=rep(seq(0, 20, by=4),times=c(8,8,2,5,8,8)),minper=20,maxper=24)

#################################################

annotations<-fread("data_files/AD_GNPS_BA/nf_output/library/library_matches_cleanwmz.csv")%>%
  mutate(Compound_Name = ifelse(grepl("^NA(?:\\|NA)+$", Compound_Name), NA, Compound_Name),
         FeatureID=as.character(FeatureID))

#read in Tg-ALF metacycle hits
tg_metacyc<-fread("metacycle/TgALF_metacycle/meta2d_gnps_quant_fecal_tgALF_log.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,annotations, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=ifelse(is.na(Compound_Name),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))%>%
  mutate(conditions="Tg-ALF")

#read in NTg-ALF metacycle hits
ntg_metacyc<-fread("metacycle/NTgALF_metacycle/meta2d_gnps_quant_fecal_ntgALF_log.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,annotations, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=ifelse(is.na(Compound_Name),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))%>%
  mutate(conditions="NTg-ALF")

#read in Tg-TRF metacycle hits
trf_metacyc<-fread("metacycle/TgTRF_metacycle/meta2d_gnps_quant_fecal_tgTRF_log.txt")%>%
  dplyr::rename(FeatureID=CycID) %>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,annotations, by="FeatureID")%>%
  arrange(JTK_adjphase)%>%
  mutate(label_name=ifelse(is.na(Compound_Name),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))%>%
  mutate(conditions="Tg-TRF")

#made summ table
all_res<-rbind(tg_metacyc,ntg_metacyc,trf_metacyc)
write.table(all_res,"metacycle/metacyc_res_log_tab6.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#################################################
#summarize the hits

#make a bar plot summarizing the number of cycling and non-cycling hits

sigTg<-tg_metacyc%>%filter(JTK_BH.Q<0.05) #963
sigNTg<-ntg_metacyc%>%filter(JTK_BH.Q<0.05) #652
sigTRF<-trf_metacyc%>%filter(JTK_BH.Q<0.05) #2322
write.table(sigTg,"metacycle/TgALF_metacycle/meta2d_TgALF_sighits_qless0.05_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"metacycle/NTgALF_metacycle/meta2d_NTgALF_sighits_qless0.05_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"metacycle/TgTRF_metacycle/meta2d_TgTRF_sighits_qless0.05_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

cyc_summ<-data.frame(condition=c("NTg-ALF","Tg-ALF","Tg-TRF"),
                     cycling_n=c(652,963,2322),
                     not_cycling_n=c(5079-652,5079-963,5079-2322))%>%
  mutate(cycling_perc=(cycling_n/5079)*100,
         not_cycling_perc=(not_cycling_n/5079)*100)
write.table(cyc_summ,"metacycle/SFR26_0707_cyclingsumm_log.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-ALF","Tg-TRF")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0)) +
  labs(y="cycling metabolites (% of n=5079)") +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))

ggsave("metacycle/SFR26_0707_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  #filter(condition!="NTg-ALF")%>%
  #filter(condition!="Tg-TRF")%>%
  filter(condition!="Tg-ALF")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) 
#X-squared = 830.93, df = NA, p-value = 0.0004998 tg-alf vs tf-trf
#X-squared = 71.211, df = NA, p-value = 0.0004998 ntg-alf vs tg-alf
#X-squared = 1326, df = NA, p-value = 0.0004998 ntg-alf vs tg-trf

#make a bar plot summarizing the number of cycling and non-cycling hits (removing compounds with no name)

sigTg<-tg_metacyc%>%filter(JTK_BH.Q<0.05,!is.na(Compound_Name)) #232
sigNTg<-ntg_metacyc%>%filter(JTK_BH.Q<0.05,!is.na(Compound_Name)) #169
sigTRF<-trf_metacyc%>%filter(JTK_BH.Q<0.05,!is.na(Compound_Name)) #385
write.table(sigTg,"metacycle/TgALF_metacycle/meta2d_TgALF_sighits_qless0.05_rmNA_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"metacycle/NTgALF_metacycle/meta2d_NTgALF_sighits_qless0.05_rmNA_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"metacycle/TgTRF_metacycle/meta2d_TgTRF_sighits_qless0.05_rmNA_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

noduf<-tg_metacyc%>%filter(!is.na(Compound_Name)) #667

cyc_summ<-data.frame(condition=c("NTg-ALF","Tg-ALF","Tg-TRF"),
                     cycling_n=c(169,232,385),
                     not_cycling_n=c(667-169,667-232,667-385))%>%
  mutate(cycling_perc=(cycling_n/667)*100,
         not_cycling_perc=(not_cycling_n/667)*100)
write.table(cyc_summ,"metacycle/SFR26_0707_cyclingsumm_rmNA_log.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-ALF","Tg-TRF")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0)) +
  labs(y="cycling metabolites (% of n=667)") +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))

ggsave("metacycle/SFR26_0707_cycling_summ_barplot_rmNA_log.pdf", plot=p,height=3, width=2.5)
#################################################

#not removing NA
list_venn <- list(NTgALF = sigNTg$FeatureID,
                  TgALF = sigTg$FeatureID,
                  TgTRF = sigTRF$FeatureID)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", FeatureID="X..i..")%>%
  left_join(.,annotations, by="FeatureID")
write.table(shared_cyc,"metacycle/cyclingmtb_vennoverlap_summ_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

ntgalf<-shared_cyc%>%filter(pattern=="NTgALF")
ntgalftrf<-shared_cyc%>%filter(pattern=="NTgALF:TgTRF")

sigTg_unique<-sigTg%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,shared_cyc,by="FeatureID")
write.table(sigTg_unique,"metacycle/TgALF_metacycle/meta2d_TgALF_sighits_qless0.05_woverlapinfo_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigNTg_unique<-sigNTg%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,shared_cyc,by="FeatureID")
write.table(sigNTg_unique,"metacycle/NTgALF_metacycle/meta2d_NTgALF_sighits_qless0.05_woverlapinfo_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigTRF_unique<-sigTRF%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,shared_cyc,by="FeatureID")
write.table(sigTRF_unique,"metacycle/TgTRF_metacycle/meta2d_TgTRF_sighits_qless0.05_woverlapinfo_log.txt",sep = "\t",row.names = FALSE,quote=FALSE)

draw.venn(sigNTg$FeatureID,sigTg$FeatureID,sigTRF$FeatureID,
          title="",subtitle="",
          xtitle="NTg-ALF",ytitle="Tg-ALF",ztitle="Tg-TRF",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="darkorange",y_c="purple4",z_c="forestgreen",
          output="pdf", 
          filename="metacycle/SFR25_0129_venn.pdf")

#################################################

sem <- function(x) sd(x)/sqrt(length(x))

ntgalftrf<-shared_cyc%>%filter(pattern=="NTgALF:TgTRF")%>%
  filter(!is.na(Compound_Name))

ntgalf<-shared_cyc%>%filter(pattern=="NTgALF")%>%
  filter(!is.na(Compound_Name))

mtbZT<-fread("data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_log.csv")%>%
  dplyr::rename(FeatureID=row.ID)%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  dplyr::select(FeatureID,all_of(md$sampleid))%>%
  #filter(FeatureID %in% ntgalftrf$FeatureID)%>% 
  filter(FeatureID %in% ntgalf$FeatureID)%>% 
  gather(sampleid,mtb_counts,-FeatureID) %>%
  left_join(.,annotations, by ="FeatureID") %>%
  mutate(label_name=ifelse(is.na(Compound_Name),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))%>%
  left_join(.,md,by="sampleid") %>%
  group_by(FeatureID,label_name,ATTRIBUTE_geno_n_feedingcond,ATTRIBUTE_timepoint)%>%summarise(mn_count=mean(mtb_counts),sem=sem(mtb_counts))%>%
  mutate(ATTRIBUTE_geno_n_feedingcond=factor(ATTRIBUTE_geno_n_feedingcond, levels = c("NTg-ALF","Tg-ALF","Tg-TRF")))

ZT_dist <- function(data,FeatureID) {
  
  lab<-unique((data%>%filter(FeatureID==FeatureID))$label_name)
  
  data%>%
    ggplot(aes(x=ATTRIBUTE_timepoint, y=mn_count, color=ATTRIBUTE_geno_n_feedingcond)) +
    geom_point(alpha=1.0) + geom_line() +
    theme_pubr() +
    scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
    scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
    scale_x_continuous(breaks=c(0,4,8,12,16,20))+
    geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=ATTRIBUTE_geno_n_feedingcond),alpha=0.3,colour = NA)+
    labs(color="genotype",
         fill="genotype",
         y ="log(peak abun)",
         x ="ZT time",
         title = lab)+
         #title = stringr::str_wrap(lab, width=30))+ 
    theme(plot.title = element_text(face = "bold",hjust = 0.5))
  #ggsave(paste("metacycle/ntgalftgtrfannot28_log/",FeatureID,".pdf",sep=""), height=3.5, width=4)
  ggsave(paste("metacycle/ntgalfannot7_log/",FeatureID,".pdf",sep=""), height=3.5, width=4)
}


mtbZT_nested <- mtbZT %>% 
  group_by(FeatureID) %>% 
  nest()

mtbZT_plots <- 
  mtbZT_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ ZT_dist(.x,.y)))

####################################################

#plot amine conjugated bile acids together
# sem <- function(x) sd(x)/sqrt(length(x))
# sel_BA<-c("287_i","28682","24695","24814","628_i","24286","25780","24591","900_i","537_i","574_i")
# mtbZT<-fread("data_files/FBMN_simone/gnps_quant_collapseiin_sum_corrgroup_rclr.csv")%>%
# #mtbZT<-fread("data_files/FBMN_simone/gnps_quant_collapseiin_sum_corrgroup.csv")%>%
#   dplyr::rename(FeatureID=row.ID)%>%
#   dplyr::select(FeatureID,all_of(md$sampleid))%>%
#   filter(FeatureID %in% sel_BA)%>% 
#   gather(sampleid,mtb_counts,-FeatureID) %>%
#   left_join(.,annotations, by ="FeatureID") %>%
#   mutate(Compound_Name=ifelse(FeatureID=="537_i","Ser-CA",Compound_Name))%>%
#   mutate(label_name=ifelse(is.na(Compound_Name),
#                            paste("mtb",FeatureID, "m/z", signif(mn_mz,digits=4), sep=" "),
#                            paste("mtb",FeatureID, "m/z", signif(mn_mz,digits=4), Compound_Name, sep=" ")))%>%
#   left_join(.,md,by="sampleid") %>%
#   group_by(FeatureID,label_name,ATTRIBUTE_geno_n_feedingcond,ATTRIBUTE_timepoint)%>%summarise(mn_count=mean(mtb_counts),sem=sem(mtb_counts))%>%
#   mutate(ATTRIBUTE_geno_n_feedingcond=factor(ATTRIBUTE_geno_n_feedingcond, levels = c("NTg-ALF","Tg-ALF","Tg-TRF")))
# 
# p<-ggplot(mtbZT,aes(x=ATTRIBUTE_timepoint, y=mn_count, color=label_name)) +
#   geom_point(alpha=1.0) + geom_line() +
#   theme_pubr() +
#   scale_x_continuous(breaks=c(0,4,8,12,16,20))+
#   facet_grid(~ATTRIBUTE_geno_n_feedingcond)+
#   #geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=ATTRIBUTE_geno_n_feedingcond),alpha=0.3,colour = NA)+
#   labs(color="genotype",
#        fill="genotype",
#        y ="RCLR",
#        x ="ZT time")+
#   theme(legend.position = "right")
# ggsave("metacycle/amine_cong_BA_plt_peakabun.pdf", p, height=4, width=10)
# ggsave("metacycle/amine_cong_BA_plt_log10peakabun.pdf", p, height=4, width=10)
# ggsave("metacycle/amine_cong_BA_plt_rclr.pdf", p, height=4, width=10)
# 
