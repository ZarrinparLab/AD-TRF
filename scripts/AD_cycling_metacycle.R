setwd("~/scratch/AD_EM/metagenomics/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggvenn)
library(gplots)
library(viridis)
library(RColorBrewer)
library(BioVenn) 
library("qiime2R")
library("ggdendro")
library(ggbreak)
##########################################################
#load data

df_species<-fread("data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)

md <- fread("data_files/stool.metagenomic.metadata.with_amyloid_saczt.tsv")%>%
#md <- fread("data_files/overly_specified_metadata.tsv")%>%
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"),
         sample_name =sub("^[^.]*\\.", "", sample_name),
         mouseid = sub("\\..*", "",sample_name))

##########################################################
#just tg-ALF
mdTg<-md%>%filter(genotype_diet=="Tg.ALF")%>%
  mutate(zt_time=factor(zt_time))%>%
  arrange(zt_time)

mdTg%>%group_by(zt_collection_time)%>%summarize(n=n())

Tg<-df_species%>%dplyr::select(taxonomy,mdTg$sample_name)
write.table(Tg,"data_files/159013_wolR2_per_genome.feature-table.species.clr.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just ntg-ALF
mdNTg<-md%>%filter(genotype_diet=="Ntg.ALF")%>%
  mutate(zt_time=factor(zt_time))%>%
  arrange(zt_time)

mdNTg%>%group_by(zt_collection_time)%>%summarize(n=n())

NTg<-df_species%>%dplyr::select(taxonomy,mdNTg$sample_name)
write.table(NTg,"data_files/159013_wolR2_per_genome.feature-table.species.clr.NTgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just tg-TRF
mdTRF<-md%>%filter(genotype_diet=="Tg.TRF")%>%
  mutate(zt_time=factor(zt_time))%>%
  arrange(zt_time)

mdTRF%>%group_by(zt_collection_time)%>%summarize(n=n())

TRF<-df_species%>%dplyr::select(taxonomy,mdTRF$sample_name)
write.table(TRF,"data_files/159013_wolR2_per_genome.feature-table.species.clr.TgTRF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
##########################################################
#run metacycle

#tgalf
meta2d(infile="data_files/159013_wolR2_per_genome.feature-table.species.clr.TgALF.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="cycling_metacycle/TgALF_metacycle",timepoints=rep(seq(0, 20, by=4),times=c(9,9,8,9,9,9)),minper=20,maxper=24)

#ntgalf
meta2d(infile="data_files/159013_wolR2_per_genome.feature-table.species.clr.NTgALF.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="cycling_metacycle/NTgALF_metacycle",timepoints=rep(seq(0, 20, by=4),times=c(7,7,7,7,7,7)),minper=20,maxper=24)

#tgtrf
meta2d(infile="data_files/159013_wolR2_per_genome.feature-table.species.clr.TgTRF.txt", cycMethod=c("JTK","LS"),filestyle="txt", 
       outdir="cycling_metacycle/TgTRF_metacycle",timepoints=rep(seq(0, 20, by=4),times=c(9,9,5,9,9,9)),minper=20,maxper=24)
##########################################################

#read in Tg-ALF metacycle hits
tg_metacyc<-fread("cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.clr.TgALF.txt")%>%
  dplyr::rename(taxonomy=CycID) %>%
  arrange(JTK_adjphase)

# tg_metacyc<-fread("cycling_metacycle/meta2d_tg_alf.clr.species.csv")%>%
#   dplyr::rename(taxonomy=CycID) %>%
#   arrange(JTK_adjphase)

#read in NTg-ALF metacycle hits
ntg_metacyc<-fread("cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.clr.NTgALF.txt")%>%
  dplyr::rename(taxonomy=CycID) %>%
  arrange(JTK_adjphase)

# ntg_metacyc<-fread("cycling_metacycle/meta2d_ntg_alf.clr.species.csv")%>%
#   dplyr::rename(taxonomy=CycID) %>%
#   arrange(JTK_adjphase)

#read in Tg-TRF metacycle hits
trf_metacyc<-fread("cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.clr.TgTRF.txt")%>%
  dplyr::rename(taxonomy=CycID) %>%
  arrange(JTK_adjphase)

# trf_metacyc<-fread("cycling_metacycle/meta2d_tg_trf.clr.species.csv")%>%
#   dplyr::rename(taxonomy=CycID) %>%
#   arrange(JTK_adjphase)

##########################################################

sigTg<-tg_metacyc%>%filter(JTK_BH.Q<0.05) #11
sigNTg<-ntg_metacyc%>%filter(JTK_BH.Q<0.05) #187
sigTRF<-trf_metacyc%>%filter(JTK_BH.Q<0.05) #4371
write.table(sigTg,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.05.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.05.NTgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.05.TgTRF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigTg<-tg_metacyc%>%filter(JTK_BH.Q<0.1) #3094
sigNTg<-ntg_metacyc%>%filter(JTK_BH.Q<0.1) #339
sigTRF<-trf_metacyc%>%filter(JTK_BH.Q<0.1) #4766
write.table(sigTg,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1.NTgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1.TgTRF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigTg<-tg_metacyc%>%filter(JTK_pvalue<0.05) #3125
sigNTg<-ntg_metacyc%>%filter(JTK_pvalue<0.05) #795
sigTRF<-trf_metacyc%>%filter(JTK_pvalue<0.05) #4511
write.table(sigTg,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05.NTgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05.TgTRF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

# sigTg<-tg_metacyc%>%filter(meta2d_BH.Q<0.05) #14
# sigNTg<-ntg_metacyc%>%filter(meta2d_BH.Q<0.05) #164
# sigTRF<-trf_metacyc%>%filter(meta2d_BH.Q<0.05) #3355
# sigTg<-tg_metacyc%>%filter(meta2d_BH.Q<0.1) #14
# sigNTg<-ntg_metacyc%>%filter(meta2d_BH.Q<0.1) #222
# sigTRF<-trf_metacyc%>%filter(meta2d_BH.Q<0.1) #3750
# sigTg<-tg_metacyc%>%filter(meta2d_pvalue<0.05) #244
# sigNTg<-ntg_metacyc%>%filter(meta2d_pvalue<0.05) #480
# sigTRF<-trf_metacyc%>%filter(meta2d_pvalue<0.05) #3669

cyc_summ<-data.frame(condition=c("NTg-ALF","Tg-ALF","Tg-TRF"),
                     cycling_n=c(nrow(sigNTg),nrow(sigTg),nrow(sigTRF)),
                     not_cycling_n=c(nrow(trf_metacyc)-nrow(sigNTg),nrow(trf_metacyc)-nrow(sigTg),nrow(trf_metacyc)-nrow(sigTRF)))%>%
  mutate(cycling_perc=(cycling_n/nrow(trf_metacyc))*100,
         not_cycling_perc=(not_cycling_n/nrow(trf_metacyc))*100)
write.table(cyc_summ,"cycling_metacycle/SFR25_0507_cyclingsumm_qless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(cyc_summ,"cycling_metacycle/SFR25_0507_cyclingsumm_qless0.1.txt",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(cyc_summ,"cycling_metacycle/SFR25_0507_cyclingsumm_pless0.05.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-ALF","Tg-TRF")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0)) +
  labs(y="cycling metabolites (% of n=6303)") +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))

ggsave("cycling_metacycle/SFR25_0307_cycling_summ_barplot_qless0.05.pdf", plot=p,height=3, width=2.2)
ggsave("cycling_metacycle/SFR25_0307_cycling_summ_barplot_qless0.1.pdf", plot=p,height=3, width=2.2)
ggsave("cycling_metacycle/SFR25_0307_cycling_summ_barplot_pless0.05.pdf", plot=p,height=3, width=2.2)

AD_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  #filter(condition!="NTg-ALF")%>%
  #filter(condition!="Tg-TRF")%>%
  filter(condition!="Tg-ALF")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(AD_chi,simulate.p.value=TRUE, B=2000)
#X-squared = 944.71, df = NA, p-value = 0.0004998 tg-alf v tg-trf
#X-squared = 3038.3, df = NA, p-value = 0.0004998 ntg-alf v tg-alf
#X-squared = 6451.8, df = NA, p-value = 0.0004998 ntg-alf v tg-trf
##########################################################

#not removing NA
list_venn <- list(NTgALF = sigNTg$taxonomy,
                  TgALF = sigTg$taxonomy,
                  TgTRF = sigTRF$taxonomy)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", taxonomy="X..i..")
#write.table(shared_cyc,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05woverlap.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(shared_cyc,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1woverlap.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

ntgalf<-shared_cyc%>%filter(pattern=="NTgALF")
ntgalftrf<-shared_cyc%>%filter(pattern=="NTgALF:TgTRF")

sigTg_unique<-sigTg%>%
  left_join(.,shared_cyc,by="taxonomy")
#write.table(sigTg,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05woverlap.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTg,"cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1woverlap.TgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigNTg_unique<-sigNTg%>%
  left_join(.,shared_cyc,by="taxonomy")
#write.table(sigNTg,"cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05woverlap.NTgALFwoverlap.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1woverlap.NTgALFwoverlap.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigTRF_unique<-sigTRF%>%
  left_join(.,shared_cyc,by="taxonomy")
#write.table(sigTRF,"cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.pless0.05woverlap.TgTRF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.qless0.1woverlap.TgTRF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

draw.venn(sigNTg$taxonomy,sigTg$taxonomy,sigTRF$taxonomy,
          title="",subtitle="",
          xtitle="NTg-ALF",ytitle="Tg-ALF",ztitle="Tg-TRF",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="darkorange",y_c="purple4",z_c="forestgreen",
          output="pdf", 
          filename="cycling_metacycle/SFR25_0307_vennq0.1.pdf")

##########################################################

ntgalftgtrf_heatmap<-df_species %>%
  filter(taxonomy %in% ntgalftrf$taxonomy)%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr))%>%
  group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_clr - mean(mn_clr))/sd(mn_clr))%>%
  mutate(zt_time=factor(zt_time, levels = c("0", "4","8","12","16","20")))%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))

# Run clustering
sp_df<-df_species%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  filter(taxonomy %in% ntgalftrf$taxonomy)%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  spread(sample_name,clr)

sp_matrix<-sp_df%>%column_to_rownames("taxonomy") %>%as.matrix()
sp_dendro <- as.dendrogram(hclust(d = dist(x = sp_matrix)))
dendro_plot <- ggdendrogram(data = sp_dendro, rotate = TRUE)

sp_order <- order.dendrogram(sp_dendro)
ntgalftgtrf_heatmap$taxonomy <- factor(x = ntgalftgtrf_heatmap$taxonomy,
                                       levels = sp_df$taxonomy[sp_order], 
                                       ordered = TRUE)

plt<-ggplot(ntgalftgtrf_heatmap,aes(x=zt_time, y=taxonomy)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~genotype_diet,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 
ggsave("cycling_metacycle/ntgalftgtrf_heatmap.pdf", plot=plt,height=20, width=10)

#arrange by JTK_adjphase

tgtrf<-sigTRF_unique%>%
  filter(pattern == "TgTRF" & str_detect(taxonomy, "g__Limosilactobacillus"))
  #filter(pattern == "TgTRF" & str_detect(taxonomy, "g__Prevotella"))

ntgalftgtrf_heatmap<-df_species %>%
  filter(taxonomy %in% ntgalftrf$taxonomy)%>%
  #filter(taxonomy %in% tgtrf$taxonomy)%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr))%>%
  group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_clr - mean(mn_clr))/sd(mn_clr))%>%
  mutate(zt_time=factor(zt_time, levels = c("0", "4","8","12","16","20")))%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))

# Run clustering

sigNTg_unique_htmp<-sigNTg_unique%>%
  filter(pattern=="NTgALF:TgTRF")%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  arrange(JTK_adjphase)

sigNTg_unique_htmp<-sigTRF_unique%>%
  #filter(pattern == "TgTRF" & str_detect(taxonomy, "g__Prevotella"))%>%
  filter(pattern == "TgTRF" & str_detect(taxonomy, "g__Limosilactobacillus"))%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  arrange(JTK_adjphase)

ntgalftgtrf_heatmap$taxonomy<-factor(ntgalftgtrf_heatmap$taxonomy,levels = sigNTg_unique_htmp$taxonomy)

plt<-ggplot(ntgalftgtrf_heatmap,aes(x=zt_time, y=taxonomy)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~genotype_diet,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 
#ggsave("cycling_metacycle/ntgalftgtrf_heatmap_byphase.pdf", plot=plt,height=20, width=10)
#ggsave("cycling_metacycle/ntgalftgtrf_heatmap_clnbyphase.pdf", plot=plt,height=8, width=9)
ggsave("cycling_metacycle/tgtrflimo_heatmap_clnbyphase.pdf", plot=plt,height=2.5, width=9)
ggsave("cycling_metacycle/tgtrfprev_heatmap_clnbyphase.pdf", plot=plt,height=4, width=9)
#plot the 22 species with >9 phase ntg-alf

Ntg22sharedcyc<-sigNTg_unique_htmp%>%filter(JTK_adjphase>9)%>%
  dplyr::select(taxonomy,JTK_adjphase)%>%
  mutate(condition="NTg-ALF")

trf22sharedcyc<-sigTRF_unique%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  filter(taxonomy%in% Ntg22sharedcyc$taxonomy)%>%
  dplyr::select(taxonomy,JTK_adjphase)%>%
  mutate(condition="Tg-TRF")

alf22sharedcyc<-tg_metacyc%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  filter(taxonomy%in% Ntg22sharedcyc$taxonomy)%>%
  dplyr::select(taxonomy,JTK_adjphase)%>%
  mutate(condition="Tg-ALF")

# comb_22sharedcyc<-rbind(Ntg22sharedcyc,trf22sharedcyc,alf22sharedcyc)%>%
#   mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-ALF","Tg-TRF")))

comb_22sharedcyc<-rbind(Ntg22sharedcyc,trf22sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-TRF")))

#write.table(comb_22sharedcyc,"cycling_metacycle/ntgalftgtrf_densityphase_clust1_all.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(comb_22sharedcyc,"cycling_metacycle/ntgalftgtrf_densityphase_clust1.txt",sep = "\t",row.names = FALSE, quote=FALSE)

comb_22sharedcyc<-fread("cycling_metacycle/ntgalftgtrf_densityphase_clust1.txt")

ggplot(comb_22sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0), limits = c(0,0.5)) +
  scale_x_continuous(expand=c(0,0),limits = c(0,24)) +
  # scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  # scale_fill_manual(values=c("darkorange","purple4","forestgreen"))
  scale_color_manual(values=c("darkorange","forestgreen"))+
  scale_fill_manual(values=c("darkorange","forestgreen"))

#ggsave("cycling_metacycle/ntgalftgtrf_densityphase_clust1_all.pdf", width = 4.5, height = 4)
#ggsave("cycling_metacycle/ntgalftgtrf_densityphase_clust1.pdf", width = 4.5, height = 4)
ggsave("cycling_metacycle/ntgalftgtrf_densityphase_clust1_extscale.pdf", width = 4.5, height = 4)

#plot the other species with phase ntg-alf

Ntg115sharedcyc<-sigNTg_unique_htmp%>%filter(JTK_adjphase<10)%>%
  dplyr::select(taxonomy,JTK_adjphase)%>%
  mutate(condition="NTg-ALF")

trf115sharedcyc<-sigTRF_unique%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  filter(taxonomy%in% Ntg115sharedcyc$taxonomy)%>%
  dplyr::select(taxonomy,JTK_adjphase)%>%
  mutate(condition="Tg-TRF")

alf115sharedcyc<-tg_metacyc%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  filter(taxonomy%in% Ntg115sharedcyc$taxonomy)%>%
  dplyr::select(taxonomy,JTK_adjphase)%>%
  mutate(condition="Tg-ALF")

comb_115sharedcyc<-rbind(Ntg115sharedcyc,trf115sharedcyc,alf115sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-ALF","Tg-TRF")))

comb_115sharedcyc<-rbind(Ntg115sharedcyc,trf115sharedcyc)%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-TRF")))

write.table(comb_115sharedcyc,"cycling_metacycle/ntgalftgtrf_densityphase_clust2_all.txt",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(comb_115sharedcyc,"cycling_metacycle/ntgalftgtrf_densityphase_clust2.txt",sep = "\t",row.names = FALSE, quote=FALSE)


ggplot(comb_115sharedcyc, aes(x=JTK_adjphase, color=condition,fill=condition)) +
  geom_density(alpha=0.3)+ theme_pubr()+
  scale_y_continuous(expand=c(0,0), limits = c(0,2)) +
  scale_x_continuous(expand=c(0,0),limits = c(0,24)) +
  # scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  # scale_fill_manual(values=c("darkorange","purple4","forestgreen")) 
  scale_color_manual(values=c("darkorange","forestgreen"))+
  scale_fill_manual(values=c("darkorange","forestgreen")) 

#ggsave("cycling_metacycle/ntgalftgtrf_densityphase_clust2_all.pdf", width = 4.5, height = 4)
ggsave("cycling_metacycle/ntgalftgtrf_densityphase_clust2.pdf", width = 4.5, height = 4)

#summarize the top microbes shared

summ_shared<-sigNTg_unique_htmp<-sigNTg_unique%>%
  filter(pattern=="NTgALF:TgTRF", JTK_adjphase<10)%>%
  #filter(pattern=="NTgALF:TgTRF", JTK_adjphase>9)%>%
  dplyr::select(taxonomy)%>%
  mutate(taxonomy = sub(".*(f__.*)", "\\1",taxonomy))%>%
  tidyr::separate(taxonomy, into = c("family","genus","species"), sep = ";")%>%
  group_by(family)%>%
  summarise(n=n())%>%
  #top_n(5)%>%
  filter(n>1)%>%
  arrange(n)

summ_shared$family<-factor(summ_shared$family,levels = summ_shared$family[order(summ_shared$n)])
  

p<-ggplot(data=summ_shared, aes(x=family, y=n)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0)) +
  labs(y="num. of cycling microbes") + coord_flip()
#ggsave("cycling_metacycle/ntgalftgtrf_clust2_top5fam.pdf", width = 3.5, height = 3)
ggsave("cycling_metacycle/ntgalftgtrf_clust1_top5fam.pdf", width = 3.5, height = 3)

#################################################

#summary of unique hits per condition
sigNTg_unique_sum<-sigNTg_unique%>%
  filter(pattern=="NTgALF")%>%
  dplyr::select(taxonomy)%>%
  mutate(taxonomy = sub(".*(f__.*)", "\\1",taxonomy))%>%
  tidyr::separate(taxonomy, into = c("family","genus","species"), sep = ";")%>%
  group_by(family)%>%
  summarise(n=n())%>%
  filter(n>3)%>%
  arrange(n)

sigNTg_unique_sum$family<-factor(sigNTg_unique_sum$family,levels = sigNTg_unique_sum$family[order(sigNTg_unique_sum$n)])

p<-ggplot(data=sigNTg_unique_sum, aes(x=family, y=n)) +
  geom_bar(stat="identity", fill="darkorange") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0),limits=c(0,65),breaks=seq(0, 60, by=10)) +
  labs(y="num. of cycling species") + coord_flip()
ggsave("cycling_metacycle/uniquentgalf_summplot.pdf", width = 4, height = 3)

sigTg_unique_sum<-sigTg_unique%>%
  filter(pattern=="TgALF")%>%
  dplyr::select(taxonomy)%>%
  mutate(taxonomy = sub(".*(f__.*)", "\\1",taxonomy))%>%
  tidyr::separate(taxonomy, into = c("family","genus","species"), sep = ";")%>%
  group_by(family)%>%
  summarise(n=n())%>%
  filter(n>3)%>%
  arrange(n)

sigTg_unique_sum$family<-factor(sigTg_unique_sum$family,levels = sigTg_unique_sum$family[order(sigTg_unique_sum$n)])

p<-ggplot(data=sigTg_unique_sum, aes(x=family, y=n)) +
  geom_bar(stat="identity", fill="purple4") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0),limits=c(0,20),breaks=seq(0, 20, by=5)) +
  labs(y="num. of cycling species") + coord_flip()
ggsave("cycling_metacycle/uniquetgalf_summplot_expand.pdf", width = 4, height = 3)

sigTRF_unique_sum<-sigTRF_unique%>%
  filter(pattern=="TgTRF")%>%
  dplyr::select(taxonomy)%>%
  mutate(taxonomy = sub(".*(f__.*)", "\\1",taxonomy))%>%
  tidyr::separate(taxonomy, into = c("family","genus","species"), sep = ";")%>%
  group_by(family)%>%
  summarise(n=n())%>%
  filter(n>29)%>%
  arrange(n)

sigTRF_unique_sum$family<-factor(sigTRF_unique_sum$family,levels = sigTRF_unique_sum$family[order(sigTRF_unique_sum$n)])

p<-ggplot(data=sigTRF_unique_sum, aes(x=family, y=n)) +
  geom_bar(stat="identity", fill="forestgreen") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0),limits=c(0,170),breaks=seq(0, 170, by=10)) +
  scale_y_break(c(60,160),scales=0.2, expand=FALSE) +
  labs(y="num. of cycling species") + coord_flip()
ggsave("cycling_metacycle/uniquetgtrf_summplot.pdf", width = 4, height = 3.5)
  
#################################################

sem <- function(x) sd(x)/sqrt(length(x))

ntgalftgtrf<-sigNTg_unique%>%
  filter(pattern=="NTgALF:TgTRF")

tgalf<-sigTg_unique%>%
  filter(pattern=="TgALF")

ntgalf<-sigNTg_unique%>%
  filter(pattern == "NTgALF" & str_detect(taxonomy, "f__Lachnospiraceae"))

tgtrf<-sigTRF_unique%>%
  #filter(pattern == "TgTRF" & str_detect(taxonomy, "f__Lactobacillaceae"))
  filter(pattern == "TgTRF" & str_detect(taxonomy, "f__Bacteroidaceae"))

mtbZT<-fread("data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  #filter(taxonomy %in% ntgalftgtrf$taxonomy)%>% 
  #filter(taxonomy %in% ntgalf$taxonomy)%>% 
  #filter(taxonomy %in% tgtrf$taxonomy)%>% 
  filter(taxonomy %in% tgalf$taxonomy)%>% 
  left_join(.,md,by="sample_name") %>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%summarise(mn_count=mean(clr),sem=sem(clr))%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF","Tg.ALF","Tg.TRF")))

ZT_dist <- function(data,taxonomy) {
  
  #lab<-unique(data%>%filter(taxonomy==taxonomy))
  
  data%>%
    ggplot(aes(x=zt_time, y=mn_count, color=genotype_diet)) +
    geom_point(alpha=1.0) + geom_line() +
    theme_pubr() +
    scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
    scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
    scale_x_continuous(breaks=c(0,4,8,12,16,20))+
    geom_ribbon(aes(ymin = mn_count-sem, ymax = mn_count+sem, fill=genotype_diet),alpha=0.3,colour = NA)+
    labs(color="genotype",
         fill="genotype",
         y ="CLR",
         x ="ZT time",
         title = taxonomy)+
    #title = stringr::str_wrap(lab, width=30))+ 
    theme(plot.title = element_text(face = "bold",hjust = 0.5))
  #ggsave(paste("cycling_metacycle/ntgalftgtrf137/",taxonomy,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("cycling_metacycle/ntgalf_lachno/",taxonomy,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("cycling_metacycle/tgtrf_lacto/",taxonomy,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("cycling_metacycle/tgtrf_prev/",taxonomy,".pdf",sep=""), height=4, width=4)
  ggsave(paste("cycling_metacycle/tgalf169/",taxonomy,".pdf",sep=""), height=4, width=4)
}


mtbZT_nested <- mtbZT %>% 
  group_by(taxonomy) %>% 
  nest()

mtbZT_plots <- 
  mtbZT_nested %>% 
  mutate(plot = map2(data, taxonomy,  ~ ZT_dist(.x,.y)))

##########################################################
#create summary of cycling results

ntgalf<-fread("cycling_metacycle/NTgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.clr.NTgALF.txt")%>%
  mutate(conditions="NTg-ALF")

tgalf<-fread("cycling_metacycle/TgALF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.clr.TgALF.txt")%>%
  mutate(conditions="Tg-ALF")

tgtrf<-fread("cycling_metacycle/TgTRF_metacycle/meta2d_159013_wolR2_per_genome.feature-table.species.clr.TgTRF.txt")%>%
  mutate(conditions="Tg-TRF")

allcycres<-rbind(ntgalf,tgalf,tgtrf)%>%
  arrange(JTK_BH.Q)

write.table(allcycres,"cycling_metacycle/metag_micro_metacycle_allres.txt",sep = "\t",row.names = FALSE,quote=FALSE) #Table S2

