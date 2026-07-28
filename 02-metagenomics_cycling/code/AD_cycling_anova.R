setwd("~/scratch/AD_EM/metagenomics/")

library(MetaCycle)
library(tidyverse)
library(data.table)
library(ggvenn)
library(VennDiagram) 
library(grDevices)
library(gplots)
library(ggpubr)
library(BioVenn) 
library("ggdendro")
#################################################
#load data

df_species<-fread("data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  #mutate(taxonomy = sub(".*(s__.*)", "\\1",taxonomy))%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))

md <- fread("data_files/overly_specified_metadata.tsv")%>%
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"),
         mouseid = sub("\\..*", "",sample_name))

anova_allmicrobes<-function(df){
  mod<-aov(clr~ zt_time, data = df)
  res<-data.frame(zt_time=summary(mod)[[1]][[5]][[1]])
  return(res)
}

#################################################
# run ANOVA

#just tg-ALF
mdTg<-md%>%filter(genotype_diet=="Tg.ALF")%>%
  mutate(zt_time=factor(zt_time))%>%
  arrange(zt_time)

Tg<-df_species%>%filter(sample_name %in% mdTg$sample_name)%>%
  left_join(.,mdTg,by="sample_name")

all_Tg<-Tg %>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, anova_allmicrobes))

Tg_anova<-all_Tg%>%dplyr::select(-data)%>%unnest()
write.table(Tg_anova,"cycling/SFR25_0405_TgALFzt_anova.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just ntg-ALF
mdNTg<-md%>%filter(genotype_diet=="Ntg.ALF")%>%
  mutate(zt_time=factor(zt_time))%>%
  arrange(zt_time)

NTg<-df_species%>%filter(sample_name %in% mdNTg$sample_name)%>%
  left_join(.,mdNTg,by="sample_name")

all_NTg<-NTg %>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, anova_allmicrobes))

NTg_anova<-all_NTg%>%dplyr::select(-data)%>%unnest()
write.table(NTg_anova,"cycling/SFR25_0405_NTgALFzt_anova.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just tg-TRF
mdTRF<-md%>%filter(genotype_diet=="Tg.TRF")%>%
  mutate(zt_time=factor(zt_time))%>%
  arrange(zt_time)

TRF<-df_species%>%filter(sample_name %in% mdTRF$sample_name)%>%
  left_join(.,mdTRF,by="sample_name")

all_TRF<-TRF %>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, anova_allmicrobes))

TRF_anova<-all_TRF%>%dplyr::select(-data)%>%unnest()
write.table(TRF_anova,"cycling/SFR25_0305_TgTRFzt_anova.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#################################################
#summarize the hits

#make a bar plot summarizing the number of cycling and non-cycling hits

sigTg<-Tg_anova%>%filter(zt_time<0.05) #103
sigNTg<-NTg_anova%>%filter(zt_time<0.05) #214
sigTRF<-TRF_anova%>%filter(zt_time<0.05) #3500
write.table(sigTg,"cycling/SFR25_0305_TgALFzt_anova_sig0.05.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigNTg,"cycling/SFR25_0305_NTgALFzt_anova_sig0.05.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(sigTRF,"cycling/SFR25_0305_TgTRFzt_anova_sig0.05.txt",sep = "\t",row.names = FALSE,quote=FALSE)

cyc_summ<-data.frame(condition=c("NTg-ALF","Tg-ALF","Tg-TRF"),
                     cycling_n=c(214,103,3500),
                     not_cycling_n=c(6303-214,6303-103,6303-3500))%>%
  mutate(cycling_perc=(cycling_n/6303)*100,
         not_cycling_perc=(not_cycling_n/6303)*100)
write.table(cyc_summ,"cycling/SFR25_0305_cyclingsumm.txt",sep = "\t",row.names = FALSE, quote=FALSE)

cyc_summ_plot<-cyc_summ%>%gather(cycling_grp, value, -condition)%>%
  filter(cycling_grp=="cycling_perc"|cycling_grp=="not_cycling_perc")%>%
  mutate(condition=factor(condition,levels=c("NTg-ALF","Tg-ALF","Tg-TRF")))%>%
  filter(cycling_grp=="cycling_perc")

p<-ggplot(data=cyc_summ_plot, aes(x=condition, y=value, fill=condition)) +
  geom_bar(stat="identity") + theme_classic() +
  theme(legend.position = "top")+
  scale_y_continuous(expand=c(0,0)) +
  labs(y="species with temporal changes (% of n=6303)") +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))

ggsave("cycling/SFR25_0305_cycling_summ_barplot.pdf", plot=p,height=3, width=2.2)

mtx_chi<-cyc_summ%>%dplyr::select(1:3)%>%
  filter(condition!="NTg-ALF")%>%
  #filter(condition!="Tg-TRF")%>%
  #filter(condition!="Tg-ALF")%>%
  column_to_rownames("condition")%>%as.matrix()

chisq.test(mtx_chi,simulate.p.value=TRUE, B=2000) 
#X-squared = 4484.5, df = NA, p-value = 0.0004998 tg-alf vs tf-trf
#X-squared = 39.87, df = NA, p-value = 0.0004998 ntg-alf vs tg-alf
#X-squared = 4121.6, df = NA, p-value = 0.0004998 ntg-alf vs tg-trf
#################################################

list_venn <- list(NTgALF = sigNTg$taxonomy,
                  TgALF = sigTg$taxonomy,
                  TgTRF = sigTRF$taxonomy)

ItemsList <- venn(list_venn, show.plot = FALSE)
all<-attributes(ItemsList)$intersections

shared_cyc<-plyr::ldply (all, data.frame)%>%
  dplyr::rename(pattern=".id", taxonomy="X..i..")
write.table(shared_cyc,"cycling/cyclingmtb_vennoverlap_summ.txt",sep = "\t",row.names = FALSE,quote=FALSE)

ntgalf<-shared_cyc%>%filter(pattern=="NTgALF")
ntgalftrf<-shared_cyc%>%filter(pattern=="NTgALF:TgTRF")

sigTg_unique<-sigTg%>%
  left_join(.,shared_cyc,by="taxonomy")
write.table(sigTg_unique,"cycling/TgALF_sighits_qless0.05_woverlapinfo.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigNTg_unique<-sigNTg%>%
  left_join(.,shared_cyc,by="taxonomy")
write.table(sigNTg_unique,"cycling/NTgALFsighits_qless0.05_woverlapinfo.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sigTRF_unique<-sigTRF%>%
  left_join(.,shared_cyc,by="taxonomy")
write.table(sigTRF_unique,"cycling/TgTRF_sighits_qless0.05_woverlapinfo.txt",sep = "\t",row.names = FALSE,quote=FALSE)

draw.venn(sigNTg$taxonomy,sigTg$taxonomy,sigTRF$taxonomy,
          title="",subtitle="",
          xtitle="NTg-ALF",ytitle="Tg-ALF",ztitle="Tg-TRF",
          xt_s=3,yt_s=3,zt_s=3, nr_s=3,
          x_c="darkorange",y_c="purple4",z_c="forestgreen",
          output="pdf", 
          filename="cycling/SFR25_0305_venn.pdf")

#################################################

sem <- function(x) sd(x)/sqrt(length(x))

ntgalftrf<-shared_cyc%>%filter(pattern=="NTgALF:TgTRF")
tgtrf<-fread("cycling/TgTRF_sighits_qless0.05_woverlapinfo.txt")%>%
  filter(pattern=="TgTRF")%>%
  #filter(str_detect(taxonomy,"Limosilactobacillus"))
  #filter(str_detect(taxonomy,"Ligilactobacillus"))
  filter(str_detect(taxonomy,"Prevotella"))
mtbZT<-fread("data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  #filter(taxonomy %in% ntgalftrf$taxonomy)%>% 
  #filter(taxonomy %in% tgtrf$taxonomy)%>% 
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
  #ggsave(paste("cycling/ntgalftgtrf82/",taxonomy,".pdf",sep=""), height=4, width=4)
  ggsave(paste("cycling/tgtrf/",taxonomy,".pdf",sep=""), height=4, width=4)
}


mtbZT_nested <- mtbZT %>% 
  group_by(taxonomy) %>% 
  nest()

mtbZT_plots <- 
  mtbZT_nested %>% 
  mutate(plot = map2(data, taxonomy,  ~ ZT_dist(.x,.y)))

####################################################
ntgalftgtrf_heatmap<-df_species %>%
  filter(taxonomy %in% ntgalftrf$taxonomy)%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr))%>%
  group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_clr - mean(mn_clr))/sd(mn_clr))%>%
  mutate(zt_time=factor(zt_time, levels = c("0", "4","8","12","16","20")))%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))

# Run clustering
sp_df<-df_species%>%
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
ggsave("cycling/ntgalftgtrf_heatmap.pdf", plot=plt,height=10, width=10)

