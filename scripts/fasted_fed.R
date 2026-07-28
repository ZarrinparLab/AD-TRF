setwd("~/scratch/AD_EM")

library(tidyverse)
library(data.table)
library(qiime2R)
library(ggpubr)
library(cowplot)
library(agricolae)
library(ggridges)
#####################################################

subset_dat<-function(df,microbe){
  dat<-df%>%
    filter(taxonomy==microbe)%>%
    left_join(.,md,by="sample_name")
  return(dat)
}

fasted_fed_plt<-function(df,microbe){
  
  p<-ggplot(df, aes(x=genotype_diet, y=clr, fill=genotype_diet)) +
    #geom_violin(alpha=0.3,trim = FALSE) + 
    geom_boxplot(alpha=0.3) +
    geom_dotplot(binaxis='y', stackdir='center',
                 position=position_dodge(1)) +
    scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
    facet_grid(~fasted)+
    theme_pubr()+
    labs(x="condition",y="CLR abundance", title=microbe)+
    theme(legend.position = "none")
  
  return(p)
}


fasted_fed_pval<-function(df){
  df_fasted<-df%>%filter(fasted=="fasted")
  x<-pairwise.wilcox.test(df_fasted$clr, df_fasted$genotype_diet,p.adjust.method="fdr")
  df_fasted_pvals<-data.frame(NtgALFvTgALF=x[["p.value"]][1],
                              NtgALFvTgTRF=x[["p.value"]][2],
                              TgALFvTgTRF=x[["p.value"]][4],
                              fasted="fasted")
  df_fed<-df%>%filter(fasted=="fed")
  y<-pairwise.wilcox.test(df_fed$clr, df_fed$genotype_diet,p.adjust.method="fdr")
  df_fed_pvals<-data.frame(NtgALFvTgALF=y[["p.value"]][1],
                              NtgALFvTgTRF=y[["p.value"]][2],
                              TgALFvTgTRF=y[["p.value"]][4],
                              fasted="fed")
  pvals<-rbind(df_fasted_pvals,df_fed_pvals)
  return(pvals)
}

fasted_fed_pval2<-function(df){
  df_ntgalf<-df%>%filter(genotype_diet=="Ntg.ALF")
  x<-pairwise.wilcox.test(df_ntgalf$clr, df_ntgalf$fasted,p.adjust.method="fdr")
  df_ntgalf_pvals<-data.frame(fastedvfed=x[["p.value"]][1],
                              genotype_diet="Ntg.ALF")
  df_tgalf<-df%>%filter(genotype_diet=="Tg.ALF")
  y<-pairwise.wilcox.test(df_tgalf$clr, df_tgalf$fasted,p.adjust.method="fdr")
  df_tgalf_pvals<-data.frame(fastedvfed=y[["p.value"]][1],
                              genotype_diet="Tg.ALF")
  df_tgtrf<-df%>%filter(genotype_diet=="Tg.TRF")
  y<-pairwise.wilcox.test(df_tgtrf$clr, df_tgtrf$fasted,p.adjust.method="fdr")
  df_tgtrf_pvals<-data.frame(fastedvfed=y[["p.value"]][1],
                              genotype_diet="Tg.TRF")
  pvals<-rbind(df_ntgalf_pvals,df_tgalf_pvals,df_tgtrf_pvals)
  return(pvals)
}

get_ratio_dat<-function(df){
  dat<-df%>% 
    group_by(mouseid,fasted)%>%
    reframe(mn_counts = mean(counts))%>%
    spread(fasted,mn_counts)%>%
    mutate(ratioff=fasted/fed,
           log_ratioff=log2((fasted/fed)+0.01))%>%
    left_join(.,md_subject,by="mouseid")
  return(dat)
}

fastedfed_ratio_plt<-function(df, microbe){
  p<-ggplot(df, aes(x=genotype_diet, y=log_ratioff, fill=genotype_diet)) +
    geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                           position=position_dodge(1)) +
    scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
    theme_pubr()+
    labs(x="condition",y="log2(fasted/fed)", title=microbe)+
    theme(legend.position = "none")
  return(p)
}

anova_allmicrobes<-function(df){
  mod<-aov(clr~ genotype_diet*zt_time, data = df)
  res<-data.frame(genotype_diet=summary(mod)[[1]][[5]][1],
                  zt_time=summary(mod)[[1]][[5]][2],
                  genotype_diet_zt_time=summary(mod)[[1]][[5]][3])
  return(res)
}
#####################################################
md <- fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"),
         mouseid = sub("\\..*", "",sample_name))

md_subject<-md%>%filter(zt_time==0)

#family

df_family<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.family.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  mutate(taxonomy = sub(".*(f__.*)", "\\1",taxonomy))%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))

df_lacto<-subset_dat(df_family, "f__Lactobacillaceae")

df_lacto%>%group_by(genotype_diet,zt_collection_time)%>%
  summarise(n=n())

p<-fasted_fed_plt(df_lacto,"Lactobacillaceae")
pvals<-fasted_fed_pval(df_lacto)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.6371151 0.0004156008 4.405410e-05 fasted
# 2    0.1112974 0.0020649600 9.894045e-07    fed

lacto_mod <- aov(clr~ genotype_diet*zt_time, data = df_lacto)
summary(lacto_mod)

LSD_Test_lacto<- (LSD.test(lacto_mod, c("genotype_diet", "zt_time")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("genotype_diet", "zt_time"),sep=":")
write.table(LSD_Test_lacto,"metagenomics/SFR24_1028_ADlacto_anovaLSD_tmpt.txt",sep = "\t",row.names = FALSE,quote=FALSE)









df_bacterio<-subset_dat(df_family, "f__Bacteroidaceae")

df_bacterio%>%group_by(genotype_diet,zt_collection_time)%>%
  summarise(n=n())

p<-fasted_fed_plt(df_bacterio,"Bacteroidaceae")
pvals<-fasted_fed_pval(df_bacterio)
# NtgALFvTgALF NtgALFvTgTRF TgALFvTgTRF fasted
# 1    0.6789411    0.2495927   0.3444813 fasted
# 2    0.2721927    0.2031794   0.2721927    fed

bacterio_mod <- aov(clr~ genotype_diet*zt_time, data = df_bacterio)
summary(bacterio_mod)

LSD_Test_bacterio<- (LSD.test(bacterio_mod, c("genotype_diet", "zt_time")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("genotype_diet", "zt_time"),sep=":")
write.table(LSD_Test_bacterio,"metagenomics/SFR24_1028_ADbacterio_anovaLSD_tmpt.txt",sep = "\t",row.names = FALSE,quote=FALSE)


df_gastra<-subset_dat(df_family, "f__Gastranaerophilaceae")
p<-fasted_fed_plt(df_gastra,"Gastranaerophilaceae")
pvals<-fasted_fed_pval(df_gastra)
# NtgALFvTgALF NtgALFvTgTRF TgALFvTgTRF fasted
# 1    0.4349660   0.03553745  0.08204492 fasted
# 2    0.6067473   0.35725815  0.13750453    fed

gastra_mod <- aov(clr~ genotype_diet*zt_time, data = df_gastra)
summary(gastra_mod)

LSD_Test_gastra<- (LSD.test(gastra_mod, c("genotype_diet", "zt_time")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("genotype_diet", "zt_time"),sep=":")
write.table(LSD_Test_bacterio,"metagenomics/SFR24_1028_ADbacterio_anovaLSD_tmpt.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#####################################################

#genus

df_genus<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.genus.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))

df_genus_notnorm<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.genus.tsv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  mutate(taxonomy = sub(".*(g__.*)", "\\1",taxonomy))%>%
  gather(sample_name,counts,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  mutate(sample_name = sub("^[^.]*\\.", "", sample_name))
  

df_prev<-subset_dat(df_genus, "g__Prevotella")

df_prev%>%group_by(genotype_diet,fasted)%>%
  summarise(n=n())

p<-fasted_fed_plt(df_prev,"Prevotella")+
  scale_y_continuous(limits = c(8.5, 12), breaks = seq(8.5, 12, by = 0.5))
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotella.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prev)

# NtgALFvTgALF NtgALFvTgTRF TgALFvTgTRF fasted
# 1   0.93099804    0.2773095 0.277309489 fasted
# 2   0.02680473    0.3007068 0.003388289    fed

prev_mod <- aov(clr~ genotype_diet*zt_time, data = df_prev)
summary(prev_mod)

LSD_Test_prev<- (LSD.test(prev_mod, c("genotype_diet", "zt_time")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("genotype_diet", "zt_time"),sep=":")
write.table(LSD_Test_prev,"metagenomics/SFR24_1028_ADprev_anovaLSD_tmpt.txt",sep = "\t",row.names = FALSE,quote=FALSE)

df_prev_notnorm<-subset_dat(df_genus_notnorm, "g__Prevotella")
df_prev_ratio<-get_ratio_dat(df_prev_notnorm)
p<-fastedfed_ratio_plt(df_prev_ratio,"Prevotella")
pairwise.t.test(df_prev_ratio$log_ratioff, df_prev_ratio$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.13    -     
#   Tg.TRF 0.68    0.1
ggsave("metagenomics/SFR24_1025_fastedfedratio_Prevotella.pdf", plot=p,height=3.5, width=3)

df_limosilac<-subset_dat(df_genus, "g__Limosilactobacillus")

df_limosilac%>%group_by(genotype_diet,fasted)%>%
  summarise(n=n())

p<-fasted_fed_plt(df_limosilac,"Limosilactobacillus")+
  scale_y_continuous(limits = c(4, 12.5), breaks = seq(4, 12.5, by = 1))
ggsave("metagenomics/SFR24_1023_fastedfed_Limosilactobacillus.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_limosilac)

# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1   0.45806583 0.0002416134 6.701025e-06 fasted
# 2   0.03663654 0.2818368405 1.456268e-05    fed

limo_mod <- aov(clr~ genotype_diet*zt_time, data = df_limosilac)
summary(limo_mod)

LSD_Test_limosilac<- (LSD.test(limo_mod, c("genotype_diet", "zt_time")))$groups%>%as.data.frame()%>%rownames_to_column("group")%>%
  separate(group,c("genotype_diet", "zt_time"),sep=":")
write.table(LSD_Test_limosilac,"metagenomics/SFR24_1028_ADlimosilac_anovaLSD_tmpt.txt",sep = "\t",row.names = FALSE,quote=FALSE)

df_limosilac_notnorm<-subset_dat(df_genus_notnorm, "g__Limosilactobacillus")
df_limosilac_ratio<-get_ratio_dat(df_limosilac_notnorm)
p<-fastedfed_ratio_plt(df_limosilac_ratio,"Limosilactobacillus")
pairwise.t.test(df_limosilac_ratio$log_ratioff, df_limosilac_ratio$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.85    -     
#   Tg.TRF 0.74    0.74 

#removing outlier
# Ntg.ALF Tg.ALF
# Tg.ALF 0.040   -     
#   Tg.TRF 0.081   0.751 
ggsave("metagenomics/SFR24_1025_fastedfedratio_Limosilactobacillus.pdf", plot=p,height=3.5, width=3)



df_ligila<-subset_dat(df_genus, "g__Ligilactobacillus")

p<-fasted_fed_plt(df_ligila,"Ligilactobacillus")+
  scale_y_continuous(limits = c(7, 16), breaks = seq(8, 16, by = 2))
ggsave("metagenomics/SFR24_1023_fastedfed_Ligilactobacillus.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_ligila)

# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.5312539 0.0005126662 8.133223e-05 fasted
# 2    0.1788012 0.0026081665 1.130189e-06    fed


df_ligila_notnorm<-subset_dat(df_genus_notnorm, "g__Ligilactobacillus")
df_ligila_ratio<-get_ratio_dat(df_ligila_notnorm)
p<-fastedfed_ratio_plt(df_ligila_ratio,"Ligilactobacillus")
pairwise.t.test(df_ligila_ratio$ratioff, df_ligila_ratio$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.74    -     
#   Tg.TRF 0.74    0.74  

ggsave("metagenomics/SFR24_1025_fastedfedratio_Ligilactobacillus.pdf", plot=p,height=3.5, width=3)
#####################################################

#species

df_species<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  mutate(taxonomy = sub(".*(s__.*)", "\\1",taxonomy))%>%
  gather(sample_name,clr,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))

df_species_notnorm<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv")%>%
  dplyr::rename(taxonomy=`OTU ID`)%>%
  mutate(taxonomy = sub(".*(s__.*)", "\\1",taxonomy))%>%
  gather(sample_name,counts,-taxonomy)%>%
  filter(!str_detect(sample_name, "BLANK"))%>%
  mutate(sample_name = sub("^[^.]*\\.", "", sample_name))

all_limosilac<-df_species %>%
  filter(str_detect(taxonomy, "Limosilactobacillus"))%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, anova_allmicrobes))

limosilac_anova<-all_limosilac%>%dplyr::select(-data)%>%unnest()
write.table(limosilac_anova,"metagenomics/SFR24_1028_ADlimosilac_sp_anova_condtime_pvals.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sig_limo<-limosilac_anova%>%filter(zt_time<0.1)

limo_heatmap<-df_species %>%
#limo_heatmap<-df_species_notnorm %>%
  filter(taxonomy %in% sig_limo$taxonomy)%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr))%>%
  group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_clr - mean(mn_clr))/sd(mn_clr))%>%
  #dplyr::summarise(mn_counts=mean(counts))%>%
  #group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_counts - mean(mn_counts))/sd(mn_counts))%>%
  mutate(zt_time=factor(zt_time, levels = c("0", "4","8","12","16","20")))
  
plt<-ggplot(limo_heatmap,aes(x=zt_time, y=taxonomy)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~genotype_diet,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 
#ggsave("metagenomics/SFR24_1028_5limosig_heatmap.pdf", plot=plt,height=2, width=7)
ggsave("metagenomics/SFR24_1028_5limosig_heatmap_usingnotnormcounts.pdf", plot=plt,height=1.7, width=7)

sem <- function(x) sd(x)/sqrt(length(x))

limo_lnplt<-df_species %>%
  filter(taxonomy %in% sig_limo$taxonomy)%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr),
                   sem=sem(clr))

p <- ggplot(limo_lnplt, aes(x = zt_time, y = mn_clr, color = genotype_diet)) +
  geom_point(alpha = 1.0) +
  geom_line() +
  geom_ribbon(aes(ymin = mn_clr - sem, ymax = mn_clr + sem, fill = genotype_diet), alpha = 0.3, color = NA) +
  facet_grid(~taxonomy) +
  theme_pubr() +
  scale_color_manual(values = c("darkorange", "purple4", "forestgreen")) +
  scale_fill_manual(values = c("darkorange", "purple4", "forestgreen")) +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20))
ggsave("metagenomics/SFR24_1028_5limosig_facelinplt_condi.pdf", plot=p,height=3.5, width=10)

p <- ggplot(limo_lnplt, aes(x = zt_time, y = mn_clr, color = taxonomy)) +
  geom_point(alpha = 1.0) +
  geom_line() +
  geom_ribbon(aes(ymin = mn_clr - sem, ymax = mn_clr + sem, fill = taxonomy), alpha = 0.3, color = NA) +
  facet_grid(~genotype_diet) +
  theme_pubr() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20))
ggsave("metagenomics/SFR24_1028_5limosig_facelinplt_taxa.pdf", plot=p,height=3.5, width=10)

limo_LDplt<-df_species %>%
  filter(taxonomy %in% sig_limo$taxonomy)%>%
  left_join(.,md,by="sample_name")

p<-ggplot(limo_LDplt, aes(x=fasted, y=clr, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3)+
  #geom_violin(alpha=0.3, trim=FALSE) +
  geom_point(shape=21, color="black", position=position_jitterdodge(jitter.width=0.2), size=1) +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              position=position_dodge(1)) +
  facet_grid(~taxonomy)+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr()+
  labs(x="condition",y="CLR abundance")+
  theme(legend.position = "none")
#ggsave("metagenomics/SFR24_1028_5limosig_clrabun_byfastedfed_cond_taxa.pdf", plot=p,height=3.5, width=10)
ggsave("metagenomics/SFR24_1028_5limosig_clrabun_bycond_fasted_taxa.pdf", plot=p,height=3.5, width=10)


limo5_pvals<-limo_LDplt%>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, fasted_fed_pval))%>%
  dplyr::select(-data)%>%
  unnest()

limo5_pvals2<-limo_LDplt%>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, fasted_fed_pval2))%>%
  dplyr::select(-data)%>%
  unnest()

limo_LDplt_notnorm<-df_species_notnorm %>%
  filter(taxonomy %in% sig_limo$taxonomy)%>%
  left_join(.,md,by="sample_name")

limo_LDpltratio<-limo_LDplt_notnorm%>% 
  group_by(taxonomy,mouseid,fasted)%>%
  reframe(mn_counts = mean(counts))%>%
  spread(fasted,mn_counts)%>%
  mutate(ratioff=fasted/fed,
         log_ratioff=log2(fasted/fed))%>%
  left_join(.,md_subject,by="mouseid")%>%
  mutate(log_ratioff=ifelse(log_ratioff==Inf,5,ifelse(log_ratioff==-Inf,-7,log_ratioff)))

p<-ggplot(limo_LDpltratio, aes(x=genotype_diet, y=log_ratioff, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) +
  #geom_violin(alpha=0.3, trim=FALSE) +
  #geom_point(shape=21, color = "black", position=position_jitterdodge())+
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1)) +
  facet_wrap(~taxonomy,scales="free_y",nrow=1)+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr()+
  labs(x="condition",y="log2(fasted/fed)")+
  theme(legend.position = "none")
ggsave("metagenomics/SFR24_1028_5limosig_log2ratiofastedfed.pdf", plot=p,height=3.5, width=10)

leq<-limo_LDpltratio%>%filter(taxonomy=="s__Limosilactobacillus equigenerosi")
pairwise.t.test(leq$log_ratioff, leq$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.93    -     
#   Tg.TRF 0.19    0.19
lre<-limo_LDpltratio%>%filter(taxonomy=="s__Limosilactobacillus reuteri")
pairwise.t.test(lre$log_ratioff, lre$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.87    -     
#   Tg.TRF 0.49    0.49
lfe<-limo_LDpltratio%>%filter(taxonomy=="s__Limosilactobacillus fermentum")
pairwise.t.test(lfe$log_ratioff, lfe$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.95    -     
#   Tg.TRF 0.95    0.95  
lmu<-limo_LDpltratio%>%filter(taxonomy=="s__Limosilactobacillus mucosae")
pairwise.t.test(lmu$log_ratioff, lmu$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.58    -     
#   Tg.TRF 0.15    0.20 
lse<-limo_LDpltratio%>%filter(taxonomy=="s__Limosilactobacillus secaliphilus")
pairwise.t.test(lse$log_ratioff, lse$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.038   -     
#   Tg.TRF 0.038   0.806

df_lan<-subset_dat(df_species, "s__Limosilactobacillus antri")

df_lan%>%group_by(genotype_diet,fasted)%>%
  summarise(n=n())

p<-fasted_fed_plt(df_lan,"Limosilactobacillus antri")
ggsave("metagenomics/SFR24_1023_fastedfed_Limosilactobacillusantri.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_lan)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.7508011 0.0003736382 0.0003736382 fasted
# 2    0.2775330 0.5600830496 0.0170086888    fed

df_lan_ratio<-get_ratio_dat(df_lan)
p<-fastedfed_ratio_plt(df_lan_ratio,"Limosilactobacillus antri")
ggsave("metagenomics/SFR24_1023_fastedfedratio_Limosilactobacillusantri_rmoutlier.pdf", plot=p,height=3.5, width=3)
pairwise.t.test(df_lan_ratio$ratioff, df_lan_ratio$genotype_diet,p.adjust.method="none")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.042   -     
#   Tg.TRF 0.510   0.175 

df_lsp<-subset_dat(df_species, "s__Limosilactobacillus sp900156885")
p<-fasted_fed_plt(df_lsp,"Limosilactobacillus sp900156885")
ggsave("metagenomics/SFR24_1023_fastedfed_Limosilactobacillussp900156885.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_lsp)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.5062283  0.002872502 0.0006584341 fasted
# 2    0.2502305  0.560083050 0.0216911905    fed

######################################################################################

all_prevotella<-df_species %>%
  filter(str_detect(taxonomy, "Prevotella"))%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, anova_allmicrobes))

prevotella_anova<-all_prevotella%>%dplyr::select(-data)%>%unnest()
write.table(prevotella_anova,"metagenomics/SFR24_1028_ADprevotella_sp_anova_condtime_pvals.txt",sep = "\t",row.names = FALSE,quote=FALSE)

sig_prevotella<-prevotella_anova%>%filter(genotype_diet<0.1 & zt_time<0.1)
sig_prevotella<-fread("metagenomics/ancombc/ancombbc_primaryresfastedtrf_B10.txt")%>%
  filter(`p_genotype_dietTg.TRF:fastedfed`<0.05 & `passed_ss_genotype_dietTg.TRF:fastedfed`==1)%>%
  filter(str_detect(taxon,"Prevotella"))%>%
  mutate(taxon = sub(".*; (s__.*)", "\\1", taxon)) 

prev_heatmap<-df_species %>%
#prev_heatmap<-df_species_notnorm %>%
  filter(taxonomy %in% sig_prevotella$taxon)%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr))%>%
  group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_clr - mean(mn_clr))/sd(mn_clr))%>%
  #dplyr::summarise(mn_counts=mean(counts))%>%
  #group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_counts - mean(mn_counts))/sd(mn_counts))%>%
  mutate(zt_time=factor(zt_time, levels = c("0", "4","8","12","16","20")))

plt<-ggplot(prev_heatmap,aes(x=zt_time, y=taxonomy)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~genotype_diet,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 
#ggsave("metagenomics/SFR24_1028_35prevsig_heatmap.pdf", plot=plt,height=7, width=7)
ggsave("metagenomics/ancombc/SFR25_0304_12prevsigfastedtrf_heatmap.pdf", plot=plt,height=3.5, width=7)
#ggsave("metagenomics/SFR24_1028_35prevsig_heatmap_usingnotnormcounts.pdf", plot=plt,height=7, width=7)

#ridgeplot
sem <- function(x) sd(x)/sqrt(length(x))
prev_ridge<-df_species %>%
  filter(taxonomy %in% sig_prevotella$taxon)%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  summarise(meanclr=mean(clr),
            semclr=sem(clr))

p <- ggplot(prev_ridge, aes(x=zt_time, y=meanclr,colour = genotype_diet)) +
  geom_point(alpha = 1.0) +
  geom_line() +
  geom_ribbon(aes(ymin = meanclr - semclr, ymax = meanclr + semclr, fill = genotype_diet), alpha = 0.3, color = NA) +
  facet_wrap(~taxonomy,scales="free_y",nrow=1) +
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20))
ggsave("metagenomics/ancombc/prevsigwtimefastedtrf_clrplt.pdf", plot=p,height=2.5, width=30)

p<-ggplot(prev_ridge, aes(x = zt_time, y = genotype_diet, fill = genotype_diet)) +
  geom_density_ridges(y=prev_ridge$clr, alpha=0.3)+
  facet_grid(fasted~taxonomy, scales="free_x")+
  theme_bw()+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme(legend.position = "bottom")
ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_bycond_fastedfed_taxa.ridgeplot.pdf", plot=p,height=2.5, width=10)
ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_bycond_fastedfed_taxa.ridgeplot_ext.pdf", plot=p,height=3, width=10)


sem <- function(x) sd(x)/sqrt(length(x))

prev_lnplt<-df_species %>%
  filter(taxonomy %in% sig_prevotella$taxonomy)%>%
  left_join(.,md,by="sample_name")%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  dplyr::summarise(mn_clr=mean(clr),
                   sem=sem(clr))

p <- ggplot(prev_lnplt, aes(x = zt_time, y = mn_clr, color = genotype_diet)) +
  geom_point(alpha = 1.0) +
  geom_line() +
  geom_ribbon(aes(ymin = mn_clr - sem, ymax = mn_clr + sem, fill = genotype_diet), alpha = 0.3, color = NA) +
  facet_wrap(~taxonomy) +
  theme_pubr() +
  scale_color_manual(values = c("darkorange", "purple4", "forestgreen")) +
  scale_fill_manual(values = c("darkorange", "purple4", "forestgreen")) +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20))
ggsave("metagenomics/SFR24_1028_35prevsig_facelinplt_condi.pdf", plot=p,height=10, width=10)

p <- ggplot(prev_lnplt, aes(x = zt_time, y = mn_clr, color = taxonomy)) +
  geom_point(alpha = 1.0) +
  geom_line() +
  geom_ribbon(aes(ymin = mn_clr - sem, ymax = mn_clr + sem, fill = taxonomy), alpha = 0.3, color = NA) +
  facet_grid(~genotype_diet) +
  theme_pubr() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20))
ggsave("metagenomics/SFR24_1028_35prevsig_facelinplt_taxa.pdf", plot=p,height=4, width=10)

prev_LDplt<-df_species %>%
  filter(taxonomy %in% sig_prevotella$taxonomy)%>%
  left_join(.,md,by="sample_name")

p<-ggplot(prev_LDplt, aes(x=genotype_diet, y=clr, fill=fasted)) +
  geom_boxplot(alpha=0.3)+
  #geom_violin(alpha=0.3, trim=FALSE) +
  geom_point(shape=21, color="black", position=position_jitterdodge(jitter.width=0.2), size=1) +
  # geom_dotplot(binaxis='y', stackdir='center',
  #              position=position_dodge(1)) +
  facet_wrap(~taxonomy, scales="free_y")+
  # scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  # scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr()+
  labs(x="condition",y="CLR abundance")+
  theme(legend.position = "none")
ggsave("metagenomics/SFR24_1028_35prevsig_clrabun_byfastedfed_cond_taxa.pdf", plot=p,height=10, width=12)
#ggsave("metagenomics/SFR24_1028_35prevsig_clrabun_bycond_fasted_taxa.pdf", plot=p,height=10, width=12)

prev35_pvals<-prev_LDplt%>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, fasted_fed_pval))%>%
  dplyr::select(-data)%>%
  unnest()

prev35_pvals2<-prev_LDplt%>%
  group_by(taxonomy)%>%
  nest()%>%
  mutate(pvals = map(data, fasted_fed_pval2))%>%
  dplyr::select(-data)%>%
  unnest()

#just sig with fasted/fed diff 

sig_heatmap<-df_species_notnorm %>%
  filter(taxonomy %in% c("s__Limosilactobacillus secaliphilus",
                         "s__Limosilactobacillus reuteri",
                         "s__Limosilactobacillus mucosae",
                         "s__Prevotella baroniae",
                         "s__Prevotella sp002251295",
                         "s__Prevotella sp003447235",
                         "s__Prevotella sp004792655",
                         "s__Prevotella sp900318625"))%>%
  left_join(.,md,by="sample_name")%>%
  #group_by(taxonomy,genotype_diet,fasted)%>%
  group_by(taxonomy,genotype_diet,zt_time)%>%
  #dplyr::summarise(mn_clr=mean(clr))%>%
  #group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_clr - mean(mn_clr))/sd(mn_clr))%>%
  dplyr::summarise(mn_counts=mean(counts))%>%
  group_by(taxonomy,genotype_diet)%>%mutate(Zscore=(mn_counts - mean(mn_counts))/sd(mn_counts))#%>%
  mutate(zt_time=factor(zt_time, levels = c("0", "4","8","12","16","20")))

#plt<-ggplot(sig_heatmap,aes(x=fasted, y=taxonomy)) +theme_classic()+
plt<-ggplot(sig_heatmap,aes(x=zt_time, y=taxonomy)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+
  facet_grid(~genotype_diet,scales="free",space="free")+
  theme(axis.ticks.y=element_blank(),strip.text.y = element_text(angle = 0))+
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") 
ggsave("metagenomics/SFR24_1030_sigspec_heatmap.pdf", plot=plt,height=1.8, width=7)

sig_LDplt<-df_species %>%
  filter(taxonomy %in% c("s__Limosilactobacillus secaliphilus",
                         "s__Limosilactobacillus reuteri",
                         "s__Limosilactobacillus mucosae",
                         "s__Prevotella baroniae",
                         "s__Prevotella sp002251295",
                         "s__Prevotella sp003447235",
                         "s__Prevotella sp004792655",
                         "s__Prevotella sp900318625"))%>%
  left_join(.,md,by="sample_name")

p<-ggplot(sig_LDplt, aes(x=genotype_diet, y=clr, fill=genotype_diet)) +
  #geom_boxplot(alpha=0.3)+
  geom_violin(alpha=0.3, trim=FALSE) +
  #geom_point(shape=21, color="black", position=position_jitterdodge(jitter.width=0.2), size=1) +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1)) +
  facet_grid(taxonomy~fasted, scales="free_y")+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_bw()+
  labs(x="condition",y="CLR abundance")+
  theme(legend.position = "none")
#ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_byfastedfed_cond_taxa.pdf", plot=p,height=8, width=4.5)
ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_byfastedfed_cond_taxa.violin.pdf", plot=p,height=8, width=4.5)

#ridgeplot?
p<-ggplot(sig_LDplt, aes(x = clr, y = genotype_diet, fill = genotype_diet)) +
  geom_density_ridges(alpha=0.3)+
  facet_grid(fasted~taxonomy, scales="free_x")+
  theme_bw()+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme(legend.position = "bottom")
ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_bycond_fastedfed_taxa.ridgeplot.pdf", plot=p,height=2.5, width=10)
ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_bycond_fastedfed_taxa.ridgeplot_ext.pdf", plot=p,height=3, width=10)


p<-ggplot(sig_LDplt, aes(x = clr, y = genotype_diet, fill = genotype_diet)) +
  geom_density_ridges(alpha=0.3)+
  facet_grid(fasted~taxonomy, scales="free_x")+
  theme_bw()+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))
ggsave("metagenomics/SFR24_1028_sigmicrobs_clrabun_bycond_fastedfed_taxa.ridgeplot.pdf", plot=p,height=2.5, width=10)


df_prevsp<-subset_dat(df_species, "s__Prevotella sp002372235")
p<-fasted_fed_plt(df_prevsp,"Prevotella sp002372235")
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotellasp002372235.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prevsp)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.4349660 8.042873e-05 0.0001587145 fasted
# 2    0.3155795 5.467207e-03 0.0002472083    fed

df_prevsp_ratio<-get_ratio_dat(df_prevsp)
p<-fastedfed_ratio_plt(df_prevsp_ratio,"Prevotella sp002372235")
ggsave("metagenomics/SFR24_1023_fastedfedratio_Prevotellasp002372235.pdf", plot=p,height=3.5, width=3)
pairwise.t.test(df_prevsp_ratio$ratioff, df_prevsp_ratio$genotype_diet,p.adjust.method="none")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.45    -     
#   Tg.TRF 0.83    0.59  

df_prevloe<-subset_dat(df_species, "s__Prevotella loescheii")
p<-fasted_fed_plt(df_prevloe,"Prevotella loescheii")
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotellaloescheii.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prevloe)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.4125411 0.0004618167 4.618167e-04 fasted
# 2    0.8006307 0.0002841360 1.036849e-05    fed

df_prevnig<-subset_dat(df_species, "s__Prevotella nigrescens")
p<-fasted_fed_plt(df_prevnig,"Prevotella nigrescens")
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotellanigrescens.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prevnig)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.5699121 0.0004618167 3.928062e-04 fasted
# 2    0.5807963 0.0009104129 2.522188e-05    fed

df_prevnig<-subset_dat(df_species, "s__Prevotella sp000688375")
p<-fasted_fed_plt(df_prevnig,"Prevotella sp000688375")
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotellasp000688375.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prevnig)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.6931098 0.0008312016 0.0009506346 fasted
# 2    0.6067473 0.0017618516 0.0013583017    fed

df_prevnig<-subset_dat(df_species, "s__Prevotella ruminicola")
p<-fasted_fed_plt(df_prevnig,"Prevotella ruminicola")
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotellaruminicola.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prevnig)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.5186645 0.0001723139 0.0001723139 fasted
# 2    0.8442984 0.0030384119 0.0005969445    fed

df_prevnig<-subset_dat(df_species, "s__Prevotella sp900119285")
p<-fasted_fed_plt(df_prevnig,"Prevotella sp900119285")
ggsave("metagenomics/SFR24_1023_fastedfed_Prevotellasp900119285.pdf", plot=p,height=3.5, width=4.5)
pvals<-fasted_fed_pval(df_prevnig)
# NtgALFvTgALF NtgALFvTgTRF  TgALFvTgTRF fasted
# 1    0.1491350 1.342459e-05 0.0051472478 fasted
# 2    0.3247643 2.514543e-04 0.0002514543    fed

df_prevnig_ratio<-get_ratio_dat(df_prevnig)
p<-fastedfed_ratio_plt(df_prevnig_ratio,"Prevotella sp900119285")
ggsave("metagenomics/SFR24_1023_fastedfedratio_Prevotellasp900119285r.pdf", plot=p,height=3.5, width=3)
pairwise.t.test(df_prevnig_ratio$ratioff, df_prevnig_ratio$genotype_diet,p.adjust.method="none")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.48    -     
#   Tg.TRF 0.62    0.23  
