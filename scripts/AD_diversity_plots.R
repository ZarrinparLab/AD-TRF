setwd("~/scratch/AD_EM")

library(tidyverse)
library(data.table)
library(qiime2R)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(emmeans)
#####################################################
md <- fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  mutate(mouse_id = str_extract(sample_name, "^[^.]+"))

alphafpd <- read_qza("metagenomics/diversity_results/faith_pd_vector.qza")$data %>%
  dplyr::rename(sample_name=V1,faithpd=V2)%>%
  left_join(.,md,by="sample_name")

alphafpd%>%group_by(genotype_diet)%>%
  summarise(n=n())

p<-ggplot(alphafpd, aes(x=genotype_diet, y=faithpd, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr()+
  labs(x="condition",y="Faith's PD", title="alpha diversity")+
  theme(legend.position = "none")

ggsave("metagenomics/diversity_results/SFR24_1018_AD_alphafaithpd.pdf", plot=p,height=4, width=3.5)

# m <- lmer(faithpd ~ genotype_diet + zt_time+ (1|mouse_id), data = alphafpd)
# summary(m)
# emmeans(m, pairwise ~ genotype_diet, adjust = "fdr")$contrasts
# contrast         estimate   SE   df t.ratio p.value
# Ntg.ALF - Tg.ALF     4.33 8.25 21.8   0.525  0.6049
# Ntg.ALF - Tg.TRF   -12.17 8.52 22.1  -1.429  0.2507
# Tg.ALF - Tg.TRF    -16.50 8.28 22.1  -1.992  0.1765

pairwise.wilcox.test(alphafpd$faithpd, alphafpd$genotype_diet,p.adjust.method="fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.34976 -      
#   Tg.TRF 0.01725 0.00065

#####################################################
#RPCA
ord <- read_qza("metagenomics/diversity_results/wol2_genome_none.gemelli.biplot.qza") #check deicode.sh
md <- fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  mutate(mouse_id = str_extract(sample_name, "^[^.]+"))%>%
  dplyr::rename(SampleID=sample_name)%>%
  mutate(fasted=ifelse(zt_collection_time %in% c("ZT08","ZT12","ZT16"),"fasted","fed"))

rpcaAD<-ord$data$Vectors%>%
  dplyr::select(SampleID, PC1, PC2)%>%
  left_join(.,md,by="SampleID")

p<-rpcaAD %>%
  ggplot(aes(x=PC1, y=PC2, fill=genotype_diet)) +
  geom_point(alpha=1.0, shape=21, size=2) +
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(color = genotype_diet))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="conditions")+ theme(plot.title = element_text(face = "bold"))

ggsave("metagenomics/diversity_results/SFR24_1018_rpca_ADgenotypediet.pdf", plot=p,height=3.5, width=3.5)

cond_rpca <- ggplot(rpcaAD, aes(x =genotype_diet, y = PC1,fill=genotype_diet)) +
  geom_boxplot(alpha=0.9) + #geom_point(colour="black", aes(fill=genotype_diet),
                                       #position=position_jitter(width=0.1,height = 0.1),
                                       #size=1.5,pch=21) + 
  labs(x = NA, y = "PC1") +
  scale_y_continuous(breaks=c(-0.2,-0.1,0.,0.1,0.2,0.3))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  # theme_void()
  scale_x_discrete(limits=rev)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") + coord_flip()

# m <- lmer(PC1 ~ genotype_diet + zt_time+ (1|mouse_id), data = rpcaAD)
# summary(m)
# emmeans(m, pairwise ~ genotype_diet, adjust = "fdr")$contrasts
# contrast         estimate     SE   df t.ratio p.value
# Ntg.ALF - Tg.ALF  -0.0164 0.0228 21.7  -0.719  0.4800
# Ntg.ALF - Tg.TRF   0.0881 0.0236 22.1   3.738  0.0017
# Tg.ALF - Tg.TRF    0.1045 0.0229 22.3   4.554  0.0005

#pairwise.t.test(rpcaAD$PC1, rpcaAD$genotype_diet, p.adjust.method = "fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.25    -      
# Tg.TRF 2.3e-08 6.5e-11

pairwise.wilcox.test(rpcaAD$PC1, rpcaAD$genotype_diet, p.adjust.method = "fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.28    -      
#   Tg.TRF 1.8e-07 6.6e-09

final_plot <- plot_grid(p,cond_rpca,nrow=2,rel_heights  = c(4, 1))
ggsave("metagenomics/diversity_results/SFR24_1018_rpca_ADgenotypediet_addbxplt.pdf", plot=final_plot,height=5, width=4)

#to show no female male differences
p<-rpcaAD %>%
  ggplot(aes(x=PC1, y=PC2, fill=sex)) +
  geom_point(alpha=1.0, shape=21, size=2) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +#stat_ellipse(type = "t", linetype = 2,aes(color = genotype_diet))+
  scale_fill_manual(values=c("#CC0000","blue3"))+
  scale_color_manual(values=c("#CC0000","blue3"))+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="conditions")+ theme(plot.title = element_text(face = "bold"))

ggsave("metagenomics/diversity_results/SFR25_0128_rpca_ADsex.pdf", plot=p,height=3.5, width=3.5)

#need distance matrix
dist_mat<-fread("metagenomics/diversity_results/wol2_genome_none.gemelli.distance.tsv")%>%column_to_rownames("V1")
permanova_result <- adonis2(dist_mat ~ sex, data = md , permutations = 999, by="terms")
print(permanova_result)

# adonis2(formula = dist_mat ~ sex, data = md, permutations = 999)
# Df SumOfSqs      R2      F Pr(>F)
# sex        1    1.723 0.00591 0.8505  0.422
# Residual 143  289.708 0.99409              
# Total    144  291.431 1.00000 

#to show time effect
means <- rpcaAD %>%
  mutate(zt_time=factor(zt_time))%>%
  group_by(zt_time) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

p<-rpcaAD %>%
  ggplot(aes(x=PC1, y=PC2, fill=factor(zt_time), shape=genotype_diet)) +
  geom_point(alpha=1.0,size=2) +
  theme_pubr() +
  #stat_ellipse(type = "t", linetype = 2,aes(color = ATTRIBUTE_ZT_collection_time))+
  scale_shape_manual(values = c(21, 22, 23)) + 
  scale_fill_brewer(palette = "RdYlBu")+theme(legend.position = "right")+
  scale_color_brewer(palette = "RdYlBu")+theme(legend.position = "right")+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="ZT time") +
  geom_point(data = means, aes(x = PC1, y = PC2, color = zt_time),
             size = 4.5, shape = 8)

dist_mat<-fread("metagenomics/diversity_results/wol2_genome_none.gemelli.distance.tsv")%>%column_to_rownames("V1")
permanova_result <- adonis2(dist_mat ~ zt_time, data = md , permutations = 999, by="terms")
print(permanova_result)
# adonis2(formula = dist_mat ~ zt_time, data = md, permutations = 999)
# Df SumOfSqs      R2      F Pr(>F)
# zt_time    1    4.736 0.01625 2.3625  0.1
# Residual 143  286.695 0.98375              
# Total    144  291.431 1.00000 
ggsave("metagenomics/diversity_results/SFR25_0304_rpca_ADzt.pdf", plot=p,height=3.5, width=5.5)

permanova_result <- adonis2(dist_mat ~ zt_time*genotype_diet, data = md , permutations = 999, by="terms")
# adonis2(formula = dist_mat ~ zt_time * genotype_diet, data = md, permutations = 999)
# Df SumOfSqs      R2      F Pr(>F)
# zt_time                 1    4.736 0.01625 2.3612  0.116
# genotype_diet           2    5.052 0.01734 1.2593  0.301
# zt_time:genotype_diet   2    2.811 0.00964 0.7006  0.597
# Residual              139  278.832 0.95677              
# Total                 144  291.431 1.00000

permanova_result <- adonis2(dist_mat ~ fasted*genotype_diet, data = md , permutations = 999, by="terms")
# adonis2(formula = dist_mat ~ fasted * genotype_diet, data = md, permutations = 999, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)
# fasted                 1    0.923 0.00317 0.4560  0.650
# genotype_diet          2    4.929 0.01691 1.2173  0.311
# fasted:genotype_diet   2    4.200 0.01441 1.0374  0.412
# Residual             139  281.379 0.96551              
# Total                144  291.431 1.00000 
#####################################################
#RPCA-ntg-alf and tg-alf
ord <- read_qza("metagenomics/diversity_results/wol2_genome_none.alf_only.biplot.qza") 
md <- fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  dplyr::rename(SampleID=sample_name)

rpcaAD<-ord$data$Vectors %>%
  select(SampleID, PC1, PC2)%>%
  left_join(.,md,by="SampleID")

p<-rpcaAD %>%
  ggplot(aes(x=PC1, y=PC2, fill=genotype_diet)) +
  geom_point(alpha=1.0, shape=21, size=2) + #alpha controls transparency and helps when points are overlapping
  theme_pubr() +stat_ellipse(type = "t", linetype = 2,aes(color = genotype_diet))+
  scale_fill_manual(values=c("darkorange","purple4"))+
  scale_color_manual(values=c("darkorange","purple4"))+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="conditions")+ theme(plot.title = element_text(face = "bold"))

ggsave("metagenomics/diversity_results/SFR24_1018_rpca_ADgenotypediet_ntgalftgalf.pdf", plot=p,height=3.5, width=3.5)
ggsave("metagenomics/diversity_results/SFR24_1018_rpca_ADgenotypediet_ntgalftgalf_wellipse.pdf", plot=p,height=3.5, width=3.5)

#####################################################
#Prevotella/Limosilactobacillus

md <- fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  mutate(mouse_id = str_extract(sample_name, "^[^.]+"))%>%
  dplyr::rename(SampleID=sample_name)%>%
  mutate(fasted=ifelse(zt_collection_time %in% c("ZT08","ZT12","ZT16"),"fasted","fed"))

natlog_prevlimo<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_prevotella_limosilactobacillus.tsv")%>%dplyr::select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")))

natlog_prevlimo%>%group_by(genotype_diet)%>%
  summarise(n=n())

prevlimo<-ggplot(natlog_prevlimo, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nPrevotella/Limosilactobacillus")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR24_1017_deicode_prevlimo_natlog.pdf",prevlimo, height=3, width=3.5)

natlog_prevlimo<-natlog_prevlimo%>%
  dplyr::rename(SampleID= `Sample ID`)%>%
  left_join(.,md,by="SampleID")

m <- lmer(Current_Natural_Log_Ratio ~ genotype_diet.x + zt_time+ (1|mouse_id), data = natlog_prevlimo)
summary(m)
emmeans(m, pairwise ~ genotype_diet.x, adjust = "fdr")$contrasts
# contrast         estimate    SE   df t.ratio p.value
# Ntg.ALF - Tg.ALF   -0.721 0.527 21.8  -1.367  0.1855
# Ntg.ALF - Tg.TRF    1.158 0.544 22.1   2.128  0.0671
# Tg.ALF - Tg.TRF     1.879 0.529 22.2   3.548  0.0054

pairwise.wilcox.test(natlog_prevlimo$Current_Natural_Log_Ratio, natlog_prevlimo$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF
# Tg.ALF 0.02212 -     
# Tg.TRF 0.00057 2e-09 

#Prevotella/Ligilactobacillus
natlog_prevligi<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_ligiprev.tsv")%>%select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")))

natlog_prevligi%>%group_by(genotype_diet)%>%
  summarise(n=n())

prevligi<-ggplot(natlog_prevligi, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nPrevotella/Ligilactobacillus")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR25_0304_deicode_prevligi_natlog.pdf",prevligi, height=3, width=3.5)

pairwise.wilcox.test(natlog_prevligi$Current_Natural_Log_Ratio, natlog_prevligi$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.14    -      
#   Tg.TRF 3.5e-06 1.4e-09 

#Prevotella/Ligilactobacillus & Limosilactobacillus
natlog_prevligilimo<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_ligilimoprev.tsv")%>%select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")))

natlog_prevligilimo%>%group_by(genotype_diet)%>%
  summarise(n=n())

prevligilimo<-ggplot(natlog_prevligilimo, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nPrevotella/(Ligilactobacillus + Limosilactobacillus)")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR25_0304_deicode_prevligilimo_natlog.pdf",prevligilimo, height=3, width=3.5)

pairwise.wilcox.test(natlog_prevligilimo$Current_Natural_Log_Ratio, natlog_prevligilimo$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.1     -      
#   Tg.TRF 3.1e-06 1.5e-09

#Prevotella/Enterococcus
natlog_preventero<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_prevotella_enterococcus.tsv")%>%select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")))

natlog_preventero%>%group_by(genotype_diet)%>%
  summarise(n=n())

preventero<-ggplot(natlog_preventero, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nPrevotella/Enterococcus")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR24_1017_deicode_preventero_natlog.pdf",preventero, height=3, width=3.5)

pairwise.wilcox.test(natlog_preventero$Current_Natural_Log_Ratio, natlog_preventero$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.59817 -      
# Tg.TRF 0.00280 0.00079

#Limosilactobacillus/Enterococcus
natlog_limoentero<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_limosilactobacillus_enterococcus.tsv")%>%select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")))

natlog_limoentero%>%group_by(genotype_diet)%>%
  summarise(n=n())

limoentero<-ggplot(natlog_limoentero, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nLimosilactobacillus/Enterococcus")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR24_1017_deicode_limoentero_natlog.pdf",limoentero, height=3, width=3.5)

pairwise.wilcox.test(natlog_limoentero$Current_Natural_Log_Ratio, natlog_limoentero$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.07466 -      
# Tg.TRF 0.20229 0.00021

#Bacillus/Limosilactobacillus

natlog_bacilimo<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_bacillus_limosilactobacillus.tsv")%>%select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

bacilimo<-ggplot(natlog_bacilimo, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nBacillus/Limosilactobacillus")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR24_1017_deicode_bacilimo_natlog.pdf",bacilimo, height=3, width=3.5)

pairwise.wilcox.test(natlog_bacilimo$Current_Natural_Log_Ratio, natlog_bacilimo$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF
# Tg.ALF 0.82    -     
# Tg.TRF 0.82    0.82 

#Bacillus/Enterococcus
natlog_bacientero<-fread("metagenomics/rpca_results_all/qurro_results/sample_plot_data_bacillus_enteroccus.tsv")%>%select(1:3)%>%
  mutate(genotype_diet=factor(genotype_diet, levels = c("Ntg.ALF", "Tg.ALF","Tg.TRF")),
         Current_Natural_Log_Ratio=as.numeric(Current_Natural_Log_Ratio))

bacientero<-ggplot(natlog_bacientero, aes(x=genotype_diet, y=Current_Natural_Log_Ratio,fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1)) +
  theme_minimal()+ scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(x="condition",y="Natural Log Ratio\nBacillus/Enterococcus")+
  theme(legend.position = "none")

ggsave("metagenomics/rpca_results_all/qurro_results//SFR24_1017_deicode_bacientero_natlog.pdf", bacientero, height=3, width=3.5)

pairwise.wilcox.test(natlog_bacientero$Current_Natural_Log_Ratio, natlog_bacientero$genotype_die,
                     p.adjust.method="fdr")

# Ntg.ALF Tg.ALF 
# Tg.ALF 0.096   -      
# Tg.TRF 0.057   7.9e-05

#plot ntg-alf vs tg-alf example
md<-fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  mutate(mouse_id = str_extract(sample_name, "^[^.]+"))%>%
  mutate(fasted=ifelse(zt_collection_time %in% c("ZT08","ZT12","ZT16"),"fasted","fed"))

micro_count<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(taxon=`OTU ID`)%>%
  dplyr::select(taxon,any_of(md$sample_name))%>%
  gather(sample_name,clr,-taxon)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(memdx=factor(genotype_diet, levels=c("Ntg.ALF","Tg.ALF","Tg.TRF")))

tsing<-micro_count%>%filter(str_detect(taxon, "s__Turicibacter sanguinis"))%>%
  filter(genotype_diet!="Tg.TRF")

p<-ggplot(tsing, aes(x=genotype_diet, y=clr, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 0.7) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="condition",y="CLR abundance", title="Turicibacter sanguinis")+
  theme(legend.position = "none")

ggsave("metagenomics/diversity_results/SFR25_0128_Turicibactersang_clrabun.pdf", plot=p,height=3.5, width=3.2)

m <- lmer(clr ~ genotype_diet + zt_time+ (1|mouse_id), data = tsing)
summary(m)
emmeans(m, pairwise ~ genotype_diet, adjust = "fdr")$contrasts
# contrast         estimate    SE df t.ratio p.value
# Ntg.ALF - Tg.ALF      2.6 0.987 15   2.633  0.0188
#pairwise.wilcox.test(tsing$clr,tsing$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF
# Tg.ALF 7.2e-05

tsp<-micro_count%>%filter(str_detect(taxon, "s__Turicibacter sp001543345"))%>%
  filter(genotype_diet!="Tg.TRF")

p<-ggplot(tsp, aes(x=genotype_diet, y=clr, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 0.7) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="condition",y="CLR abundance", title="Turicibacter sp001543345")+
  theme(legend.position = "none")

ggsave("metagenomics/diversity_results/SFR25_0128_Turicibactersp345_clrabun.pdf", plot=p,height=3.5, width=3.2)
m <- lmer(clr ~ genotype_diet + zt_time+ (1|mouse_id), data = tsp)
summary(m)
emmeans(m, pairwise ~ genotype_diet, adjust = "fdr")$contrasts
# contrast         estimate    SE df t.ratio p.value
# Ntg.ALF - Tg.ALF     2.93 0.825 15   3.549  0.0029
#pairwise.wilcox.test(tsp$clr,tsp$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF
# Tg.ALF 7.1e-10

#plot tg-alf vs tg-trf example
md<-fread("metagenomics/data_files/overly_specified_metadata.tsv")%>%
  mutate(mouse_id = str_extract(sample_name, "^[^.]+"))%>%
  mutate(fasted=ifelse(zt_collection_time %in% c("ZT08","ZT12","ZT16"),"fasted","fed"))

micro_count<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.genus.clr.csv")%>%
  dplyr::rename(taxon=`OTU ID`)%>%
  dplyr::select(taxon,any_of(md$sample_name))%>%
  gather(sample_name,clr,-taxon)%>%
  left_join(.,md,by="sample_name")%>%
  mutate(memdx=factor(genotype_diet, levels=c("Ntg.ALF","Tg.ALF","Tg.TRF")))

prevotella<-micro_count%>%filter(str_detect(taxon, "g__Prevotella")) #%>%
#   filter(genotype_diet!="Ntg.ALF")

p<-ggplot(prevotella, aes(x=genotype_diet, y=clr, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 0.7) +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  #scale_fill_manual(values=c("purple4","forestgreen"))+
  theme_pubr()+
  labs(x="condition",y="CLR abundance", title="Prevotella")+
  theme(legend.position = "none")

ggsave("metagenomics/diversity_results/SFR25_0226_Prevotella_clrabun.pdf", plot=p,height=2.8, width=3.6)
m <- lmer(clr ~ genotype_diet + zt_time+ (1|mouse_id), data = prevotella)
summary(m)
emmeans(m, pairwise ~ genotype_diet, adjust = "fdr")$contrasts
# contrast         estimate    SE   df t.ratio p.value
# Ntg.ALF - Tg.ALF   -0.140 0.150 21.8  -0.934  0.4195
# Ntg.ALF - Tg.TRF    0.127 0.154 22.1   0.823  0.4195
# Tg.ALF - Tg.TRF     0.267 0.150 22.2   1.775  0.2687
#pairwise.wilcox.test(prevotella$clr,prevotella$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF
# Tg.ALF 0.0835  -     
#   Tg.TRF 0.0835  0.0036

limo<-micro_count%>%filter(str_detect(taxon, "g__Limosilactobacillus"))#%>%
  #filter(genotype_diet!="Ntg.ALF")

p<-ggplot(limo, aes(x=genotype_diet, y=clr, fill=genotype_diet)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 0.7) +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  #scale_fill_manual(values=c("purple4","forestgreen"))+
  theme_pubr()+
  labs(x="condition",y="CLR abundance", title="Limosilactobacillus")+
  theme(legend.position = "none")

ggsave("metagenomics/diversity_results/SFR25_0226_Limosilactobacillus_clrabun.pdf", plot=p,height=2.8, width=3.6)
m <- lmer(clr ~ genotype_diet + zt_time+ (1|mouse_id), data = limo)
summary(m)
emmeans(m, pairwise ~ genotype_diet, adjust = "fdr")$contrasts
# contrast         estimate    SE   df t.ratio p.value
# Ntg.ALF - Tg.ALF    0.581 0.466 21.8   1.247  0.2258
# Ntg.ALF - Tg.TRF   -1.031 0.481 22.1  -2.143  0.0651
# Tg.ALF - Tg.TRF    -1.612 0.468 22.2  -3.443  0.0069
pairwise.wilcox.test(limo$clr,limo$genotype_diet,p.adjust.method="fdr")
# Ntg.ALF Tg.ALF 
# Tg.ALF 0.0341  -      
#   Tg.TRF 0.0015  2.3e-09

x <- c("darkorange","purple4","forestgreen")
col2hex <- function(x, alpha = FALSE) {
  args <- as.data.frame(t(col2rgb(x, alpha = alpha)))
  args <- c(args, list(names = x, maxColorValue = 255))
  do.call(rgb, args)
}
col2hex(x)
