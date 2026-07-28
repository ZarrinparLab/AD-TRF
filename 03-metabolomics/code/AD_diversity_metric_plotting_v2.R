setwd("~/scratch/AD_EM/metabolomics/")

library("tidyverse")
library("data.table")
library("qiime2R")
library("ggpubr")
library("VennDiagram")
library("ape")
library("ggsci")
library("phyloseq")
library("vegan")

##################################################
#just feces
md<-fread("data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces")%>%
  mutate(ATTRIBUTE_sex=ifelse(ATTRIBUTE_sex=="F","female",ATTRIBUTE_sex),
         genofedZT=paste(ATTRIBUTE_ZT_collection_time, ATTRIBUTE_geno_n_feedingcond, sep = "_"))
write.table(md,"data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-forrpca.txt",sep = "\t",row.names = FALSE,quote=FALSE)

mtb<-fread("data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean.csv")%>%
  dplyr::rename(FeatureID=row.ID)%>%
  select(FeatureID, all_of(md$sampleid))
write.table(mtb,"data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt",sep = "\t",row.names = FALSE,quote=FALSE)

##################################################

#plot the alpha diversity

shannon<-fread("diversity_metrics/alpha_shannon_fecal_all/alpha/alpha-diversity.tsv")%>% 
  rename(sampleid=V1) %>%
  left_join(.,md, by="sampleid")%>%
  mutate(ATTRIBUTE_ZT_collection_time=factor(ATTRIBUTE_ZT_collection_time,c("ZT0","ZT4","ZT8","ZT12","ZT16","ZT20")))

p<-ggplot(shannon, aes(x=ATTRIBUTE_sex, y=shannon_entropy,fill=ATTRIBUTE_sex)) + geom_violin(alpha=0.3, trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  scale_fill_manual(values=c("blue","red"))+
  labs(y="shannon estimate",x="")+theme_minimal()+ theme(legend.position = "top")
ggsave("diversity_metrics/alpha_shannon_fecal_all/shannon_violinplot_sex.pdf", plot=p,height=3.5, width=5)
pairwise.wilcox.test(shannon$shannon_entropy,shannon$ATTRIBUTE_sex,p.adjust.method = "fdr")
# female
# male 0.083 
p<-ggplot(shannon, aes(x=ATTRIBUTE_geno_n_feedingcond, y=shannon_entropy,fill=ATTRIBUTE_geno_n_feedingcond)) + geom_violin(alpha=0.3, trim=FALSE) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  labs(y="shannon estimate",x="")+theme_minimal()+ theme(legend.position = "top")
ggsave("diversity_metrics/alpha_shannon_fecal_all/shannon_violinplot_genofeedcond.pdf", plot=p,height=3.5, width=5)
pairwise.wilcox.test(shannon$shannon_entropy,shannon$ATTRIBUTE_geno_n_feedingcond,p.adjust.method = "fdr")
# NTg-ALF Tg-ALF
# Tg-ALF 0.85    -     
#   Tg-TRF 0.85    0.85 
p <- ggplot(shannon, aes(x=ATTRIBUTE_ZT_collection_time, y=shannon_entropy, fill=ATTRIBUTE_geno_n_feedingcond)) + 
  geom_boxplot(alpha=0.3, position=position_dodge(width=0.75)) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75), size=1.5, aes(color=ATTRIBUTE_geno_n_feedingcond)) +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen")) +
  scale_color_manual(values=c("darkorange","purple4","forestgreen")) +
  labs(y="shannon estimate", x="") + 
  theme_minimal() + 
  theme(legend.position = "top")
ggsave("diversity_metrics/alpha_shannon_fecal_all/shannon_boxplot_ZT_fclgroup.pdf", plot=p,height=3.5, width=8)

##just fecal shannon
shannonTRF<-shannon%>%filter(ATTRIBUTE_geno_n_feedingcond=="Tg-TRF")
pairwise.wilcox.test(shannonTRF$shannon_entropy,shannonTRF$ATTRIBUTE_ZT_collection_time,p.adjust.method = "fdr")
# ZT0    ZT4    ZT8    ZT12   ZT16  
# ZT4  0.4476 -      -      -      -     
#   ZT8  0.2667 0.5000 -      -      -     
#   ZT12 0.1127 0.1748 0.9890 -      -     
#   ZT16 0.0521 0.1055 1.0000 1.0000 -     
#   ZT20 0.2173 0.1392 0.1127 0.0233 0.0093

shannonALF<-shannon%>%filter(ATTRIBUTE_geno_n_feedingcond=="Tg-ALF")
pairwise.wilcox.test(shannonALF$shannon_entropy,shannonALF$ATTRIBUTE_ZT_collection_time,p.adjust.method = "fdr")
# ZT0  ZT4  ZT8  ZT12 ZT16
# ZT4  0.92 -    -    -    -   
#   ZT8  0.92 1.00 -    -    -   
#   ZT12 1.00 0.92 0.92 -    -   
#   ZT16 0.92 0.92 0.92 0.92 -   
#   ZT20 0.92 0.92 0.92 0.92 0.92

shannonNTg<-shannon%>%filter(ATTRIBUTE_geno_n_feedingcond=="NTg-ALF")
pairwise.wilcox.test(shannonNTg$shannon_entropy,shannonNTg$ATTRIBUTE_ZT_collection_time,p.adjust.method = "fdr")
# ZT0  ZT4  ZT8  ZT12 ZT16
# ZT4  0.28 -    -    -    -   
#   ZT8  0.10 1.00 -    -    -   
#   ZT12 0.42 0.67 0.67 -    -   
#   ZT16 1.00 0.43 0.47 0.74 -   
#   ZT20 0.42 0.14 0.10 0.14 0.49
##################################################
#RPCA

set.seed(1234)
ord <- read_qza("diversity_metrics/rpca_fecal_all/ordination.qza")
dis_mat<-fread("diversity_metrics/rpca_fecal_all/distance_matrix/distance-matrix.tsv")%>%column_to_rownames("V1")

rpca_fcl<-ord$data$Vectors %>%
  select(SampleID, PC1, PC2,PC3) %>%
  dplyr::rename(sampleid=SampleID)%>%
  left_join(md,by="sampleid")%>%
  mutate(ATTRIBUTE_ZT_collection_time=factor(ATTRIBUTE_ZT_collection_time,c("ZT0","ZT4","ZT8","ZT12","ZT16","ZT20")))

p<-rpca_fcl %>%
  ggplot(aes(x=PC1, y=PC2, fill=BiologicalSex)) +
  geom_point(alpha=1.0, shape=21, size=2) +
  theme_pubr() +
  scale_fill_manual(values=c("#CC0000","blue3"))+
  scale_color_manual(values=c("#CC0000","blue3"))+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="conditions")+ theme(plot.title = element_text(face = "bold"))

ggsave("diversity_metrics/rpca_fecal_all/SFR26_0708_rpca_fecal_ADsex.pdf", plot=p,height=3.5, width=3.5)
permanova_result <- adonis2(dis_mat ~ BiologicalSex, data = md , permutations = 999) #0.253

p<-rpca_fcl %>%
  ggplot(aes(x=PC1, y=PC2, fill=ATTRIBUTE_geno_n_feedingcond)) +
  geom_point(alpha=1.0, shape=21, size=2) +
  theme_pubr() +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen")) +
  scale_color_manual(values=c("darkorange","purple4","forestgreen")) +
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="feeding and genotype groups")
ggsave("diversity_metrics/rpca_fecal_all/SFR26_0708_RPCA_fecal_genofeedinggrp.pdf", plot=p,height=3.5, width=3.5)

means <- rpca_fcl %>%
  group_by(ATTRIBUTE_ZT_collection_time) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

p<-rpca_fcl %>%
  ggplot(aes(x=PC1, y=PC2, fill=ATTRIBUTE_ZT_collection_time, shape=ATTRIBUTE_geno_n_feedingcond)) +
  geom_point(alpha=1.0,size=3) +
  theme_pubr() +
  scale_shape_manual(values = c(21, 22, 23)) + 
  scale_fill_brewer(palette = "YlGnBu")+theme(legend.position = "right")+
  scale_color_brewer(palette = "YlGnBu")+theme(legend.position = "right")+
  labs(x =paste("PC1 (",round(ord$data$ProportionExplained$PC1*100,digits=2),"%)",sep=""),
       y =paste("PC2 (",round(ord$data$ProportionExplained$PC2*100,digits=2),"%)",sep=""),
       color="ZT time") +
  geom_point(data = means, aes(x = PC1, y = PC2, color = ATTRIBUTE_ZT_collection_time),
               size = 4.5, shape = 8)
ggsave("diversity_metrics/rpca_fecal_all/SFR24_1101_RPCA_fecal_condZTtime_wmeanpt.pdf", plot=p,height=3.5, width=7)

#calc permanova
md_dis<-md%>%mutate(fasted=ifelse(ATTRIBUTE_ZT_collection_time %in% c("ZT8","ZT12","ZT16"),"fasted","fed"))
permanova_result <- adonis2(dis_mat ~ ATTRIBUTE_geno_n_feedingcond+ATTRIBUTE_ZT_collection_time, data = md_dis , permutations = 999,by = "terms")
print(permanova_result)
# adonis2(formula = dis_mat ~ ATTRIBUTE_geno_n_feedingcond + ATTRIBUTE_ZT_collection_time, data = md_dis, permutations = 999, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)    
# ATTRIBUTE_geno_n_feedingcond   2   37.234 0.13139 14.894  0.001 ***
#   ATTRIBUTE_ZT_collection_time   5   87.410 0.30844 13.986  0.001 ***
#   Residual                     127  158.749 0.56017                  
# Total                        134  283.392 1.0000

permanova_result <- adonis2(dis_mat ~ ATTRIBUTE_geno_n_feedingcond+fasted, data = md_dis , permutations = 999,by="terms")
print(permanova_result)
# adonis2(formula = dis_mat ~ ATTRIBUTE_geno_n_feedingcond + fasted, data = md_dis, permutations = 999, by = "terms")
# Df SumOfSqs      R2      F Pr(>F)    
# ATTRIBUTE_geno_n_feedingcond   2   37.234 0.13139 14.138  0.001 ***
#   fasted                         1   73.663 0.25993 55.943  0.001 ***
#   Residual                     131  172.495 0.60868                  
# Total                        134  283.392 1.00000 

#show PC1 over time

p<-rpca_fcl %>%
  ggplot(aes(x=ATTRIBUTE_ZT_collection_time, y=PC1, fill=ATTRIBUTE_geno_n_feedingcond)) +
  geom_boxplot(alpha=0.3)+
  theme_pubr() +
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(0.75),
               dotsize=0.75)
ggsave("diversity_metrics/rpca_fecal_all/SFR26_0708_RPCA_fecal_PC1vZTtime.pdf", plot=p,height=4, width=6)
