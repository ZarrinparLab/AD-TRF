setwd("~/scratch/AD_EM")

library(tidyverse)
library(data.table)
library(ggpubr)
#####################################################
#aldex2 results

#taxon<-fread("metagenomics/data_files/lineages.split.tsv")
taxon<-fread("metagenomics/data_files/lineages.txt")%>%
  dplyr::rename(Featureid=featureid)
#res_aldex<-fread("metagenomics/human_BEAM/aldex2_preb.nc.mci.csv")%>%
res_aldex<-fread("metagenomics/human_BEAM/aldex2_prea.nc.mci.csv")%>%
  dplyr::rename(Featureid=V1)%>%
  left_join(.,taxon,by="Featureid")%>%
  dplyr::select(taxonomy,everything())%>%
  filter()
#write.table(res_aldex, "metagenomics/human_BEAM/aldex2_preb.nc.mci.clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(res_aldex, "metagenomics/human_BEAM/aldex2_prea.nc.mci.clean.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

significant_lfc <- res_aldex %>%
  #filter(abs(effect) >= 0.5, abs(diff.btw) >= 1) %>%
  arrange(diff.btw)


md<-fread("metagenomics/human_BEAM/qiita_sample_metadata.preab.tsv")%>%
  filter(visit_type=="prea-cog")%>%
  #filter(visit_type=="preb-cog")#%>%
  filter(!(sample_name %in% c("s_13662.5556874", "s_13662.5556880", "s_13662.5556877", "s_13662.5556892")))

counts_tab<-fread("metagenomics/human_BEAM/199894_none.clr.csv")%>%
  dplyr::rename(Featureid=V1)%>%
  dplyr::select(Featureid,all_of(md$sample_name))%>%
  gather(sample_name,clr,-Featureid)%>%
  left_join(.,md,by="sample_name")%>%
  left_join(.,taxon,by="Featureid")%>%
  mutate(memdx=factor(memdx, levels=c("nc","mci")))

md %>% count(memdx, visit_type)

#md_sub<-md%>%filter(sample_name %in% unique(counts_tab$sample_name))
#md%>%group_by(visit_type,memdx)%>%summarise(n=n())

prev235<-counts_tab%>%
  filter(Featureid=="G900539795")

p<-ggplot(prev235, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  scale_y_continuous(limits = c(-2, 14), breaks = seq(-2, 14, by = 2.5)) +
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Prevotella sp003447235 (G900539795)")+
  theme(legend.position = "none")

#ggsave("metagenomics/human_BEAM/SFR24_1206_prevotellasp003447235_clrabun_preb.pdf", plot=p,height=3.5, width=4)
ggsave("metagenomics/human_BEAM/SFR24_1206_prevotellasp003447235_clrabun_prea.pdf", plot=p,height=3.5, width=4)
pairwise.t.test(prev235$clr,prev235$memdx,p.adjust.method="fdr")
# nc    
# mci 0.0046
pairwise.wilcox.test(prev235$clr,prev235$memdx,p.adjust.method="fdr")
# nc    
# mci 0.0088


lreu<-counts_tab%>%
  filter(Featureid=="G000016825")

p<-ggplot(lreu, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Limosilactobacillus reuteri (G000016825)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Limosilactobacillusreuteri_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(lreu$clr,lreu$memdx,p.adjust.method="fdr")

# nc  
# mci 0.44

lmuco<-counts_tab%>%
  filter(Featureid=="G000248095")

p<-ggplot(lmuco, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Limosilactobacillus mucosae (G000248095)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Limosilactobacillusmucosae_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(lmuco$clr,lmuco$memdx,p.adjust.method="none")

# nc  
# mci 0.11

prevbaro<-counts_tab%>%
  filter(Featureid=="G000426585")

p<-ggplot(prevbaro, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Prevotella baroniae (G000426585)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Prevotellabaroniae_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(prevbaro$clr,prevbaro$memdx,p.adjust.method="none")

# nc   
# mci 0.098

prev365<-counts_tab%>%
  filter(Featureid=="G002251365")

p<-ggplot(prev365, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Prevotella sp002251365 (G002251365)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Prevotellasp002251365_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(prev365$clr,prev365$memdx,p.adjust.method="none")

# nc   
# mci 0.055

prev655<-counts_tab%>%
  filter(Featureid=="G004792655")

p<-ggplot(prev655, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Prevotella sp004792655 (G004792655)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Prevotellasp004792655_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(prev655$clr,prev655$memdx,p.adjust.method="none")
# nc  
# mci 0.78

prev295<-counts_tab%>%
  filter(Featureid=="G002251205")

p<-ggplot(prev295, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Prevotella sp002251295 (G002251205)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Prevotellasp002251295_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(prev295$clr,prev295$memdx,p.adjust.method="none")
# nc  
# mci 0.12

prev625<-counts_tab%>%
  filter(Featureid=="G902778655")

p<-ggplot(prev625, aes(x=memdx, y=clr, fill=memdx)) +
  geom_boxplot(alpha=0.3) + geom_dotplot(binaxis='y', stackdir='center',
                                         position=position_dodge(1), dotsize = 1.5) +
  scale_fill_manual(values=c("darkorange","purple4"))+
  theme_pubr()+
  labs(x="memdx",y="CLR abundance", title="Prevotella sp900318625 (G902778655)")+
  theme(legend.position = "none")

ggsave("metagenomics/human_BEAM/SFR24_1213_Prevotellasp900318625_clrabun.pdf", plot=p,height=3.5, width=4)
pairwise.wilcox.test(prev625$clr,prev625$memdx,p.adjust.method="none")
# nc  
# mci 0.71
