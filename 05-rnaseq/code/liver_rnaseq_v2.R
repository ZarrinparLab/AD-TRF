setwd("/mnt/zarrinpar/scratch/sfloresr/AD_EM/rnaseq/")

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(EnhancedVolcano)
library("ggdendro")
library(ggpubr)
library(ggbreak) 
library(cowplot)
############################################################

get_genelist<-function(df){
  annot<-fread("stringtie_annotations.tsv")%>%
    dplyr::select(stringtie_id,ncbi_id)
  df<-df%>%dplyr::rename(stringtie_id=V1)%>%left_join(.,annot,by="stringtie_id")
  
  original_gene_list <- df$log2FoldChange
  names(original_gene_list) <- df$ncbi_id
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

gse_nash_plot<-function(df){
  set.seed(1234)
  gene_list<-get_genelist(df)
  kk2 <- gseKEGG(geneList     = gene_list,
                 organism     = "mmu",
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")
  
  p<-dotplot(kk2, showCategory=15, split=".sign") + facet_grid(.~.sign)
  return(p)
}

gse_compareclust_allcond<-function(df1,df2,df3,df4,cond1,cond2,cond3,cond4){
  gene_list1<-get_genelist(df1)
  gene_list2<-get_genelist(df2)
  gene_list3<-get_genelist(df3)
  gene_list4<-get_genelist(df4)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"a-NTg-ALF","b-Tg-ALF"),
           condition=cond1)
  
  mydf2 <- data.frame(Entrez=names(gene_list2), FC=gene_list2)%>%
    mutate(direction=ifelse(FC<0,"a-NTg-ALF","b-Tg-ALF"),
           condition=cond2)
  
  mydf3 <- data.frame(Entrez=names(gene_list3), FC=gene_list3)%>%
    mutate(direction=ifelse(FC<0,"a-Tg-ALF","b-Tg-TRF"),
           condition=cond3)
  
  mydf4 <- data.frame(Entrez=names(gene_list4), FC=gene_list4)%>%
    mutate(direction=ifelse(FC<0,"a-Tg-ALF","b-Tg-TRF"),
           condition=cond4)
  
  mydf<-rbind(mydf1,mydf2,mydf3,mydf4)
  #mydf<-rbind(mydf1,mydf2)
  #formula_res <- compareCluster(Entrez~direction+condition, OrgDb = "org.Mm.eg.db", data=mydf,
  #                              fun="enrichGO",pvalueCutoff=0.05,pAdjustMethod="fdr")
  # formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
  #                               fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  formula_res <- compareCluster(Entrez~condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.1,pAdjustMethod="fdr")
  return(formula_res)
}

gse_compareclust_tgtrf<-function(df1,df2,cond1,cond2){
  gene_list1<-get_genelist(df1)
  gene_list2<-get_genelist(df2)
  
  mydf1 <- data.frame(Entrez=names(gene_list1), FC=gene_list1)%>%
    mutate(direction=ifelse(FC<0,"a-Tg-ALF","b-Tg-TRF"),
           condition=cond1)
  
  mydf2 <- data.frame(Entrez=names(gene_list2), FC=gene_list2)%>%
    mutate(direction=ifelse(FC<0,"a-Tg-ALF","b-Tg-TRF"),
           condition=cond2)

  mydf<-rbind(mydf1,mydf2)
  formula_res <- compareCluster(Entrez~direction+condition, organism = "mmu", data=mydf,
                                fun="enrichKEGG",pvalueCutoff=0.05,pAdjustMethod="fdr")
  return(formula_res)
}

#overall
deq_NTgTgALF<-fread('danielcarlin_rnaseq_files_noZT/liver_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF<-fread('danielcarlin_rnaseq_files_noZT/liver_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

gse_nash_plot(deq_NTgTgALF)+ggtitle("NTg-ALF vs Tg-ALF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_noZT/liver_TGvsNTGinALF-pval_05.kegg.pdf",height=6, width=10)

gse_nash_plot(deq_TgALFTgTRF)+ggtitle("Tg-ALF vs Tg-TRF")+
  scale_y_discrete(position = "right") 
ggsave("danielcarlin_rnaseq_files_noZT/liver_TRFvsALFinTG-pval_05.kegg.pdf",height=6, width=10)

#ZT0

deq_NTgTgALF<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_tgDE_in_ALF/ZT0_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

gse_nash_plot(deq_NTgTgALF)+ggtitle("ZT0 NTg-ALF vs Tg-ALF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/liver_zt_specific_tgDE_in_ALF/ZT0_TGvsNTGinALF-pval_05.kegg.pdf",height=6, width=10)

gse_nash_plot(deq_TgALFTgTRF)+ggtitle("ZT0 Tg-ALF vs Tg-TRF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG-pval_05.kegg.pdf",height=6, width=10)

#ZT12

deq_NTgTgALF<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_tgDE_in_ALF/ZT12_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

gse_nash_plot(deq_NTgTgALF)+ggtitle("ZT12 NTg-ALF vs Tg-ALF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/liver_zt_specific_tgDE_in_ALF/ZT12_TGvsNTGinALF-pval_05.kegg.pdf",height=6, width=12)

gse_nash_plot(deq_TgALFTgTRF)+ggtitle("ZT12 Tg-ALF vs Tg-TRF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG-pval_05.kegg.pdf",height=6, width=12)

deq_NTgTgALF_ZT0<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_tgDE_in_ALF/ZT0_TGvsNTGinALF.txt')%>%filter(padj<0.1)
deq_TgALFTgTRF_ZT0<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG.txt')%>%filter(padj<0.1)
deq_NTgTgALF_ZT12<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_tgDE_in_ALF/ZT12_TGvsNTGinALF.txt')%>%filter(padj<0.1)
deq_TgALFTgTRF_ZT12<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG.txt')%>%filter(padj<0.1)

allcond_AD<-gse_compareclust_allcond(deq_NTgTgALF_ZT0,deq_NTgTgALF_ZT12,deq_TgALFTgTRF_ZT0,deq_TgALFTgTRF_ZT12,
                                     "a-ZT0_compar1","b-ZT12_compar1","c-ZT0_compar2","d-ZT12_compar2")

allcond_AD_df<-allcond_AD@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
allcond_AD@compareClusterResult = allcond_AD@compareClusterResult[allcond_AD@compareClusterResult$ID %in% allcond_AD_df$ID, ]
#dotplot(allcond_AD,x="direction", showCategory=200) +
dotplot(allcond_AD,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondcompar_compclust.nohumandisease.top30.pdf",height=8, width=12)

allcond_AD<-gse_compareclust_allcond(deq_NTgTgALF_ZT0,deq_NTgTgALF_ZT12,deq_TgALFTgTRF_ZT0,deq_TgALFTgTRF_ZT12,
                                     "a-ZT0_geno","b-ZT12_geno","c-ZT0_treatment","d-ZT12_treatment")

allcond_AD_df<-allcond_AD@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))#%>%
  #filter(p.adjust<0.1)
allcond_AD@compareClusterResult = allcond_AD@compareClusterResult[allcond_AD@compareClusterResult$ID %in% allcond_AD_df$ID, ]
dotplot(allcond_AD,x="condition",size="p.adjust",color="Count",showCategory=15) +
#dotplot(allcond_AD,x="condition", showCategory=30) +
  facet_grid(category~., scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  scale_size_continuous(trans = "reverse")+
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(angle=90, vjust = 0.5))+
  labs(title="Liver")
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.nohumandisease.top30.pdf",height=9, width=8)
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.nohumandisease.top15.pdf",height=8, width=8)
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust_genecount.nohumandisease.top15.pdf",height=8, width=8)
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust_genecount.pcutoff0.1.nohumandisease.top15.pdf",height=8, width=8)

deq_TgALFTgTRF_ZT12<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF_ZT0<-fread('danielcarlin_rnaseq_files_byZT/liver_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

allcond_AD<-gse_compareclust_tgtrf(deq_TgALFTgTRF_ZT0,deq_TgALFTgTRF_ZT12,
                                     "ZT0","ZT12")

allcond_AD_df<-allcond_AD@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
allcond_AD@compareClusterResult = allcond_AD@compareClusterResult[allcond_AD@compareClusterResult$ID %in% allcond_AD_df$ID, ]
#dotplot(allcond_AD,x="direction",showCategory=15) +#size="p.adjust",color="Count", 
dotplot(allcond_AD,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  scale_size_continuous(trans = "reverse")+
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))+
  labs(title="Liver")
ggsave("danielcarlin_rnaseq_files_byZT/liver_enrichkegg_tgALFtgTRF_compar_compclust.nohumandisease.top30.pdf",height=8, width=9)
ggsave("danielcarlin_rnaseq_files_byZT/liver_enrichkegg_tgALFtgTRF_compar_compclust.nohumandisease.top15.pdf",height=8, width=9)
ggsave("danielcarlin_rnaseq_files_byZT/liver_enrichkegg_tgALFtgTRF_compar_compclust_genecount.nohumandisease.top15.pdf",height=8, width=9)

#look at bile secretion 

annot<-fread("stringtie_annotations.tsv")%>%
  dplyr::select(stringtie_id,ncbi_id)

#TRF_bileZT12<-allcond_AD_df%>%filter(Cluster=="d-ZT12_treatment" & grepl("Bile secretion",Description))%>%
#TRF_bileZT12<-allcond_AD_df%>%filter(Cluster=="a-ZT0_geno" & grepl("Fatty acid degradation",Description))%>%
#TRF_bileZT12<-allcond_AD_df%>%filter(Cluster=="a-ZT0_geno" & grepl("PPAR signaling pathway",Description))%>%
TRF_bileZT12<-allcond_AD_df%>%filter(Cluster=="c-ZT0_treatment" & grepl("Glutathione metabolism",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  #mutate(Kegg="Bile secretion")%>%
  #mutate(Kegg="Fatty acid degradation")%>%
  #mutate(Kegg="PPAR signaling pathway")%>%
  mutate(Kegg="Glutathione metabolism")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

#TRF_bileZT0<-allcond_AD_df%>%filter(Cluster=="c-ZT0_treatment" & grepl("Bile secretion",Description))%>%
#TRF_bileZT0<-allcond_AD_df%>%filter(Cluster=="a-ZT0_geno" & grepl("Fatty acid degradation",Description))%>%
#TRF_bileZT0<-allcond_AD_df%>%filter(Cluster=="a-ZT0_geno" & grepl("PPAR signaling pathway",Description))%>%
TRF_bileZT0<-allcond_AD_df%>%filter(Cluster=="c-ZT0_treatment" & grepl("Glutathione metabolism",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  #mutate(Kegg="Bile secretion")%>%
  #mutate(Kegg="Fatty acid degradation")%>%
  #mutate(Kegg="PPAR signaling pathway")%>%
  mutate(Kegg="Glutathione metabolism")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

md <- fread("/mnt/zarrinpar/scratch/sfloresr/AD_EM/rnaseq/per_sample.metadata.tsv")%>%
  filter(organ_code=="TI")%>%
  mutate(organ_code="LVR",
         sample_id = str_replace(sample_id, "_TI$", "_LVR"),
         geno_diet=paste(genotype,feeding_condition,sep="-"),
         fasted=ifelse(zt_time==12, "fasted","fed"))%>%
  filter(geno_diet!="NTg-TRF")

#setdiff(md$sample_id,genecounts$sample_id) #off by 2

#plot heatmap of these specific hits (Z-score)
genecounts<-fread("data_used_for_de/lvr.protein_coding.stringtie.gene_tpm_matrix.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% c(TRF_bileZT0$gene_id,TRF_bileZT12$gene_id))%>%
  gather(sample_id,TPM_counts,-gene_id) %>%
  inner_join(.,md,by="sample_id") %>%
  group_by(gene_id,geno_diet,zt_time)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id,zt_time)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>%
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)

# Run clustering
genecounts_df<-fread("data_used_for_de/lvr.protein_coding.stringtie.gene_tpm_matrix.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% c(TRF_bileZT0$gene_id,TRF_bileZT12$gene_id))%>%
  dplyr::select(gene_id,any_of(md$sample_id))

genecounts_matrix<-genecounts_df%>%column_to_rownames("gene_id") %>%as.matrix()
genecounts_dendro <- as.dendrogram(hclust(d = dist(x = genecounts_matrix)))
dendro_plot <- ggdendrogram(data = genecounts_dendro, rotate = TRUE)

genecounts_order <- order.dendrogram(genecounts_dendro)
genecounts$gene_id <- factor(x = genecounts$gene_id,
                             levels = genecounts_df$gene_id[genecounts_order], 
                             ordered = TRUE)

plt<-ggplot(genecounts,aes(x=geno_diet, y=gene_symbol)) +theme_classic()+
  geom_tile(aes(fill=Zscore))+
  scale_x_discrete(expand = c(0, 0))+scale_y_discrete(position = "right") + 
  facet_grid(~zt_time,scales="free",space="free",switch = "y")+
  theme(axis.ticks.y=element_blank(),panel.spacing.x=unit(1, "lines"),axis.text.y = element_text(size = 6),
        panel.spacing.y=unit(0.01, "lines"),strip.text.y = element_text(angle = 0))+
  #labs(title="Bile Secretion")+
  #labs(title="Biosynthesis of unsaturated fatty acids")+
  #labs(title="Fatty acid degradation")+
  #labs(title="PPAR signaling pathway")+
  labs(title="Glutathione metabolism")+
  scale_fill_viridis_c(option = "plasma")

ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsbilesecretion.heatmap_zscorecond.pdf",height=5, width=5)
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsFAdegradation.heatmap.pdf",height=3, width=5)
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsPPAR.heatmap.pdf",height=3, width=5)
ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsglutathionemetab.heatmap_zscorecond.pdf",height=3, width=5)
# ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.ZT012tgsbilesecretion.heatmap.pdf",height=5, width=5)
# ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsbilesecretion.heatmap.pdf",height=5, width=5)
# ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsFAdegradation.heatmap.pdf",height=3, width=5)
# ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsPPAR.heatmap.pdf",height=3, width=5)
# ggsave("danielcarlin_rnaseq_files_byZT/liver_allZT_enrichkegg_bycondnodirection_compclust.pcutoff0.1.ZT012tgsglutathionemetab.heatmap.pdf",height=3, width=5)

#plot indv genes

rna_indv_bxplt <- function(data,gene_id) {
  
  lab<-gene_id
  
  data%>%
    ggplot(aes(x=geno_diet, y=log10(TPM_counts+1), fill=geno_diet)) + geom_boxplot(alpha=0.3) +
    #geom_jitter(width=0.25, aes(colour = geno_diet))+
    geom_dotplot(binaxis='y', stackdir='center',
                 position=position_dodge(1)) +
    facet_wrap(~zt_time)+
    scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
    scale_colour_manual(values=c("darkorange","purple4","forestgreen"))+
    theme_pubr()+ theme(legend.position = "top")+
    labs(color="genotype",
         fill="genotype",
         y ="TPM counts",
         x ="genotype",
         title = stringr::str_wrap(lab, width=30))+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
  #ggsave(paste("danielcarlin_rnaseq_files_byZT/liver_bilesecretion/",gene_id,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("danielcarlin_rnaseq_files_byZT/FAdegradation/",gene_id,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("danielcarlin_rnaseq_files_byZT/PPAR/",gene_id,".pdf",sep=""), height=4, width=4)
  ggsave(paste("danielcarlin_rnaseq_files_byZT/GSH/",gene_id,".pdf",sep=""), height=4, width=4)
}

fasted_fed_pval<-function(df){
  df_fasted<-df%>%filter(fasted=="fasted")
  x<-pairwise.wilcox.test(df_fasted$TPM_counts, df_fasted$geno_diet,p.adjust.method="none")
  df_fasted_pvals<-data.frame(NtgALFvTgALF=x[["p.value"]][1],
                              NtgALFvTgTRF=x[["p.value"]][2],
                              TgALFvTgTRF=x[["p.value"]][4],
                              fasted="fasted")
  df_fed<-df%>%filter(fasted=="fed")
  y<-pairwise.wilcox.test(df_fed$TPM_counts, df_fed$geno_diet,p.adjust.method="none")
  df_fed_pvals<-data.frame(NtgALFvTgALF=y[["p.value"]][1],
                           NtgALFvTgTRF=y[["p.value"]][2],
                           TgALFvTgTRF=y[["p.value"]][4],
                           fasted="fed")
  pvals<-rbind(df_fasted_pvals,df_fed_pvals)
  return(pvals)
}

genecounts<-fread("data_used_for_de/lvr.protein_coding.stringtie.gene_tpm_matrix.csv")%>%
  dplyr::rename(gene_id=V1)%>%
  filter(gene_id %in% c(TRF_bileZT0$gene_id,TRF_bileZT12$gene_id))%>%
  gather(sample_id,TPM_counts,-gene_id) %>%
  inner_join(.,md,by="sample_id")

rna_nested <- genecounts %>% 
  group_by(gene_id) %>% 
  nest()

rna_plots <- 
  rna_nested %>% 
  mutate(plot = map2(data, gene_id,  ~ rna_indv_bxplt(.x,.y)))

rna_pvals <- 
  rna_nested %>% 
  mutate(pvals = map(data, fasted_fed_pval))%>%
  dplyr::select(-data)%>%
  unnest()
write.table(rna_pvals,"danielcarlin_rnaseq_files_byZT/liver_bilesecretion/liver_bilesecretion-condfastedfed-wilcox-pvals.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(rna_pvals,"danielcarlin_rnaseq_files_byZT/FAdegradation/liver_FAdegradation-condfastedfed-wilcox-pvals.txt",sep = "\t",row.names = FALSE,quote=FALSE)
write.table(rna_pvals,"danielcarlin_rnaseq_files_byZT/PPAR/liver_PPAR-condfastedfed-wilcox-pvals.txt",sep = "\t",row.names = FALSE,quote=FALSE)


#clean
acna1<-genecounts%>%
  filter(gene_id=="gene9058|Acnat1")

acnat1_plt<-acna1%>%ggplot(aes(x=geno_diet, y=log10(TPM_counts+1), fill=geno_diet)) + geom_boxplot(alpha=0.3) +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1)) +
  facet_wrap(~zt_time)+
  scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
  scale_colour_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr()+ theme(legend.position = "top")+
  #scale_y_continuous(limits=c(11,14.5))+
  labs(color="genotype",
       fill="genotype",
       y ="log10(TPM counts)",
       x ="genotype",
       title = "Acnat1")+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
ggsave("danielcarlin_rnaseq_files_byZT/liver_bilesecretion/acnat1_cln.pdf", height=4.5, width=5)

# baat_fasted<-genecounts%>%
#   filter(gene_id=="gene9061|Baat" & fasted=="fasted")
# baat_fed<-genecounts%>%
#   filter(gene_id=="gene9061|Baat" & fasted=="fed")
# 
# plotA<-baat_fasted%>%ggplot(aes(x=geno_diet, y=TPM_counts, fill=geno_diet)) + geom_boxplot(alpha=0.3) +
#   geom_dotplot(binaxis='y', stackdir='center',
#                position=position_dodge(1)) +
#   facet_wrap(~fasted,scales="free_y")+
#   scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
#   scale_colour_manual(values=c("darkorange","purple4","forestgreen"))+
#   scale_y_continuous(limits=c(11.5,13))+
#   theme_pubr()+ theme(legend.position = "top")+
#   labs(color="genotype",
#        fill="genotype",
#        y ="TPM counts",
#        x ="genotype",
#        title = "baat")+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
# ggsave("liver_bilesecretion/baatfasted_cln.pdf", height=4.5, width=3)

# plotB<-baat_fed%>%ggplot(aes(x=geno_diet, y=TPM_counts, fill=geno_diet)) + geom_boxplot(alpha=0.3) +
#   geom_dotplot(binaxis='y', stackdir='center',
#                position=position_dodge(1),dotsize = 0.20) +
#   facet_wrap(~fasted,scales="free_y")+
#   scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
#   scale_colour_manual(values=c("darkorange","purple4","forestgreen"))+
#   scale_y_continuous(limits=c(4,13),breaks=c(4,4.5,5,11.75,12,12.25, 12.5,12.75,13))+
#   scale_y_break(c(5, 11.5),scales=7, expand=FALSE) +
#   theme_pubr()+ theme(legend.position = "top")+
#   labs(color="genotype",
#        fill="genotype",
#        y ="TPM counts",
#        x ="genotype",
#        title="baat")+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
# ggsave("liver_bilesecretion/baatfed_cln.pdf", height=4.5, width=3.5)
