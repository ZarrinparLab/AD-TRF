setwd("/mnt/zarrinpar/scratch/sfloresr/AD_EM/rnaseq/")

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(EnhancedVolcano)
library("ggdendro")
library(ggpubr)
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
                                fun="enrichKEGG",pvalueCutoff=1,pAdjustMethod="fdr")
  return(formula_res)
}

#overall

deq_NTgTgALF<-fread('danielcarlin_rnaseq_files_noZT/ti_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF<-fread('danielcarlin_rnaseq_files_noZT/ti_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

gse_nash_plot(deq_NTgTgALF)+ggtitle("NTg-ALF vs Tg-ALF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_noZT/TI_TGvsNTGinALF-pval_05.kegg.pdf",height=6, width=10)

gse_nash_plot(deq_TgALFTgTRF)+ggtitle("Tg-ALF vs Tg-TRF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_noZT/TI_TRFvsALFinTG-pval_05.kegg.pdf",height=6, width=12)

#ZT0

deq_NTgTgALF<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_tgDE_in_ALF/ZT0_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

gse_nash_plot(deq_NTgTgALF)+ggtitle("ZT0 NTg-ALF vs Tg-ALF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/ti_zt_specific_tgDE_in_ALF/ZT0_TGvsNTGinALF-pval_05.kegg.pdf",height=6, width=10)

gse_nash_plot(deq_TgALFTgTRF)+ggtitle("ZT0 Tg-ALF vs Tg-TRF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG-pval_05.kegg.pdf",height=6, width=10)

#ZT12

deq_NTgTgALF<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_tgDE_in_ALF/ZT12_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

#not significant
# gse_nash_plot(deq_NTgTgALF)+ggtitle("ZT12 NTg-ALF vs Tg-ALF")+
#   scale_y_discrete(position = "right")

gse_nash_plot(deq_TgALFTgTRF)+ggtitle("ZT12 Tg-ALF vs Tg-TRF")+
  scale_y_discrete(position = "right")
ggsave("danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG-pval_05.kegg.pdf",height=6, width=10)

deq_NTgTgALF_ZT0<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_tgDE_in_ALF/ZT0_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF_ZT0<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)
deq_NTgTgALF_ZT12<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_tgDE_in_ALF/ZT12_TGvsNTGinALF.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF_ZT12<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

allcond_AD<-gse_compareclust_allcond(deq_NTgTgALF_ZT0,deq_NTgTgALF_ZT12,deq_TgALFTgTRF_ZT0,deq_TgALFTgTRF_ZT12,
                                     "a-ZT0_compar1","b-ZT12_compar1","c-ZT0_compar2","d-ZT12_compar2")

allcond_AD_df<-allcond_AD@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
allcond_AD@compareClusterResult = allcond_AD@compareClusterResult[allcond_AD@compareClusterResult$ID %in% allcond_AD_df$ID, ]
dotplot(allcond_AD,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("danielcarlin_rnaseq_files_byZT/ti_allZT_enrichkegg_bycondcompar_compclust.nohumandisease.top30.pdf",height=8, width=12)

allcond_AD<-gse_compareclust_allcond(deq_NTgTgALF_ZT0,deq_NTgTgALF_ZT12,deq_TgALFTgTRF_ZT0,deq_TgALFTgTRF_ZT12,
                                     "a-ZT0_geno","b-ZT12_geno","c-ZT0_treatment","d-ZT12_treatment")

allcond_AD_df<-allcond_AD@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
allcond_AD@compareClusterResult = allcond_AD@compareClusterResult[allcond_AD@compareClusterResult$ID %in% allcond_AD_df$ID, ]
dotplot(allcond_AD,x="condition",size="p.adjust",color="Count",showCategory=15) +
#dotplot(allcond_AD,x="condition", showCategory=30) +
  facet_grid(category~., scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  scale_size_continuous(trans = "reverse")+
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(angle=90, vjust = 0.5))+
  labs(title="Ileum")
ggsave("danielcarlin_rnaseq_files_byZT/ti_allZT_enrichkegg_bycondnodirection_compclust.nohumandisease.top30.pdf",height=9, width=8)
ggsave("danielcarlin_rnaseq_files_byZT/ti_allZT_enrichkegg_bycondnodirection_compclust.nohumandisease.top15.pdf",height=8, width=8)
ggsave("danielcarlin_rnaseq_files_byZT/ti_allZT_enrichkegg_bycondnodirection_compclust_genecount.nohumandisease.top15.pdf",height=8, width=8)
ggsave("danielcarlin_rnaseq_files_byZT/ti_allZT_enrichkegg_bycondnodirection_compclust_genecount.pcutoff0.1.nohumandisease.top15.pdf",height=8, width=8)


deq_TgALFTgTRF_ZT12<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT12_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)
deq_TgALFTgTRF_ZT0<-fread('danielcarlin_rnaseq_files_byZT/ti_zt_specific_treatDE_in_TG/ZT0_TRFvsALFinTG.txt')%>%filter(pvalue<0.05)

allcond_AD<-gse_compareclust_tgtrf(deq_TgALFTgTRF_ZT0,deq_TgALFTgTRF_ZT12,
                                   "ZT0","ZT12")

allcond_AD_df<-allcond_AD@compareClusterResult%>%
  filter(!(category=="Human Diseases" | is.na(category)))
allcond_AD@compareClusterResult = allcond_AD@compareClusterResult[allcond_AD@compareClusterResult$ID %in% allcond_AD_df$ID, ]
dotplot(allcond_AD,x="direction",size="p.adjust",color="Count",showCategory=15) +# 
#dotplot(allcond_AD,x="direction", showCategory=30) +
  facet_grid(category~condition, scales="free",space='free_y') + 
  scale_y_discrete(position = "right") +
  scale_size_continuous(trans = "reverse")+
  theme(axis.text.y = element_text(size=8),
        legend.position = "left",
        axis.text.x=element_text(vjust = 0.5))
ggsave("danielcarlin_rnaseq_files_byZT/ti_enrichkegg_tgALFtgTRF_compar_compclust.nohumandisease.top30.pdf",height=8, width=9)
ggsave("danielcarlin_rnaseq_files_byZT/ti_enrichkegg_tgALFtgTRF_compar_compclust.nohumandisease.top15.pdf",height=8, width=9)
ggsave("danielcarlin_rnaseq_files_byZT/ti_enrichkegg_tgALFtgTRF_compar_compclust_genecount.nohumandisease.top15.pdf",height=8, width=9)

#look at bile secretion 

annot<-fread("stringtie_annotations.tsv")%>%
  dplyr::select(stringtie_id,ncbi_id)

TRF_bileZT12<-allcond_AD_df%>%filter(Cluster=="d-ZT12_treatment" & grepl("Biosynthesis of unsaturated fatty acids",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="Biosynthesis of unsaturated fatty acids")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

TRF_bileZT0<-allcond_AD_df%>%filter(Cluster=="c-ZT0_treatment" & grepl("Biosynthesis of unsaturated fatty acids",Description))%>%
  separate_rows(geneID)%>% group_by(geneID)%>%summarise(Description = paste(Description, collapse="|"))%>%
  dplyr::rename(ncbi_id=geneID)%>%mutate(ncbi_id=as.numeric(ncbi_id))%>%
  left_join(.,annot,by="ncbi_id")%>%
  mutate(Kegg="Bile secretion")%>%
  dplyr::rename(gene_id=stringtie_id)%>%distinct(gene_id, .keep_all = TRUE)

md <- fread("/mnt/zarrinpar/scratch/sfloresr/AD_EM/rnaseq/per_sample.metadata.tsv")%>%
  filter(organ_code=="TI")%>%
  mutate(geno_diet=paste(genotype,feeding_condition,sep="-"),
         fasted=ifelse(zt_time==12, "fasted","fed"))%>%
  filter(geno_diet!="NTg-TRF")

#setdiff(md$sample_id,genecounts$sample_id) #off by 2

#plot heatmap of these specific hits (Z-score)
genecounts<-fread("data_used_for_de/ti.ad_trf_ti.stringtie.gene.tpm.csv")%>%
  filter(gene_id %in% c(TRF_bileZT0$gene_id,TRF_bileZT12$gene_id))%>%
  gather(sample_id,TPM_counts,-gene_id) %>%
  mutate(sample_id = str_replace_all(sample_id, c("-" = "_", "1$" = "")))%>%
  inner_join(.,md,by="sample_id") %>%
  group_by(gene_id,geno_diet,zt_time)%>%dplyr::summarise(mn_TPM=mean(TPM_counts))%>%
  group_by(gene_id,zt_time)%>%mutate(Zscore=(mn_TPM - mean(mn_TPM))/sd(mn_TPM))%>% #cant do zscore at zt_time level
  separate(gene_id, c("gene_number","gene_symbol"), remove=FALSE)

# Run clustering
genecounts_df<-fread("data_used_for_de/ti.ad_trf_ti.stringtie.gene.tpm.csv")%>%
  filter(gene_id %in% c(TRF_bileZT0$gene_id,TRF_bileZT12$gene_id))%>%
  rename_with(~ str_replace_all(.x, c("-" = "_", "1$" = "")))%>%
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
  labs(title="Biosynthesis of unsaturated fatty acids")+
  scale_fill_viridis_c(option = "plasma")

ggsave("danielcarlin_rnaseq_files_byZT/TI_allZT_enrichkegg_bycondnodirection_compclust.ZT012tgsbilesecretion.heatmap.pdf",height=3, width=5)
ggsave("danielcarlin_rnaseq_files_byZT/TI_allZT_enrichkegg_bycondnodirection_compclust.pvalcutoff0.1.ZT012tgsbilesecretion.heatmap.pdf",height=3, width=5)


