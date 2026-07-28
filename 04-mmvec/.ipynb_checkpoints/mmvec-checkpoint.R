setwd("/mnt/zarrinpar/scratch/sfloresr/AD_EM/")

library(tidyverse)
library(data.table)
library("qiime2R")
library("Biostrings")
library(ggcorrplot)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(scatterplot3d)
library(VennDiagram)
library(gplots)
library(qiime2R)
library(ggpubr)
library(ggbiplot)
library(ggsci)
######################################################################

#load mtb
mtb_md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces")%>%
  mutate(ATTRIBUTE_sex=ifelse(ATTRIBUTE_sex=="F","female",ATTRIBUTE_sex))%>%
  dplyr::rename(sample_name=SubjectIdentifierAsRecorded)%>%
  mutate(zt_time=as.numeric(sub("ZT","",ATTRIBUTE_ZT_collection_time)))%>%
  arrange(zt_time)%>%
  arrange(ATTRIBUTE_host_subject_id)%>%
  dplyr::select(sample_name,everything())
write.table(mtb_md, "mmvec/mmvec_mtb_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_collapseiin_sum_corrgroup.csv")%>%
  dplyr::rename(FeatureID=row.ID)%>%
  dplyr::select(FeatureID,all_of(mtb_md$sampleid))

names(mtb) <-c("FeatureID",mtb_md$sample_name)
write.table(mtb, "mmvec/mmvec_mtb_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#load mtg
mtg_md<-fread("jointRPCA/stool.metagenomic.metadata.with_amyloid.tsv")%>%
  dplyr::rename(sampleid=sample_name)%>%
  mutate(sample_name = sub("^[^.]+\\.", "", sampleid))%>%
  #filter(sample_name %in% mtb_md$sample_name)%>%
  arrange(zt_time)%>%
  arrange(host_subject_id)%>%
  dplyr::select(sample_name,everything())
write.table(mtg_md, "mmvec/mmvec_mtg_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv")%>%
  dplyr::rename(FeatureID='OTU ID')%>%
  dplyr::select(FeatureID,all_of(mtg_md$sampleid))%>%
  rename_with(~str_replace_all(., "14708.", ""))

write.table(mtg, "mmvec/mmvec_mtg_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#shared metadata with training-test column
shared<-intersect(mtg_md$sample_name,mtb_md$sample_name)

md<-mtg_md%>%
  filter(sample_name %in% shared)%>%
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"))%>%
  dplyr::select(sample_name, host_subject_id,genotype_diet,zt_collection_time,fasted)

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sample_name %in% testsamps$sample_name, "Test","Train"))
write.table(md, "mmvec/mmvec_comb_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just keep columns that are shared
mtg_shared <- mtg%>%
  dplyr::select(FeatureID,any_of(shared))
write.table(mtg_shared, "mmvec/mmvec_mtg_shared_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtb_shared <- mtb%>%
  dplyr::select(FeatureID,any_of(shared))
write.table(mtb_shared, "mmvec/mmvec_mtb_shared_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#plot by fasted and fed
md_fasted<-md%>%filter(fasted=="fasted")
md_fed<-md%>%filter(fasted=="fed")

mtg_shared_fasted <- mtg%>%
  dplyr::select(FeatureID,any_of(md_fasted$sample_name))
write.table(mtg_shared_fasted, "mmvec/mmvec_mtg_shared_fasted_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtg_shared_fed <- mtg%>%
  dplyr::select(FeatureID,any_of(md_fed$sample_name))
write.table(mtg_shared_fed, "mmvec/mmvec_mtg_shared_fed_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb_shared_fasted <- mtb%>%
  dplyr::select(FeatureID,any_of(md_fasted$sample_name))
write.table(mtb_shared_fasted, "mmvec/mmvec_mtb_shared_fasted_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtb_shared_fed <- mtb%>%
  dplyr::select(FeatureID,any_of(md_fed$sample_name))
write.table(mtb_shared_fed, "mmvec/mmvec_mtb_shared_fed_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#taxon and mtb metadata

mtb_id <- fread("metabolomics//data_files/AD_GNPS_BA/nf_output/library/library_matches_iin_corrgroup.csv") %>%
  mutate(across(everything(), ~ str_extract(.x, "(?<!NA)\\b[^|]+\\b")))%>%
  mutate(across(everything(), ~ ifelse(.x == "NA" | .x == "N/A", NA, .x)))%>%
  mutate(Compound_Name = str_replace(Compound_Name, " \\(delta.*", ""),
         subclass = ifelse(str_detect(Compound_Name,"bile acid"),"Bile acids, alcohols and derivatives",subclass))
  
mtb_id_cleaned <- as.data.frame(lapply(mtb_id, function(x) gsub('"', '', x)))
write.table(mtb_id_cleaned, "mmvec/mmvec_mtb_feat_annotation.txt",sep = "\t",row.names = FALSE, quote=FALSE)

taxon<-mtg%>%
  dplyr::select(FeatureID)%>%
  separate(FeatureID,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="; ", remove=FALSE)
taxon_cleaned <- as.data.frame(lapply(taxon, function(x) gsub('"', '', x)))
write.table(taxon, "mmvec/mmvec_mtg_feat_annotation.txt",sep = "\t",row.names = FALSE, quote=FALSE)


#create biplot
micro_int<-c("s__Turicibacter sanguinis","s__Turicibacter sp001543345",
             "s__Limosilactobacillus antri","s__Limosilactobacillus sp900156885",
             "s__Prevotella sp900119285","s__Prevotella ruminicola","s__Prevotella sp000688375","s__Prevotella nigrescens","s__Prevotella loescheii","s__Prevotella sp002372235","s__Prevotella sp003447235",
             "s__Dorea formicigenerans","Dorea_A longicatena",
             "s__Lawsonibacter asaccharolyticus","s__Lawsonibacter sp000492175","s__Lawsonibacter sp000177015","s__Lawsonibacter sp014287875",
             "s__Oscillibacter valericigenes","s__Oscillibacter ruminantium","s__Prevotella sp900314915","s__Faecalibacterium prausnitzii_C",
             "s__Roseburia intestinalis","s__Blautia_A avistercoris","s__Eubacterium_Q ruminantium_A",
             "s__Ligilactobacillus equi",
             "s__Limosilactobacillus timonensis","s__Limosilactobacillus sp900156885","s__Limosilactobacillus vaginalis","s__Limosilactobacillus pontis","s__Limosilactobacillus panis",
             "s__Limosilactobacillus mucosae","s__Limosilactobacillus ingluviei","s__Limosilactobacillus antri","s__Limosilactobacillus secaliphilus","s__Limosilactobacillus reuteri",
             "s__Limosilactobacillus coleohominis_A","s__Limosilactobacillus fermentum","s__Limosilactobacillus equigenerosi",
             "s__Prevotella sp000436915","s__Prevotella sp900317685","s__Prevotella melaninogenica","s__Prevotella disiens","s__Prevotella corporis","s__Prevotella sp900318625",
             "s__Prevotella sp004555245","s__Prevotella oris","s__Prevotella sp902801425","s__Prevotella sp900199655","s__Prevotella sp900110895","s__Prevotella sp900110085",
             "s__Prevotella sp900100635","s__Prevotella sp000433175","s__Prevotella ruminicola_B","s__Prevotella maculosa","s__Prevotella brevis","s__Prevotella sp900313615",
             "s__Prevotella scopos","s__Prevotella bergensis","s__Prevotella sp000424185","s__Prevotella multisaccharivorax")
micro_int2<-c("s__Dorea formicigenerans","s__Lawsonibacter sp000492175","s__Lawsonibacter sp014287875",
             "s__Limosilactobacillus reuteri","s__Limosilactobacillus fermentum",
             "s__Prevotella disiens","s__Prevotella corporis")
micro_int3<-c("s__Lawsonibacter sp000492175","s__Limosilactobacillus reuteri","s__Limosilactobacillus fermentum")
# micro_int<-c("s__Limosilactobacillus mucosae","s__Limosilactobacillus reuteri","s__Limosilactobacillus secaliphilus","s__Prevotella baroniae",
#              "s__Prevotella sp002251295","s__Prevotella sp003447235","s__Prevotella sp004792655","s__Prevotella sp900318625","s__Prevotella sp002251365")


ord <- read_qza("mmvec/AD_summary_shared/conditional_biplot.qza")
#ord <- read_qza("mmvec/AD_summary_shared_fasted/conditional_biplot.qza")
#ord <- read_qza("mmvec/AD_summary_shared_fed/conditional_biplot.qza")
mtb_id<-fread("mmvec/mmvec_mtb_feat_annotation.txt")%>%
  mutate(subclass_abr=ifelse(subclass %in% c("Bile acids, alcohols and derivatives","Amino acids, peptides, and analogues",
                                               "Fatty acids and conjugates",NA),subclass,"Other"))

#mtb_id%>%group_by(subclass)%>%summarise(n=n())%>%arrange(desc(n))

rpca_mmvec<-ord$data$Vectors %>%
  select(SampleID, PC1, PC2)%>%
  dplyr::rename(FeatureID=SampleID)%>%
  left_join(.,mtb_id,by="FeatureID")

vectors <-ord$data$Species%>%
  mutate(subclass_abr=NA,
         FeatureID=sub(".*(g__.*)", "\\1", FeatureID))%>%
  arrange(desc(abs(PC1)))%>%
  filter(str_detect(FeatureID, paste(micro_int, collapse = "|"))) %>% 
  #filter(str_detect(FeatureID, paste(micro_int2, collapse = "|"))) %>% 
  #filter(str_detect(FeatureID, paste(micro_int3, collapse = "|"))) %>%
  #slice_head(n = 10) %>%
  #mutate(PC1=PC1/3300,PC2=PC2/3300,PC3=PC3/3300)
  mutate(PC1=PC1/3700,PC2=PC2/3700,PC3=PC3/3700)      

p <- rpca_mmvec %>%
  ggplot(aes(x = PC1, y = PC2, fill = subclass_abr)) +
  geom_point(data = subset(rpca_mmvec, is.na(subclass_abr)), aes(x = PC1, y = PC2, fill = subclass_abr),
             alpha = 1.0, shape = 21, size = 2) +
  geom_point(data = subset(rpca_mmvec, !is.na(subclass_abr)), aes(x = PC1, y = PC2, fill = subclass_abr),
             alpha = 1.0, shape = 21, size = 2) +
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(type = "closed", length = unit(0.05, "inches")),
               color = "black", size = 0.5) +
  geom_label(data = vectors, aes(x = PC1, y = PC2, label = FeatureID),
             size = 3, color = "black", fill = "white", label.size = 0.2, 
             label.padding = unit(0.15, "lines"), vjust = -0.5) + 
  scale_fill_manual(values = c("#00AEEF", "#EC008C", "#FFF200", "darkgray", "lightgray"),
                    na.value = "lightgray") +  
  #xlim(-0.1,0.05)+
  #xlim(-0.12,0.05)+
  theme_pubr() +
  theme(legend.position = "right")+
  labs(x = paste("PC1 (", round(ord$data$ProportionExplained$PC1 * 100, digits = 2), "%)", sep = ""),
       y = paste("PC2 (", round(ord$data$ProportionExplained$PC2 * 100, digits = 2), "%)", sep = ""),
       color = "conditions")+
       #title="Fasted") +
       #title="Fed") +
  theme(plot.title = element_text(face = "bold"))

ggsave("mmvec/AD_summary_shared_fed/SFR25_0402_mmvecbiplot_top60micro_top3subclassmtb_fed.pdf", plot=p,height=4, width=8.5)
ggsave("mmvec/AD_summary_shared_fasted/SFR25_0402_mmvecbiplot_top60micro_top3subclassmtb_fasted.pdf", plot=p,height=4, width=8.5)
ggsave("mmvec/AD_summary_shared/SFR25_0402_mmvecbiplot_60microint_top3subclassmtb.pdf", plot=p,height=4, width=8.5)
ggsave("mmvec/AD_summary_shared/SFR25_0402_mmvecbiplot_7microint_top3subclassmtb.pdf", plot=p,height=4, width=8.5)
ggsave("mmvec/AD_summary_shared/SFR25_0402_mmvecbiplot_3microint_top3subclassmtb.pdf", plot=p,height=4, width=8.5)

# ggsave("mmvec/AD_summary_shared_fed/SFR24_1018_mmvecbiplot_3micro_top3subclassmtb_fed.pdf", plot=p,height=4, width=8.5)
# ggsave("mmvec/AD_summary_shared_fasted/SFR24_1018_mmvecbiplot_3micro_top3subclassmtb_fasted.pdf", plot=p,height=4, width=8.5)

#create conditional heatmap

#Microbes of interest

#cond<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv")%>%
#cond<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv")%>%
cond<- fread("mmvec/AD_summary_shared/conditionals/conditionals.tsv")%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  filter(Species %in% micro_int)%>%
  #filter(Species %in% micro_int2)%>%
  #filter(Species %in% micro_int3)%>%
  dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
  column_to_rownames("Species")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(!is.na(Compound_Name))%>%
  mutate(label=paste(FeatureID," ",signif(mn_mz,5)," ",Compound_Name))%>%
  dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-subclass_abr,-FeatureID)
write.table(cond, "mmvec/AD_summary_shared/conditionals/conditionals_microint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

cond_mat<-cond%>%column_to_rownames("label")%>%as.matrix()
p<-pheatmap(cond_mat,main="Microbes of Int all Mtbs")
#ggsave("mmvec/AD_summary_shared_fed/mmvec_heatmap_60microb_allmtb.pdf", plot=p,height=25, width=35)
#ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_60microb_allmtb.pdf", plot=p,height=25, width=35)
#ggsave("mmvec/AD_summary_shared/mmvec_heatmap_60microb_allmtb.pdf", plot=p,height=25, width=35)
#ggsave("mmvec/AD_summary_shared/mmvec_heatmap_7microb_allmtb.pdf", plot=p,height=25, width=28)
ggsave("mmvec/AD_summary_shared/mmvec_heatmap_3newmicrob_allmtb.pdf", plot=p,height=25, width=26)

#add just mtb of interest
bdm_fasted<-fread("metabolomics/birdman/collapseiin-log/fecal-bdm-justcredannot-tgTRFtgALF-fold0.3-fasted-log.txt")
bdm_fed<-fread("metabolomics/birdman/collapseiin-log/fecal-bdm-justcredannot-tgTRFtgALF-fold0.3-fed-log.txt")

#cond<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv")%>%
#cond<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv")%>%
cond<- fread("mmvec/AD_summary_shared/conditionals/conditionals.tsv")%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  #filter(Species %in% micro_int)%>%
  #filter(Species %in% micro_int2)%>%
  #filter(Species %in% micro_int3)%>%
  dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
  column_to_rownames("Species")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  #filter(FeatureID %in% bdm_fasted$FeatureID)%>%
  filter(FeatureID %in% c(bdm_fed$FeatureID,bdm_fasted$FeatureID))%>%
  mutate(label=paste(FeatureID," ",signif(mn_mz,5)," ",Compound_Name))%>%
  dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-subclass_abr,-FeatureID)

cond_mat<-cond%>%column_to_rownames("label")%>%as.matrix()
p<-pheatmap(cond_mat,main="microbe int vs. fed/fasted mtb")
#ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_60microb_fastedfedmtb.pdf", plot=p,height=10, width=20)
#ggsave("mmvec/AD_summary_shared_fed/mmvec_heatmap_60microb_fastedfedmtb.pdf", plot=p,height=10, width=20)
#ggsave("mmvec/AD_summary_shared/mmvec_heatmap_60microb_fastedfedmtb.pdf", plot=p,height=10, width=20)
ggsave("mmvec/AD_summary_shared/mmvec_heatmap_7microb_fastedfedmtb.pdf", plot=p,height=10, width=15)
ggsave("mmvec/AD_summary_shared/mmvec_heatmap_3newmicrob_fastedfedmtb.pdf", plot=p,height=10, width=14)

# p<-pheatmap(cond_mat,main="microbe int vs. fasted mtb")
# ggsave("mmvec/AD_summary_shared/mmvec_heatmap_9microb_fastedmtb.pdf", plot=p,height=8, width=10)
# 
# p<-pheatmap(cond_mat,main="microbe int vs. fed mtb")
# ggsave("mmvec/AD_summary_shared/mmvec_heatmap_9microb_fedmtb.pdf", plot=p,height=8, width=15)


#just show the fasted and fed results for 7 microbes of interest (subset 3)

cond_fasted<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv")%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  #filter(Species %in% micro_int2)%>%
  #filter(Species %in% micro_int3)%>%
  #dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
  dplyr::select(FeatureID,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-Species)%>%
  #column_to_rownames("Species")%>%
  column_to_rownames("FeatureID")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(FeatureID %in% c(bdm_fed$FeatureID,bdm_fasted$FeatureID))%>%
  mutate(label=paste(FeatureID," ",mn_mz," ",Compound_Name))%>%
  #mutate(label=paste(FeatureID," ",signif(mn_mz,5)," ",Compound_Name))%>%
  dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-FeatureID)
  # dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
  #               -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-subclass_abr,-FeatureID)%>%
  #filter(label!="13644   387.2   3,5,5-Trimethyl-4beta-hydroxy-4-[3-(beta-D-glucopyranosyloxy)-1-butenyl]-2-cyclohexene-1-one")
write.table(cond_fasted, "mmvec/AD_summary_shared_fasted/conditionals/conditionals_fasted_all.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#cond_fasted[cond_fasted$label == "13644   387.2   3,5,5-Trimethyl-4beta-hydroxy-4-[3-(beta-D-glucopyranosyloxy)-1-butenyl]-2-cyclohexene-1-one", "s__Limosilactobacillus reuteri"] <- -4

cond_fed<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv")%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  #filter(Species %in% micro_int2)%>%
  #filter(Species %in% micro_int3)%>%
  #dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
  dplyr::select(FeatureID,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-Species)%>%
  #column_to_rownames("Species")%>%
  column_to_rownames("FeatureID")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(FeatureID %in% c(bdm_fed$FeatureID,bdm_fasted$FeatureID))%>%
  mutate(label=paste(FeatureID," ",mn_mz," ",Compound_Name))%>%
  #mutate(label=paste(FeatureID," ",signif(mn_mz,5)," ",Compound_Name))%>%
  dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-FeatureID)
# dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
#               -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-subclass_abr,-FeatureID)%>%
#filter(label!="13644   387.2   3,5,5-Trimethyl-4beta-hydroxy-4-[3-(beta-D-glucopyranosyloxy)-1-butenyl]-2-cyclohexene-1-one")
write.table(cond_fed, "mmvec/AD_summary_shared_fed/conditionals/conditionals_fed_all.tsv",sep = "\t",row.names = FALSE, quote=FALSE)


comb_fasted<-left_join(cond_fasted,cond_fed,by="label")#%>%column_to_rownames("label")%>%as.matrix()
write.table(comb_fasted, "mmvec/AD_summary_shared_fasted/conditionals/conditionals_fastedfed_microint3.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

p<-pheatmap(comb_fasted,main="microbe int vs. fastedfed mtb", cluster_cols=FALSE)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_7microb_fastedfedcombmtb.pdf", plot=p,height=10, width=17)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_7microb_fastedfedcombmtb_cappingoutlier.pdf", plot=p,height=10, width=17)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_7microb_fastedfedcombmtb_rmoutlier.pdf", plot=p,height=10, width=17)

ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_3newmicrob_fastedfedcombmtb.pdf", plot=p,height=10, width=15)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_3newmicrob_fastedfedcombmtb_cappingoutlier.pdf", plot=p,height=10, width=15)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_3newmicrob_fastedfedcombmtb_rmoutlier.pdf", plot=p,height=10, width=15)
######################################################################
mtg_shared<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(FeatureID=`OTU ID`)
mtb_shared<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_collapseiin_sum_corrgroup_log.csv")%>%
  dplyr::rename(FeatureID=row.ID)%>%
  dplyr::select(FeatureID,all_of(mtb_md$sampleid))
names(mtb_shared) <-c("FeatureID",mtb_md$sample_name)
md<-fread("mmvec/mmvec_comb_metadata.txt")

lreu <- mtg_shared %>%filter(str_detect(FeatureID, "Limosilactobacillus reuteri"))%>%
  gather(sample_name,counts_x, -FeatureID)

lawsp2175 <- mtg_shared %>%filter(str_detect(FeatureID, "Lawsonibacter sp000492175"))%>%
  gather(sample_name,counts_x, -FeatureID)

spearmancorr_mmvec<-function(mtg,mtb,mtb_id){
  
  mtb <- mtb %>% filter(FeatureID==mtb_id)%>%
    gather(sample_name,counts_y, -FeatureID)
  shared<-intersect(mtg$sample_name,mtb$sample_name)
  
  mtg_mtb<-left_join(mtg,mtb, by="sample_name")%>%
    left_join(.,md, by="sample_name")%>%
    filter(sample_name %in% shared)
  p<-ggplot(mtg_mtb, aes(x=counts_x, y=counts_y, color=fasted, fill=fasted)) + geom_point()+ 
    stat_cor(method="spearman")+
    geom_smooth(method = "lm",se = TRUE,na.rm = TRUE)+
    scale_fill_manual(values=c("#BE1E2D","gray30"))+
    scale_color_manual(values=c("#BE1E2D","gray30"))+
    theme_pubr()+
    scale_y_continuous(limits=c(-1,17))+
    #labs(x="Lreu CLR",y=mtb_id)
    labs(x="Law sp. CLR",y=mtb_id)
  return(p)
}

#BAs
ids <- c("27729", "129_i", "25364", "29947", "25834", "686_i","27820","27393","26238","27872","807_i","25209")
#all60
ids<-c(bdm_fasted$FeatureID,bdm_fed$FeatureID)

for (id in ids) {
  #p <- spearmancorr_mmvec(lreu, mtb_shared, id)
  p <- spearmancorr_mmvec(lawsp2175, mtb_shared, id)
  #ggsave(filename = paste0("mmvec/spearmacorr_mmvec/lreu/BAs/", id, ".pdf"),plot = p,height = 4,width = 4)
  #ggsave(filename = paste0("mmvec/spearmacorr_mmvec/lreu/all60mtb/", id, ".pdf"),plot = p,height = 4,width = 4)
  #ggsave(filename = paste0("mmvec/spearmacorr_mmvec/lawsinosp2175/BAs/", id, ".pdf"),plot = p,height = 4,width = 4)
  ggsave(filename = paste0("mmvec/spearmacorr_mmvec/lawsinosp2175/all60mtb/", id, ".pdf"),plot = p,height = 4,width = 4)
}

#All prevotellas and limosilactobacillus and all BAs

#cond<- fread("mmvec/AD_summary_shared/conditionals/conditionals.tsv")%>%
#cond<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv")%>%
# cond<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv")%>%
#   dplyr::rename(FeatureID=V1)%>%
#   left_join(.,taxon,by="FeatureID")%>%
#   filter(Genus=="g__Prevotella"|Genus=="g__Limosilactobacillus")%>%
#   dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
#   column_to_rownames("Species")%>%
#   t()%>%as.data.frame()%>%
#   rownames_to_column("FeatureID")%>%
#   left_join(.,mtb_id,by="FeatureID")%>%
#   filter(subclass_abr %in% c("Bile acids, alcohols and derivatives","Amino acids, peptides, and analogues"))%>%
#   mutate(label=paste(FeatureID," ",signif(mn_mz,5)," ",Compound_Name))%>%
#   dplyr::select(label,everything(),-row.ID,-mn_mz,-Compound_Name,-Organism,-superclass,-class,-subclass,
#                 -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-subclass_abr,-FeatureID)
# write.table(cond, "mmvec/AD_summary_shared/mmvec_prevlimo_BA_condprob.txt",sep = "\t",row.names = FALSE, quote=FALSE)
# write.table(cond, "mmvec/AD_summary_shared_fasted/mmvec_prevlimo_BA_condprob_fasted.txt",sep = "\t",row.names = FALSE, quote=FALSE)
# write.table(cond, "mmvec/AD_summary_shared_fed/mmvec_prevlimo_BA_condprob_fed.txt",sep = "\t",row.names = FALSE, quote=FALSE)
# 
# cond_mat<-cond%>%column_to_rownames("label")%>%as.matrix()
# 
# p<-pheatmap(cond_mat,main="Prev and Limosilactobacillus with BAs")
# ggsave("mmvec/AD_summary_shared_fed/mmvec_heatmap_PrevLimosilac_wBA.pdf", plot=p,height=20, width=30)
# ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_PrevLimosilac_wBA.pdf", plot=p,height=20, width=30)
# ggsave("mmvec/AD_summary_shared/mmvec_heatmap_PrevLimosilac_wBA.pdf", plot=p,height=20, width=30)
