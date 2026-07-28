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
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"))%>%
  arrange(zt_time)%>%
  arrange(ATTRIBUTE_host_subject_id)%>%
  dplyr::select(sample_name,everything())
write.table(mtb_md, "mmvec/mmvec_mtb_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt")%>%
  dplyr::select(FeatureID, all_of(mtb_md$sampleid)) %>%
  filter(rowSums(select(., all_of(mtb_md$sampleid)) > 0) >= ceiling(0.10 * length(mtb_md$sampleid)))

names(mtb) <-c("FeatureID",mtb_md$sample_name)
write.table(mtb, "mmvec/mmvec_mtb_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#load mtg
mtg_md<-fread("jointRPCA/stool.metagenomic.metadata.with_amyloid.tsv")%>%
  dplyr::rename(sampleid=sample_name)%>%
  mutate(sample_name = sub("^[^.]+\\.", "", sampleid))%>%
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"))%>%
  arrange(zt_time)%>%
  arrange(host_subject_id)%>%
  dplyr::select(sample_name,everything())
write.table(mtg_md, "mmvec/mmvec_mtg_metadata.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv")%>%
  dplyr::rename(FeatureID='OTU ID')%>%
  dplyr::select(FeatureID,all_of(mtg_md$sampleid))%>%
  filter(rowSums(select(., all_of(mtg_md$sampleid)) > 0) >= ceiling(0.05 * length(mtg_md$sampleid)))%>%
  rename_with(~str_replace_all(., "14708.", ""))

write.table(mtg, "mmvec/mmvec_mtg_feattab.txt",sep = "\t",row.names = FALSE, quote=FALSE)

# mtg_original<-mtg
# 
# total_reads_before <- mtg_original %>% dplyr::select(all_of(mtg_md$sampleid)) %>% sum()
# total_reads_after  <- mtg %>% dplyr::select(all_of(mtg_md$sampleid)) %>% sum()

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
######################################################################
#create files for fasted and fed analysis

##fasted
mtb_md_fasted<-mtb_md%>%filter(fasted=="fasted")

mtb_fasted<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt")%>%
  dplyr::select(FeatureID, all_of(mtb_md_fasted$sampleid)) %>%
  filter(rowSums(select(., all_of(mtb_md_fasted$sampleid)) > 0) >= ceiling(0.10 * length(mtb_md_fasted$sampleid)))

names(mtb_fasted) <-c("FeatureID",mtb_md_fasted$sample_name)
write.table(mtb_fasted, "mmvec/mmvec_mtb_feattab_fasted.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg_md_fasted<-mtg_md%>%filter(fasted=="fasted")

mtg_fasted<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv")%>%
  dplyr::rename(FeatureID='OTU ID')%>%
  dplyr::select(FeatureID,all_of(mtg_md_fasted$sampleid))%>%
  filter(rowSums(select(., all_of(mtg_md_fasted$sampleid)) > 0) >= ceiling(0.05 * length(mtg_md_fasted$sampleid)))%>%
  rename_with(~str_replace_all(., "14708.", ""))

write.table(mtg_fasted, "mmvec/mmvec_mtg_feattab_fasted.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#shared metadata with training-test column
shared<-intersect(mtg_md_fasted$sample_name,mtb_md_fasted$sample_name)

md<-mtg_md_fasted%>%
  filter(sample_name %in% shared)%>%
  dplyr::select(sample_name, host_subject_id,genotype_diet,zt_collection_time,fasted)

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sample_name %in% testsamps$sample_name, "Test","Train"))
write.table(md, "mmvec/mmvec_comb_metadata_fasted.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just keep columns that are shared
mtg_shared_fasted <- mtg_fasted%>%
  dplyr::select(FeatureID,any_of(shared))
write.table(mtg_shared_fasted, "mmvec/mmvec_mtg_shared_feattab_fasted.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtb_shared_fasted <- mtb_fasted%>%
  dplyr::select(FeatureID,any_of(shared))
write.table(mtb_shared_fasted, "mmvec/mmvec_mtb_shared_feattab_fasted.txt",sep = "\t",row.names = FALSE, quote=FALSE)

##fed
mtb_md_fed<-mtb_md%>%filter(fasted=="fed")

mtb_fed<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt")%>%
  dplyr::select(FeatureID, all_of(mtb_md_fed$sampleid)) %>%
  filter(rowSums(select(., all_of(mtb_md_fed$sampleid)) > 0) >= ceiling(0.10 * length(mtb_md_fed$sampleid)))

names(mtb_fed) <-c("FeatureID",mtb_md_fed$sample_name)
write.table(mtb_fed, "mmvec/mmvec_mtb_feattab_fed.txt",sep = "\t",row.names = FALSE, quote=FALSE)

mtg_md_fed<-mtg_md%>%filter(fasted=="fed")

mtg_fed<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv")%>%
  dplyr::rename(FeatureID='OTU ID')%>%
  dplyr::select(FeatureID,all_of(mtg_md_fed$sampleid))%>%
  filter(rowSums(select(., all_of(mtg_md_fed$sampleid)) > 0) >= ceiling(0.05 * length(mtg_md_fed$sampleid)))%>%
  rename_with(~str_replace_all(., "14708.", ""))

write.table(mtg_fed, "mmvec/mmvec_mtg_feattab_fed.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#shared metadata with training-test column
shared<-intersect(mtg_md_fed$sample_name,mtb_md_fed$sample_name)

md<-mtg_md_fed%>%
  filter(sample_name %in% shared)%>%
  dplyr::select(sample_name, host_subject_id,genotype_diet,zt_collection_time,fasted)

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sample_name %in% testsamps$sample_name, "Test","Train"))
write.table(md, "mmvec/mmvec_comb_metadata_fed.txt",sep = "\t",row.names = FALSE, quote=FALSE)

#just keep columns that are shared
mtg_shared_fed <- mtg_fed%>%
  dplyr::select(FeatureID,any_of(shared))
write.table(mtg_shared_fed, "mmvec/mmvec_mtg_shared_feattab_fed.txt",sep = "\t",row.names = FALSE, quote=FALSE)
mtb_shared_fed <- mtb_fed%>%
  dplyr::select(FeatureID,any_of(shared))
write.table(mtb_shared_fed, "mmvec/mmvec_mtb_shared_feattab_fed.txt",sep = "\t",row.names = FALSE, quote=FALSE)

######################################################################
#taxon and mtb metadata

mtb_id <- fread("metabolomics//data_files/AD_GNPS_BA/nf_output/library/library_matches_cleanwmz.csv") %>%
  mutate(across(everything(), ~ str_extract(.x, "(?<!NA)\\b[^|]+\\b")))%>%
  mutate(across(everything(), ~ ifelse(.x == "NA" | .x == "N/A", NA, .x)))%>%
  mutate(Compound_Name = str_replace(Compound_Name, " \\(delta.*", ""),
         subclass = ifelse(str_detect(Compound_Name,"bile acid"),"Bile acids, alcohols and derivatives",subclass))
  
mtb_id_cleaned <- as.data.frame(lapply(mtb_id, function(x) gsub('"', '', x)))
write.table(mtb_id_cleaned, "mmvec/mmvec_mtb_feat_annotation.txt",sep = "\t",row.names = FALSE, quote=FALSE)

taxon<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.tsv")%>%
  dplyr::rename(FeatureID='OTU ID')%>%
  dplyr::select(FeatureID)%>%
  separate(FeatureID,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="; ", remove=FALSE)
taxon_cleaned <- as.data.frame(lapply(taxon, function(x) gsub('"', '', x)))
write.table(taxon, "mmvec/mmvec_mtg_feat_annotation.txt",sep = "\t",row.names = FALSE, quote=FALSE)

######################################################################
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
micro_int4<-c("s__Lawsonibacter sp000492175","s__Limosilactobacillus reuteri","s__Limosilactobacillus fermentum","s__Prevotella corporis")

#create conditional heatmap

#Microbes of interest

#cond<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv",header=TRUE)%>%
cond<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv")%>%
#cond<- fread("mmvec/AD_summary_shared/conditionals/conditionals.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  filter(Species %in% micro_int)%>%
  dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
  column_to_rownames("Species")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(!is.na(Compound_Name))%>%
  mutate(label=paste("mtb",FeatureID,"m/z", signif(as.numeric(mn_mz),7),"RT",signif(as.numeric(mn_rt),3), sep=" "))%>%
  dplyr::select(label,everything(),-mn_mz,-mn_rt,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-FeatureID)
#write.table(cond, "mmvec/AD_summary_shared/conditionals/conditionals_microint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(cond, "mmvec/AD_summary_shared_fed/conditionals/conditionals_microint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(cond, "mmvec/AD_summary_shared_fasted/conditionals/conditionals_microint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

cond_mat<-cond%>%column_to_rownames("label")%>%as.matrix()
p<-pheatmap(cond_mat,main="Microbes of Int all Mtbs")
#ggsave("mmvec/AD_summary_shared_fed/mmvec_heatmap_60microb_allmtb.pdf", plot=p,height=25, width=35)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_60microb_allmtb.pdf", plot=p,height=25, width=35)
#ggsave("mmvec/AD_summary_shared/mmvec_heatmap_60microb_allmtb.pdf", plot=p,height=25, width=35)

#add just mtb of interest
bdm_fed<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fed-sc10-cred1.25.txt")
bdm_fasted<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fasted-sc10-cred1.5.txt")

#cond<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv",header=TRUE)%>%
cond<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv",header=TRUE)%>%
#cond<- fread("mmvec/AD_summary_shared/conditionals/conditionals.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  filter(Species %in% micro_int)%>%
  dplyr::select(Species,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-FeatureID)%>%
  column_to_rownames("Species")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(FeatureID %in% c(bdm_fed$FeatureID,bdm_fasted$FeatureID))%>%
  mutate(label=paste(FeatureID,"m/z",signif(as.numeric(mn_mz),7),"RT",signif(as.numeric(mn_rt),3),Compound_Name,sep=" "))%>%
  dplyr::select(label,everything(),-mn_mz,-mn_rt,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-FeatureID)
#write.table(cond, "mmvec/AD_summary_shared/conditionals/conditionals_microint_mtbint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
write.table(cond, "mmvec/AD_summary_shared_fasted/conditionals/conditionals_microint_mtbint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)
#write.table(cond, "mmvec/AD_summary_shared_fed/conditionals/conditionals_microint_mtbint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

cond_mat<-cond%>%column_to_rownames("label")%>%as.matrix()
p<-pheatmap(cond_mat,main="microbe int vs. fed/fasted mtb")
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_60microb_fastedfedmtb.pdf", plot=p,height=12, width=25)
#ggsave("mmvec/AD_summary_shared_fed/mmvec_heatmap_60microb_fastedfedmtb.pdf", plot=p,height=12, width=25)
#ggsave("mmvec/AD_summary_shared/mmvec_heatmap_60microb_fastedfedmtb.pdf", plot=p,height=12, width=25)

#just show the fasted and fed results for 7 microbes of interest (subset 3)

cond_fasted<- fread("mmvec/AD_summary_shared_fasted/conditionals/conditionals.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  filter(Species %in% micro_int4)%>%
  dplyr::select(FeatureID,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-Species)%>%
  column_to_rownames("FeatureID")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(FeatureID %in% c(bdm_fed$FeatureID,bdm_fasted$FeatureID))%>%
  mutate(label=paste(FeatureID,"m/z",signif(as.numeric(mn_mz),7),"RT",signif(as.numeric(mn_rt),3),Compound_Name,sep=" "))%>%
  dplyr::select(label,everything(),-mn_mz,-mn_rt,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-FeatureID)
  # filter(!(label %in% c("13644 m/z 387.2017 RT 4.45 3,5,5-Trimethyl-4beta-hydroxy-4-[3-(beta-D-glucopyranosyloxy)-1-butenyl]-2-cyclohexene-1-one",
  #                       "30838 m/z 418.2958 RT 7.68 Candidate Tyrosine-C16:1")))
write.table(cond_fasted, "mmvec/AD_summary_shared_fasted/conditionals/conditionals_fasted_4micro_mtbint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

#cond_fasted[cond_fasted$label == "13644   387.2   3,5,5-Trimethyl-4beta-hydroxy-4-[3-(beta-D-glucopyranosyloxy)-1-butenyl]-2-cyclohexene-1-one", "s__Limosilactobacillus reuteri"] <- -4

cond_fed<- fread("mmvec/AD_summary_shared_fed/conditionals/conditionals.tsv",header=TRUE)%>%
  dplyr::rename(FeatureID=V1)%>%
  left_join(.,taxon,by="FeatureID")%>%
  filter(Species %in% micro_int4)%>%
  dplyr::select(FeatureID,everything(),-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-Species)%>%
  column_to_rownames("FeatureID")%>%
  t()%>%as.data.frame()%>%
  rownames_to_column("FeatureID")%>%
  left_join(.,mtb_id,by="FeatureID")%>%
  filter(FeatureID %in% c(bdm_fed$FeatureID,bdm_fasted$FeatureID))%>%
  mutate(label=paste(FeatureID,"m/z",signif(as.numeric(mn_mz),7),"RT",signif(as.numeric(mn_rt),3),Compound_Name,sep=" "))%>%
  dplyr::select(label,everything(),-mn_mz,-mn_rt,-Compound_Name,-Organism,-superclass,-class,-subclass,
                -npclassifier_superclass,-npclassifier_class,-npclassifier_pathway,-FeatureID)
# filter(!(label %in% c("13644 m/z 387.2017 RT 4.45 3,5,5-Trimethyl-4beta-hydroxy-4-[3-(beta-D-glucopyranosyloxy)-1-butenyl]-2-cyclohexene-1-one",
#        "30838 m/z 418.2958 RT 7.68 Candidate Tyrosine-C16:1")))
write.table(cond_fed, "mmvec/AD_summary_shared_fed/conditionals/conditionals_fed_4micro_mtbint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)


comb_fasted<-left_join(cond_fasted,cond_fed,by="label")
write.table(comb_fasted, "mmvec/AD_summary_shared_fasted/conditionals/conditionals_fastedfed_microint4_mtbint.tsv",sep = "\t",row.names = FALSE, quote=FALSE)

comb_fasted<-comb_fasted%>%column_to_rownames("label")%>%as.matrix()

p<-pheatmap(comb_fasted,main="microbe int vs. fastedfed mtb", cluster_cols=FALSE)
ggsave("mmvec/AD_summary_shared_fasted/mmvec_heatmap_4newmicrob_fastedfedcombmtb.pdf", plot=p,height=18, width=18)
######################################################################
mtg_shared<-fread("metagenomics/data_files/159013_wolR2_per_genome.feature-table.species.clr.csv")%>%
  dplyr::rename(FeatureID=`OTU ID`)
mtb_shared<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_log.csv")%>%
  dplyr::rename(FeatureID=row.ID)%>%
  dplyr::select(FeatureID,all_of(mtb_md$sampleid))
names(mtb_shared) <-c("FeatureID",mtb_md$sample_name)
md<-fread("mmvec/mmvec_comb_metadata.txt")

lreu <- mtg_shared %>%filter(str_detect(FeatureID, "Limosilactobacillus reuteri"))%>%
  gather(sample_name,counts_x, -FeatureID)

lferm <- mtg_shared %>%filter(str_detect(FeatureID, "Limosilactobacillus fermentum"))%>%
  gather(sample_name,counts_x, -FeatureID)

lawsp2175 <- mtg_shared %>%filter(str_detect(FeatureID, "Lawsonibacter sp000492175"))%>%
  gather(sample_name,counts_x, -FeatureID)

prevcorp <- mtg_shared %>%filter(str_detect(FeatureID, "Prevotella corporis"))%>%
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
    #scale_y_continuous(limits=c(-1,17))+
    #labs(x="Lreu CLR",y=mtb_id)
    #labs(x="Lferm CLR",y=mtb_id)
    labs(x="Pcorp CLR",y=mtb_id)
    #labs(x="Law sp. CLR",y=mtb_id)
  return(p)
}

#BAs
ids <- c("27351","27878","27766","25808","27195","28788","25834","26254",
         "23133","27820","28285","27729","25030","27069")
#all58
ids<-c(bdm_fasted$FeatureID,bdm_fed$FeatureID)

for (id in ids) {
 #p <- spearmancorr_mmvec(lreu, mtb_shared, id)
 #p <- spearmancorr_mmvec(lferm, mtb_shared, id)
 p <- spearmancorr_mmvec(prevcorp, mtb_shared, id)
 #p <- spearmancorr_mmvec(lawsp2175, mtb_shared, id)
  #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/lreu/BAint/", id, ".pdf"),plot = p,height = 4,width = 4)
  #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/lreu/all58mtb/", id, ".pdf"),plot = p,height = 4,width = 4)
 #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/lferm/BAint/", id, ".pdf"),plot = p,height = 4,width = 4)
 #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/lferm/all58mtb/", id, ".pdf"),plot = p,height = 4,width = 4)
 #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/pcorp/BAint/", id, ".pdf"),plot = p,height = 4,width = 4)
 ggsave(filename = paste0("mmvec/spearmancorr_mmvec/pcorp/all58mtb/", id, ".pdf"),plot = p,height = 4,width = 4)
  #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/lawsinosp2175/BAint/", id, ".pdf"),plot = p,height = 4,width = 4)
  #ggsave(filename = paste0("mmvec/spearmancorr_mmvec/lawsinosp2175/all58mtb/", id, ".pdf"),plot = p,height = 4,width = 4)
}

