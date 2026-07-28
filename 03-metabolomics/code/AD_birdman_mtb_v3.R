setwd("~/scratch/AD_EM/")

library(tidyverse)
library(data.table)
library(qiime2R)
library(ggpubr)

library(data.table)
library(dplyr)
library(tidyr)

##################################################
# Filter a feature table by per-group prevalence

filter_by_group_prevalence <- function(feature_table,
                                       metadata,
                                       id_col = "row.ID",
                                       sample_col = "sampleid",
                                       group_col,
                                       min_prop = 0.15,
                                       min_count = NULL,
                                       verbose = TRUE) {
  
  sample_ids <- metadata[[sample_col]]
  sample_ids <- sample_ids[sample_ids %in% names(feature_table)]
  
  # group sizes -> per-group required count
  group_sizes <- metadata %>%
    filter(.data[[sample_col]] %in% sample_ids) %>%
    count(.data[[group_col]], name = "n_samples")
  
  group_thresholds <- group_sizes %>%
    mutate(threshold = pmax(ceiling(n_samples * min_prop),
                            if (is.null(min_count)) 0 else min_count))
  
  if (verbose) {
    cat("Group sizes and required prevalence thresholds:\n")
    print(group_thresholds)
  }
  
  # long format: feature x sample -> presence/absence, joined to group
  long_presence <- feature_table %>%
    select(all_of(id_col), all_of(sample_ids)) %>%
    pivot_longer(-all_of(id_col), names_to = sample_col, values_to = "value") %>%
    left_join(metadata %>% select(all_of(sample_col), all_of(group_col)), by = sample_col) %>%
    group_by(across(all_of(c(id_col, group_col)))) %>%
    summarise(n_present = sum(value > 0), .groups = "drop") %>%
    left_join(group_thresholds %>% select(all_of(group_col), threshold), by = group_col)
  
  # a feature must pass the threshold in EVERY group to be kept
  features_to_keep <- long_presence %>%
    group_by(across(all_of(id_col))) %>%
    summarise(passes_all_groups = all(n_present >= threshold), .groups = "drop") %>%
    filter(passes_all_groups) %>%
    pull(all_of(id_col))
  
  filtered <- feature_table %>% filter(.data[[id_col]] %in% features_to_keep)
  
  if (verbose) {
    cat("\nFeatures before filter:", nrow(feature_table), "\n")
    cat("Features after filter: ", nrow(filtered), "\n")
    cat("Features dropped:      ", nrow(feature_table) - nrow(filtered), "\n")
  }
  
  filtered
}
##################################################
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces")%>%
  mutate(fasted=ifelse(ATTRIBUTE_ZT_collection_time %in% c("ZT8","ZT12","ZT16"),"fasted","fed"))

md%>%group_by(ATTRIBUTE_geno_n_feedingcond,fasted)%>%summarise(n=n())

##################################################
#subset data to make pairwise comparisons (fecal)

##fecal Tg-TRF vs. Tg-ALF
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond!="NTg-ALF")

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(md,"metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-tgTRFtgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
md_train<-md%>%filter(Testing=="Train")

ADf_TRFALFcoll_mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")%>%
  select(`row.ID`, all_of(md$sampleid))

ADf_TRFALFcoll_mtb_filt <- filter_by_group_prevalence(
  feature_table = ADf_TRFALFcoll_mtb,
  metadata      = md_train,
  group_col     = "ATTRIBUTE_geno_n_feedingcond",
  min_prop      = 0.15
)
write.table(ADf_TRFALFcoll_mtb_filt,"metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#split by fed
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond!="NTg-ALF" & ATTRIBUTE_ZT_collection_time %in% c("ZT20","ZT0","ZT4"))

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(md,"metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-tgTRFtgALF-fed.txt",sep = "\t",row.names = FALSE,quote=FALSE)
md_train<-md%>%filter(Testing=="Train")

ADf_TRFALFcoll_mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")%>%
  select(`row.ID`, all_of(md$sampleid))

ADf_TRFALFcoll_mtb_filt <- filter_by_group_prevalence(
  feature_table = ADf_TRFALFcoll_mtb,
  metadata      = md_train,
  group_col     = "ATTRIBUTE_geno_n_feedingcond",
  min_prop      = 0.15
)
write.table(ADf_TRFALFcoll_mtb_filt,"metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fed.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#split by fasted
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond!="NTg-ALF" & ATTRIBUTE_ZT_collection_time %in% c("ZT8","ZT12","ZT16"))

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(md,"metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-tgTRFtgALF-fasted.txt",sep = "\t",row.names = FALSE,quote=FALSE)
md_train<-md%>%filter(Testing=="Train")

ADf_TRFALFcoll_mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")%>%
  select(`row.ID`, all_of(md$sampleid))

ADf_TRFALFcoll_mtb_filt <- filter_by_group_prevalence(
  feature_table = ADf_TRFALFcoll_mtb,
  metadata      = md_train,
  group_col     = "ATTRIBUTE_geno_n_feedingcond",
  min_prop      = 0.15,
  min_count     = 3
)
write.table(ADf_TRFALFcoll_mtb_filt,"metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fasted.txt",sep = "\t",row.names = FALSE,quote=FALSE)

##fecal Tg-ALF vs. NTg-ALF
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond!="Tg-TRF")

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(md,"metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-tgALFntgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)
md_train<-md%>%filter(Testing=="Train")

ADf_tgntgALFcorr_mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")%>%
  select(`row.ID`, all_of(md$sampleid))

ADf_tgntgALFcorr_mtb_filt <- filter_by_group_prevalence(
  feature_table = ADf_tgntgALFcorr_mtb,
  metadata      = md_train,
  group_col     = "ATTRIBUTE_geno_n_feedingcond",
  min_prop      = 0.15
)
write.table(ADf_tgntgALFcorr_mtb_filt,"metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10-fecal-tgALFntgALF.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just fed
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond!="Tg-TRF" & ATTRIBUTE_ZT_collection_time %in% c("ZT20","ZT0","ZT4"))

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(md,"metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-tgALFntgALF-fed.txt",sep = "\t",row.names = FALSE,quote=FALSE)
md_train<-md%>%filter(Testing=="Train")

ADf_tgntgALFcorr_mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")%>%
  select(`row.ID`, all_of(md$sampleid))

ADf_tgntgALFcorr_mtb_filt <- filter_by_group_prevalence(
  feature_table = ADf_tgntgALFcorr_mtb,
  metadata      = md_train,
  group_col     = "ATTRIBUTE_geno_n_feedingcond",
  min_prop      = 0.15
)
write.table(ADf_tgntgALFcorr_mtb_filt,"metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fed.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#just fasted
md<-fread("metabolomics/data_files/FBMN/metadata_table/metadata_table-noblanks.txt")%>%
  filter(ATTRIBUTE_BodyPartName=="feces"& ATTRIBUTE_geno_n_feedingcond!="Tg-TRF" & ATTRIBUTE_ZT_collection_time %in% c("ZT8","ZT12","ZT16"))

set.seed(1454)
testsamps<-slice_sample(md, prop = 0.15)
md<-md%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(md,"metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-tgALFntgALF-fasted.txt",sep = "\t",row.names = FALSE,quote=FALSE)
md_train<-md%>%filter(Testing=="Train")

ADf_tgntgALFcorr_mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")%>%
  select(`row.ID`, all_of(md$sampleid))

ADf_tgntgALFcorr_mtb_filt <- filter_by_group_prevalence(
  feature_table = ADf_tgntgALFcorr_mtb,
  metadata      = md_train,
  group_col     = "ATTRIBUTE_geno_n_feedingcond",
  min_prop      = 0.15
)
write.table(ADf_tgntgALFcorr_mtb_filt,"metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fasted.txt",sep = "\t",row.names = FALSE,quote=FALSE)

##########################################################

#reformated birdman results

bdm_reformat_alfs<-function(mtb,mtb_id,annot,subsetcred=TRUE){
    mtb<-mtb%>%
    dplyr::rename(ratio=`C(ATTRIBUTE_geno_n_feedingcond, Treatment('NTg-ALF'))[T.Tg-ALF]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(ATTRIBUTE_geno_n_feedingcond, Treatment('NTg-ALF'))[T.Tg-ALF]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_geno_n_feedingcond, Treatment('NTg-ALF'))[T.Tg-ALF]_hdi`))%>%
    separate(`C(ATTRIBUTE_geno_n_feedingcond, Treatment('NTg-ALF'))[T.Tg-ALF]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
    left_join(.,mtb_id, by="FeatureID")%>%
    mutate(name=ifelse(is.na(Compound_Name),
                       paste("mtb",FeatureID,"m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                       paste("mtb",FeatureID,"m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))
    
  if(subsetcred==TRUE){
    mtb<-mtb%>%filter(cred=="credible")
  }

  if(annot=="just_annot"){
    mtb<-mtb%>%filter(!is.na(Compound_Name))%>%
      select(FeatureID, ratio, min, max, cred, name,Compound_Name,superclass,class,subclass,npclassifier_superclass,
             npclassifier_class,npclassifier_pathway)%>%
      arrange(ratio)
    }
  else{
    mtb<-mtb%>%
      select(FeatureID, ratio, min, max, cred,name,Compound_Name,superclass,class,subclass,npclassifier_superclass,
             npclassifier_class,npclassifier_pathway)%>%
      arrange(ratio)
  }
  return(mtb)
}

bdm_reformat_trfs<-function(mtb,mtb_id,annot,subsetcred=TRUE){
  mtb<-mtb%>%
    dplyr::rename(ratio=`C(ATTRIBUTE_geno_n_feedingcond, Treatment('Tg-ALF'))[T.Tg-TRF]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(ATTRIBUTE_geno_n_feedingcond, Treatment('Tg-ALF'))[T.Tg-TRF]_hdi`=gsub("[(]|[)]","",`C(ATTRIBUTE_geno_n_feedingcond, Treatment('Tg-ALF'))[T.Tg-TRF]_hdi`))%>%
    separate(`C(ATTRIBUTE_geno_n_feedingcond, Treatment('Tg-ALF'))[T.Tg-TRF]_hdi`,c("min","max"), sep=",")%>%
    mutate(min=as.numeric(min),
           max=as.numeric(max),
           cred=ifelse(min>0|max<0,"credible","not_credible"))%>%
    left_join(.,mtb_id, by="FeatureID")%>%
    mutate(name=ifelse(is.na(Compound_Name),
                       paste("mtb",FeatureID,"m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                       paste("mtb",FeatureID,"m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))
  
  if(subsetcred==TRUE){
    mtb<-mtb%>%filter(cred=="credible")
  }
  
  if(annot=="just_annot"){
    mtb<-mtb%>%filter(!is.na(Compound_Name))%>%
      select(FeatureID, ratio, min, max, cred,name,Compound_Name,superclass,class,subclass,npclassifier_superclass,
             npclassifier_class,npclassifier_pathway)%>%
      arrange(ratio)
  }
  else{
    mtb<-mtb%>%
      select(FeatureID, ratio, min, max, cred,name,Compound_Name,superclass,class,subclass,npclassifier_superclass,
             npclassifier_class,npclassifier_pathway)%>%
      arrange(ratio)
  }
  return(mtb)
}

bdm_plot <- function(mtb){
  mtb$name <- factor(mtb$name, levels = unique(mtb$name))
  
  p <- ggplot(mtb, aes(x = name, y = ratio, ymin = min, ymax = max)) + 
    geom_linerange(position = position_dodge(width = 0.8)) + 
    theme_classic() +
    geom_pointrange(position = position_dodge(width = 0.8)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip()
  return(p)
}

mtb_id<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/library/library_matches_cleanwmz.csv")%>%
  mutate(Compound_Name = ifelse(grepl("^NA(?:\\|NA)+$", Compound_Name), NA, Compound_Name))
##########################################################
##NTg-ALF vs. Tg-ALF (stool)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"everything",FALSE)
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-tgALFntgALF-sc10-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"everything")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgALFntgALF-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool)+labs(title="NTg-ALF vs. Tg-ALF (stool)")
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgALFntgALF-sc10.pdf",height=13, width=23)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"just_annot")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgALFntgALF-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool)+labs(title="NTg-ALF vs. Tg-ALF (stool)")
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgALFntgALF-sc10.pdf",height=7, width=23)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"just_annot")
stool_alf<-stool%>%filter(min> 0.25|max< -0.25)
write.table(stool_alf,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgALFntgALF-sc10-cred0.25.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool_alf)+labs(title="NTg-ALF vs. Tg-ALF (stool)")
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgALFntgALF-sc10-cred0.25.pdf",height=3, width=15.5)

##NTg-ALF vs. Tg-ALF (stool-fasted)
stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fasted.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"everything",FALSE)
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-tgALFntgALF-fasted-sc10-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fasted.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"everything")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgALFntgALF-fasted-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fasted.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"just_annot")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgALFntgALF-fasted-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)

##NTg-ALF vs. Tg-ALF (stool-fed)
stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fed.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"everything",FALSE)
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-tgALFntgALF-fed-sc10-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fed.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"everything")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgALFntgALF-fed-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgALFntgALF-fed.beta_var.tsv")
stool<-bdm_reformat_alfs(stool,mtb_id,"just_annot")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgALFntgALF-fed-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#########################################################

##Tg-ALF vs. Tg-TRF (stool)
stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"everything",FALSE)
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-sc10-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"everything")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgTRFtgALF-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"just_annot")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)

##Tg-ALF vs. Tg-TRF (stool-fasted)
stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fasted.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"everything",FALSE)
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fasted-sc10-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fasted.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"everything")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgTRFtgALF-fasted-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool)+labs(title="Tg-ALF vs. Tg-TRF (stool fasted)")+
  scale_x_discrete(labels = function(name) str_wrap(name, width = 300))
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgTRFtgALF-fasted-sc10.pdf",height=45, width=22)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fasted.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"just_annot")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fasted-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool)+labs(title="Tg-ALF vs. Tg-TRF (stool fasted)")+
  scale_x_discrete(labels = function(name) str_wrap(name, width = 300))
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fasted-sc10.pdf",height=12, width=22)

stool_sub_fasted<-stool%>%filter(ratio< -1.5 | ratio>1.5)
write.table(stool_sub_fasted,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fasted-sc10-cred1.5.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool_sub_fasted)+labs(title="Tg-ALF vs. Tg-TRF (stool)")+
  scale_y_continuous(breaks = seq(-12, 2, by = 2))+
  scale_x_discrete(labels = function(name) str_wrap(name, width = 300))
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fasted-sc10-cred1.5.pdf",height=7, width=15)
#stool_fasted<-bdm_reformat_trfs(stool, mtb_id,"just_annot")

##Tg-ALF vs. Tg-TRF (stool-fed)
stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fed.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"everything",FALSE)
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fed-sc10-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fed.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"everything")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgTRFtgALF-fed-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool)+labs(title="Tg-ALF vs. Tg-TRF (stool fed)")+
  scale_x_discrete(labels = function(name) str_wrap(name, width = 300))
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcred-tgTRFtgALF-fed-sc10.pdf",height=45, width=22)

stool<-fread("metabolomics/birdman/rfm_cln/featuretable_reformated_clean_sc10-fecal-tgTRFtgALF-fed.beta_var.tsv")
stool<-bdm_reformat_trfs(stool, mtb_id,"just_annot")
write.table(stool,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fed-sc10.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool)+labs(title="Tg-ALF vs. Tg-TRF (stool fed)")+
  scale_x_discrete(labels = function(name) str_wrap(name, width = 300))
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fed-sc10.pdf",height=30, width=22)

stool_sub_fed<-stool%>%filter(ratio< -1.25 | ratio>1.25)
write.table(stool_sub_fed,"metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fed-sc10-cred1.25.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(stool_sub_fed)+labs(title="Tg-ALF vs. Tg-TRF (stool)")+
  scale_x_discrete(labels = function(name) str_wrap(name, width = 300))
ggsave("metabolomics/birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fed-sc10-cred1.25.pdf",height=7, width=16)
#stool_fed<-bdm_reformat_trfs(stool, mtb_id,"just_annot")

intersect(stool_sub_fed$FeatureID,stool_sub_fasted$FeatureID)#3026 3010

md<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/metadata/metadata_table-fecal-forrpca.txt")%>%
  mutate(fasted=ifelse(ATTRIBUTE_ZT_collection_time %in% c("ZT8","ZT12","ZT16"), "fasted","fed"))

md%>%group_by(ATTRIBUTE_geno_n_feedingcond,fasted)%>%
  summarise(n=n())

mtb_id<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/library/library_matches_cleanwmz.csv")%>%
  mutate(Compound_Name = ifelse(grepl("^NA(?:\\|NA)+$", Compound_Name), NA, Compound_Name))%>%
  mutate(Compound_Name = str_extract(Compound_Name,".*(?=\\(delta mass)"))

idmatching_bdm_signifA<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fastedfed_signif_inhumansserum_sub.txt")

mtb<-fread("metabolomics/data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_log.csv")%>%
  dplyr::rename(FeatureID=row.ID)%>%
  dplyr::select(FeatureID,all_of(md$sampleid))%>%
  #filter(FeatureID %in% stool_sub_fasted$FeatureID)%>%
  #filter(FeatureID %in% stool_sub_fed$FeatureID)%>%
  #filter(FeatureID %in% alf_res$FeatureID)%>%
  filter(FeatureID %in% idmatching_bdm_signifA$Feature_organ.mouse)%>%
  gather(sampleid,mtb_counts,-FeatureID) %>%
  left_join(.,mtb_id, by ="FeatureID") %>%
  mutate(label_name=ifelse(is.na(Compound_Name),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,digits=4), sep=" "),
                           paste("mtb",FeatureID, "m/z", signif(mn_mz,digits=4), Compound_Name, sep=" ")))%>%
  left_join(.,md,by="sampleid") %>%
  filter(ATTRIBUTE_BodyPartName=="feces")%>%
  mutate(ATTRIBUTE_geno_n_feedingcond=factor(ATTRIBUTE_geno_n_feedingcond, levels = c("NTg-ALF","Tg-ALF","Tg-TRF")))

bdm_indv_bxplt <- function(data,FeatureID) {
  
  lab<-unique((data%>%filter(FeatureID==FeatureID))$label_name)
  
  data%>%
    #ggplot(aes(x=ATTRIBUTE_geno_n_feedingcond, y=mtb_counts, fill=fasted)) + geom_boxplot(alpha=0.3) +
    ggplot(aes(x=ATTRIBUTE_geno_n_feedingcond, y=mtb_counts, fill=ATTRIBUTE_geno_n_feedingcond)) + geom_boxplot(alpha=0.3) +
    geom_jitter(width=0.25, aes(colour = ATTRIBUTE_geno_n_feedingcond))+
    #geom_jitter(aes(colour = fasted), position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.75))+
    facet_wrap(~fasted,scales="free_y")+
    scale_fill_manual(values=c("darkorange","purple4","forestgreen"))+
    scale_colour_manual(values=c("darkorange","purple4","forestgreen"))+
    theme_pubr()+ theme(legend.position = "top")+
    labs(color="genotype",
         fill="genotype",
         y ="log(peak abun)",
         x ="genotype",
         title = stringr::str_wrap(FeatureID, width=30))+
    theme(plot.title = element_text(face = "bold",hjust = 0.5))
  #ggsave(paste("metabolomics/birdman/rfm_cln/stool_tgTRFtgALF_bdmsc10_hits_logpeakabun_fasted/",FeatureID,".pdf",sep=""), height=4, width=5)
  #ggsave(paste("metabolomics/birdman/rfm_cln/stool_tgTRFtgALF_bdmsc10_hits_logpeakabun_fed/",FeatureID,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("metabolomics/birdman/rfm_cln/stool_tgTRFtgALF_bdmsc10_hits_logpeakabun_humanmatch/",FeatureID,".pdf",sep=""), height=4, width=4)
  #ggsave(paste("metabolomics/birdman/rfm_cln/stool_ntgALFtgALF_bdmsc10_hits_logpeakabun/",FeatureID,".pdf",sep=""), height=4, width=4)
}

mtbbdm_nested <- mtb %>% 
  group_by(FeatureID) %>% 
  nest()

mtbbdm_plots <- 
  mtbbdm_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ bdm_indv_bxplt(.x,.y)))

#cln plots for certain mtbs
mtb26254<-bdm_indv_bxplt(mtb%>%filter(FeatureID==26254),26254)+
  scale_y_continuous(limits=c(-1,15))+theme(legend.position = "none")
ggsave("metabolomics/birdman/rfm_cln/stool_ntgALFtgALF_bdmsc10_hits_logpeakabun/26254_cln.pdf", height=3.5, width=4)

mtb27477<-bdm_indv_bxplt(mtb%>%filter(FeatureID==27477),27477)+
  scale_y_continuous(limits=c(-1,15))+theme(legend.position = "none")
ggsave("metabolomics/birdman/rfm_cln/stool_ntgALFtgALF_bdmsc10_hits_logpeakabun/27477_cln.pdf", height=3.5, width=4)

mtb25030<-bdm_indv_bxplt(mtb%>%filter(FeatureID==25030),25030)+
  scale_y_continuous(limits=c(7,16))+theme(legend.position = "none")
ggsave("metabolomics/birdman/rfm_cln/stool_tgTRFtgALF_bdmsc10_hits_logpeakabun_humanmatch/25030_cln.pdf", height=3, width=4.5)

mtb27069<-bdm_indv_bxplt(mtb%>%filter(FeatureID==27069),27069)+
  scale_y_continuous(limits=c(8.5,18))+theme(legend.position = "none")
ggsave("metabolomics/birdman/rfm_cln/stool_tgTRFtgALF_bdmsc10_hits_logpeakabun_humanmatch/27069_cln.pdf", height=3, width=4.5)

##################################################
#get summary table of birdman results

mzs<-mtb_id%>%dplyr::select(1:3)

#combined 
ntgtgalf<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgALFntgALF-sc10-all.txt")%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)%>%
  mutate(comparison="NTg-ALF_v_Tg-ALF")

ntgalftrf<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-sc10-all.txt")%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)%>%
  mutate(comparison="Tg-ALF_v_Tg-TRF")

bdm_comm_summ<-rbind(ntgtgalf,ntgalftrf)%>%
  arrange(ratio)
write.table(bdm_comm_summ,"metabolomics/birdman/rfm_cln/birdman_mtb_comb_res.txt",sep = "\t",row.names = FALSE,quote=FALSE) #Table S5


#by fasted/fed
ntgtgalfl<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgALFntgALF-fasted-sc10-all.txt")%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)%>%
  mutate(comparison="NTg-ALF_v_Tg-ALF",
         feeding_state="fasted")

ntgalftrfl<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fasted-sc10-all.txt")%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)%>%
  mutate(comparison="Tg-ALF_v_Tg-TRF",
         feeding_state="fasted")

ntgtgalfd<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgALFntgALF-fed-sc10-all.txt")%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)%>%
  mutate(comparison="NTg-ALF_v_Tg-ALF",
         feeding_state="fed")

ntgalftrfd<-fread("metabolomics/birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fed-sc10-all.txt")%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)%>%
  mutate(comparison="Tg-ALF_v_Tg-TRF",
         feeding_state="fed")

bdm_comm_ff_summ<-rbind(ntgtgalfl,ntgalftrfl,ntgtgalfd,ntgalftrfd)%>%
  arrange(ratio)
write.table(bdm_comm_ff_summ,"metabolomics/birdman/rfm_cln/birdman_mtb_ff_comb_res.txt",sep = "\t",row.names = FALSE,quote=FALSE) #Table S7
