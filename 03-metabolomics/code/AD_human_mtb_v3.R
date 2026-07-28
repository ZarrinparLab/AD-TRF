setwd("~/scratch/AD_EM/metabolomics/")

library(tidyverse)
library(data.table)
library(qiime2R)
library(ggpubr)

#################################################
filter_by_group_prevalence <- function(feature_table,
                                       metadata,
                                       id_col = "FeatureID",
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
#################################################

mtb<-fread("data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10.csv")%>%
  dplyr::select(-contains("Blank"),-contains("Mix"))%>%
  rename_with(~ gsub("^plate\\d+_", "", .))%>%
  dplyr::rename(FeatureID=row.ID)
write.table(mtb, file = "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10_noblanks.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#pre-a
mdA<-fread("data_files/human_BEAM_GNPS_BA/nf_output/metadata/BEAM_metadata_KAW.csv")%>%
  filter(Visit.Type.Stool=="Pre 1" & Sample.Type == "Stool")%>%
  dplyr::select(filename,everything())%>%
  dplyr::rename(sampleid=filename)%>%
  mutate(MemDx_cln = gsub("-.*$", "", MemDx))

set.seed(2346)
testsamps<-slice_sample(mdA, prop = 0.15)
mdA<-mdA%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(mdA,"data_files/human_BEAM_GNPS_BA/nf_output/metadata/BEAM_metadata_KAW_preAstool_v2_bdm.txt",sep = "\t",row.names = FALSE,quote=FALSE)
mdA_train <-mdA%>%filter(Testing=="Train")
mtb_A<-fread("data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10_noblanks.txt")%>%
  select(FeatureID, all_of(mdA$sampleid))

mtb_A_filt <- filter_by_group_prevalence(
  feature_table = mtb_A,
  metadata      = mdA_train,
  group_col     = "MemDx_cln",
  min_prop      = 0.20
)

write.table(mtb_A_filt,"data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10_preA.txt",sep = "\t",row.names = FALSE,quote=FALSE)

#pre-b
mdB<-fread("data_files/human_BEAM_GNPS_BA/nf_output/metadata/BEAM_metadata_KAW.csv")%>%
  filter(Visit.Type.Stool=="Pre 2" & Sample.Type == "Stool")%>%
  dplyr::select(filename,everything())%>%
  dplyr::rename(sampleid=filename)%>%
  mutate(MemDx_cln = gsub("-.*$", "", MemDx))

set.seed(1454)
testsamps<-slice_sample(mdB, prop = 0.15)
mdB<-mdB%>%mutate(Testing=ifelse(sampleid %in% testsamps$sampleid, "Test","Train"))
write.table(mdB,"data_files/human_BEAM_GNPS_BA/nf_output/metadata/BEAM_metadata_KAW_preBstool_v2_bdm.txt",sep = "\t",row.names = FALSE,quote=FALSE)
mdB_train <-mdB%>%filter(Testing=="Train")

mtb_B<-fread("data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10_noblanks.txt")%>%
  select(FeatureID, all_of(mdB$sampleid))

mtb_B_filt <- filter_by_group_prevalence(
  feature_table = mtb_B,
  metadata      = mdB_train,
  group_col     = "MemDx_cln",
  min_prop      = 0.20
)
write.table(mtb_B_filt,"data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10_noblanks_preB.txt",sep = "\t",row.names = FALSE,quote=FALSE)
#################################################
#birdman results

bdm_reformat<-function(mtb,mtb_id,annot,subsetcred=TRUE){
  mtb<-mtb%>%
    dplyr::rename(ratio=`C(MemDx_cln, Treatment('NC'))[T.MCI]_mean`,
                  FeatureID=Feature)%>%
    mutate(`C(MemDx_cln, Treatment('NC'))[T.MCI]_hdi`=gsub("[(]|[)]","",`C(MemDx_cln, Treatment('NC'))[T.MCI]_hdi`),
           FeatureID=as.character(FeatureID))%>%
    separate(`C(MemDx_cln, Treatment('NC'))[T.MCI]_hdi`,c("min","max"), sep=",")%>%
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
      select(FeatureID, ratio, min, max, cred, name,Compound_Name,superclass,class,subclass,npclassifier_superclass,
             npclassifier_class,npclassifier_pathway)%>%
      arrange(ratio)
  }
  return(mtb)
}

bdm_plot<-function(mtb){
  mtb$name <- factor(mtb$name,levels = mtb$name)
  
  p<-ggplot(mtb, aes(x =name , y = ratio, ymin = min, ymax = max)) + 
    geom_linerange(position = position_dodge(width = 0.8)) + theme_classic()+
    geom_pointrange(position = position_dodge(width = 0.8)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip()
  return(p)
}


mtb_id<-fread("data_files/human_BEAM_GNPS_BA/nf_output/library/BEAM_library_matches_clean.csv")%>%
  mutate(Compound_Name = ifelse(grepl("^NA(?:\\|NA)+$", Compound_Name), NA, Compound_Name))%>%
  mutate(FeatureID=as.character(FeatureID))

a<-fread("birdman/human_BEAM/BEAM_gnps_quant_clean_sc10_preA.beta_var.tsv")
a<-bdm_reformat(a,mtb_id,"everything",FALSE)
write.table(a,"birdman/human_BEAM/fecal-bdm-NCMCIpreAstool-all.txt",sep = "\t",row.names = FALSE,quote=FALSE)

a<-fread("birdman/human_BEAM/BEAM_gnps_quant_clean_sc10_preA.beta_var.tsv")
a<-bdm_reformat(a,mtb_id,"everything")
write.table(a,"birdman/human_BEAM/fecal-bdm-justcred-NCMCIpreAstool.txt",sep = "\t",row.names = FALSE,quote=FALSE)
bdm_plot(a)+labs(title="BEAM_preA")

a<-fread("birdman/human_BEAM/BEAM_gnps_quant_clean_sc10_preA.beta_var.tsv")
a<-bdm_reformat(a,mtb_id,"just_annot")
bdm_plot(a)+labs(title="BEAM_preA stool")
ggsave("birdman/human_BEAM//fecal-bdm-justcredannot-NCMCIpreAstool.pdf",height=10, width=23)

a_sub<-a%>%filter(min> 0.25|max< -0.25)
bdm_plot(a_sub)+labs(title="BEAM_preA stool")
ggsave("birdman/human_BEAM/fecal-bdm-justcredannot-NCMCIpreAstool-cred0.25.pdf",height=6, width=23.5)

#summary table
mzs<-mtb_id%>%dplyr::select(1:3)
a<-fread("birdman/human_BEAM/fecal-bdm-NCMCIpreAstool-all.txt")%>%
  mutate(FeatureID=as.character(FeatureID))%>%
  left_join(.,mzs,by="FeatureID")%>%
  dplyr::select(FeatureID,mn_mz,mn_rt,Compound_Name,ratio,min,max,cred)
write.table(a,"birdman/human_BEAM/fecal-bdm-human-res.txt",sep = "\t",row.names = FALSE,quote=FALSE) #Table S8
##################################################

#plot metabolites shared between mouse and humans
bdm_fasted<-fread("birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fasted-sc10.txt")
bdm_fed<-fread("birdman/rfm_cln/fecal-bdm-justcredannot-tgTRFtgALF-fed-sc10.txt")

idmatching<-fread("data_files/ID_matching_human_mouse/nf_output/clustering/IDmatching_human_mouse.tsv")%>%
  mutate(Feature_organ.x=gsub("human_","",Feature_organ.x),
         Feature_organ.y=gsub("mouse_","",Feature_organ.y))

idmatching_bdm_signifA<-idmatching%>%
  dplyr::rename(Feature_organ.mouse=Feature_organ.y,
                Feature_organ.human=Feature_organ.x)%>%
  filter(Feature_organ.mouse %in% c(bdm_fasted$FeatureID,bdm_fed$FeatureID))%>%
  #filter(Feature_organ.human %in% a$FeatureID)
  filter(Feature_organ.human %in% a_sub$FeatureID)
#write.table(idmatching_bdm_signifA, file = "birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fastedfed_signif_inhumansserum_all.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(idmatching_bdm_signifA, file = "birdman/rfm_cln/fecal-bdm-tgTRFtgALF-fastedfed_signif_inhumansserum_sub.txt", sep = "\t", row.names = FALSE, quote = FALSE)

md<-fread("data_files/human_BEAM_GNPS_BA/nf_output/metadata/BEAM_metadata_KAW_preA_bdm.txt")%>%
  mutate(MemDx=gsub("-.*$","",MemDx))%>%
  dplyr::rename(filename=sampleid)

mtb_id<-fread("data_files/human_BEAM_GNPS_BA/nf_output/library/BEAM_library_matches_clean.csv")%>%
  mutate(Compound_Name = ifelse(grepl("^NA(?:\\|NA)+$", Compound_Name), NA, Compound_Name))

mtb_log<-fread("data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_log.csv")%>%
  dplyr::select(-contains("Blank"),-contains("Mix"))%>%
  rename_with(~ gsub("^plate\\d+_", "", .))%>%
  dplyr::rename(FeatureID=row.ID)
write.table(mtb_log, file = "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_log_noblanks.txt", sep = "\t", row.names = FALSE, quote = FALSE)

mtbA<-fread("data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_log_noblanks.txt")%>%
  dplyr::select(FeatureID, all_of(md$filename))%>%
  filter((FeatureID %in% idmatching_bdm_signifA$Feature_organ.human))%>%
  gather(filename,mtb_counts,-FeatureID) %>%
  left_join(.,mtb_id, by ="FeatureID")%>%
  mutate(label_name=ifelse(is.na(Compound_Name),
                     paste("mtb",FeatureID,"m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), sep=" "),
                     paste("mtb",FeatureID,"m/z", signif(mn_mz,7),"RT",signif(mn_rt,3), Compound_Name, sep=" ")))%>%
  right_join(.,md,by="filename")%>%
  filter(Visit.Type.Stool=="Pre 1" & Sample.Type == "Stool")%>%
  mutate(MemDx=factor(MemDx,levels=c("NC","MCI")))

bdm_indv_bxplt <- function(data,FeatureID) {
  
  lab<-unique((data%>%filter(FeatureID==FeatureID))$label_name)
  
  data%>%
    ggplot(aes(x=MemDx, y=mtb_counts, fill=MemDx)) + geom_boxplot(alpha=0.3) +
    geom_jitter(width=0.25, aes(colour = MemDx))+
    scale_fill_manual(values=c("darkorange","purple4"))+
    scale_colour_manual(values=c("darkorange","purple4"))+
    #facet_grid(~Sample.Type)+
    theme_pubr()+ theme(legend.position = "top")+
    labs(color="cognitive status",
         fill="cognitive status",
         y ="log(peak abun)",
         x ="cognitive status",
         title = stringr::str_wrap(FeatureID, width=30))+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
         #title = stringr::str_wrap(lab, width=30))+ theme(plot.title = element_text(face = "bold",hjust = 0.5))
  ggsave(paste("birdman/human_BEAM/match_mtbs_plot_preAdatastool/",FeatureID,".pdf",sep=""), height=3.5, width=3)
}

mtbbdm_nested <- mtbA %>% 
  group_by(FeatureID) %>% 
  nest()

mtbbdm_plots <- 
  mtbbdm_nested %>% 
  mutate(plot = map2(data, FeatureID,  ~ bdm_indv_bxplt(.x,.y)))
