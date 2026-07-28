setwd("~/scratch/AD_EM/metabolomics/")

library(tidyverse)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

## Cleanup feature quant table from MZmine ###########################################################################
#mouse

#norm <- read.csv("data_files/FBMN_simone/gnps_quant.csv")
norm <- read.csv("data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")
#norm <- norm %>% dplyr::select(-X)
norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ".mzML Peak area"))
#norm <- norm %>% dplyr::select(-(row.m.z:row.CCS), -(annotation.network.number:neutral.M.mass))
norm <- norm %>% dplyr::select(-(row.m.z:ccs), -(annotation.network.number:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID:correlation.group.ID, as.character))

## Merge ions from the same molecule (iimn)
iimn <- norm %>%
  group_by(correlation.group.ID) %>%
  summarise(across(where(is.numeric), ~sum(., na.rm = TRUE)),
            collapserowID = paste(row.ID, collapse = '|')) %>%
  filter(!is.na(correlation.group.ID)) %>%
  rename(row.ID = correlation.group.ID)%>%
  dplyr::select(row.ID,collapserowID,everything())
labels <- norm %>% dplyr::select(row.ID, correlation.group.ID)
iimn$row.ID <- paste0(iimn$row.ID, "_i")
norm_iimn <- norm %>% filter(is.na(correlation.group.ID)) %>% dplyr::select(-correlation.group.ID)%>%
  mutate(collapserowID=NA)%>%
  dplyr::select(row.ID,collapserowID,everything())%>%
  dplyr::select(-Unnamed..250)
norm_iimn <- rbind(norm_iimn, iimn)

#setdiff(colnames(norm_iimn),colnames((iimn)))

write_csv(norm_iimn, "data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_collapseiin_sum_corrgroup.csv")

norm_iimn_fornorm<-norm_iimn%>%dplyr::select(-collapserowID)%>%
  column_to_rownames("row.ID")%>%as.matrix()

norm_iimn_rclr <- norm_iimn_fornorm %>% decostand(method = "rclr")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_rclr, "data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_collapseiin_sum_corrgroup_rclr.csv")

# norm_iimn_alr <- norm_iimn_fornorm %>% decostand(method = "alr")%>%
#   as.data.frame()%>%rownames_to_column("row.ID")
# 
# norm_iimn_clr <- norm_iimn_fornorm %>% decostand(method = "clr",pseudocount = 1)%>%
#   as.data.frame()%>%rownames_to_column("row.ID")

norm_iimn_log <- norm_iimn_fornorm %>% decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_log, "data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_collapseiin_sum_corrgroup_log.csv")

norm_iimn_tic <- decostand(norm_iimn_fornorm, method = "total") * median(rowSums(norm_iimn_fornorm))
norm_iimn_log2 <- decostand(norm_iimn_tic, method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_log2, "data_files/AD_GNPS_BA/nf_output/clustering/gnps_quant_collapseiin_sum_corrgroup_logtic.csv")

## Read in FBMN library IDs from GNPS
mzs <- read.csv("data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")%>%
  dplyr::select(1:3)%>%
  mutate(row.ID=as.character(row.ID))


#library_ID <- read.delim("data_files/FBMN_simone/fbmn_wcanopus.txt")%>%
library_ID <- read.delim("data_files/AD_GNPS_BA/nf_output/library/merged_results_with_gnps.tsv")%>%
  #dplyr::rename(row.ID=FeatureID)
  dplyr::rename(row.ID=`X.Scan.`)
library_ID$row.ID <- as.character(library_ID$row.ID)
library_ID <- left_join(mzs, library_ID, by = "row.ID", multiple = "all")
library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")

result_list_iin <- library_ID %>%
  filter(!is.na(correlation.group.ID)) %>%
  group_by(correlation.group.ID) %>%
  summarise(
    mn_mz = mean(row.m.z, na.rm = TRUE),
    mn_rt = mean(row.retention.time, na.rm = TRUE),
    across(-row.m.z, ~paste(., collapse = "|")))%>%
  dplyr::rename(FeatureID=correlation.group.ID)%>%
  dplyr::select(c(1,4:2,19,36,44:49))

result_list_iin$FeatureID <- paste0(result_list_iin$FeatureID, "_i")

result_list <- library_ID %>%
  filter(is.na(correlation.group.ID)) %>%
  #dplyr::select(c(1,4:2,5:10))%>%
  dplyr::select(c(1,4:2,18,35,43:48))%>%
  dplyr::rename(FeatureID=row.ID,x=correlation.group.ID, mn_mz=row.m.z,mn_rt=row.retention.time)%>%
  dplyr::rename(row.ID=x)

result_list<-rbind(result_list_iin,result_list)
write_csv(result_list, "data_files/AD_GNPS_BA/nf_output/library/library_matches_iin_corrgroup.csv")

#note manually removing some annotations based on repeats
result_list_cln<-result_list%>%
  mutate(Compound_Name=ifelse(FeatureID=="612_i","NA",Compound_Name))%>%
  mutate(Compound_Name=ifelse(FeatureID=="58_i","NA",Compound_Name))
write_csv(result_list_cln, "data_files/AD_GNPS_BA/nf_output/library/library_matches_iin_corrgroup_cln.csv")

## Cleanup feature quant table from MZmine ###########################################################################
#human

norm <- read.csv("data_files/human_BEAM_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")
norm <- norm %>% rename_with(~str_replace_all(., ".mzXML.Peak.area", ".mzXML"))
norm <- norm %>% dplyr::select(-(row.m.z:ccs), -(annotation.network.number:neutral.M.mass))
#norm <- norm %>% dplyr::select(-(row.m.z:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID:correlation.group.ID, as.character))
#norm <- norm %>% mutate(across(row.ID, as.character))
write_csv(norm, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant.csv")

## Merge ions from the same molecule (iimn)
iimn <- norm %>%
  group_by(correlation.group.ID) %>%
  summarise(across(where(is.numeric), ~sum(., na.rm = TRUE)),
            collapserowID = paste(row.ID, collapse = '|')) %>%
  filter(!is.na(correlation.group.ID)) %>%
  rename(row.ID = correlation.group.ID)%>%
  dplyr::select(row.ID,collapserowID,everything())
labels <- norm %>% dplyr::select(row.ID, correlation.group.ID)
iimn$row.ID <- paste0(iimn$row.ID, "_i")
norm_iimn <- norm %>% filter(is.na(correlation.group.ID)) %>% dplyr::select(-correlation.group.ID)%>%
  mutate(collapserowID=NA)%>%
  dplyr::select(row.ID,collapserowID,everything())%>%
  dplyr::select(-Unnamed..431)
norm_iimn <- rbind(norm_iimn, iimn)

#setdiff(colnames(norm_iimn),colnames((iimn)))

write_csv(norm_iimn, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quantcollapseiin_sum_corrgroup.csv")

norm_iimn_fornorm<-norm_iimn%>%dplyr::select(-collapserowID)%>%
  column_to_rownames("row.ID")%>%as.matrix()

norm_iimn_rclr <- norm_iimn_fornorm %>% decostand(method = "rclr")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_rclr, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_collapseiin_sum_corrgroup_rclr.csv")

# norm_iimn_alr <- norm_iimn_fornorm %>% decostand(method = "alr")%>%
#   as.data.frame()%>%rownames_to_column("row.ID")
# 
# norm_iimn_clr <- norm_iimn_fornorm %>% decostand(method = "clr",pseudocount = 1)%>%
#   as.data.frame()%>%rownames_to_column("row.ID")

norm_log <- norm %>% column_to_rownames("row.ID")%>%decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_log, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_log.csv")

norm_iimn_log <- norm_iimn_fornorm %>% decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_log, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_collapseiin_sum_corrgroup_log.csv")

norm_iimn_tic <- decostand(norm_iimn_fornorm, method = "total") * median(rowSums(norm_iimn_fornorm))
norm_iimn_log2 <- decostand(norm_iimn_tic, method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_iimn_log2, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_collapseiin_sum_corrgroup_logtic.csv")

## Read in FBMN library IDs from GNPS
mzs <- read.csv("data_files/human_BEAM_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")%>%
  dplyr::select(1:3)%>%
  mutate(row.ID=as.character(row.ID))

library_ID <- read.delim("data_files/human_BEAM_GNPS_BA/nf_output/library/merged_results_with_gnps.tsv")%>%
  dplyr::rename(row.ID=`X.Scan.`)
library_ID$row.ID <- as.character(library_ID$row.ID)
library_ID <- left_join(mzs, library_ID, by = "row.ID", multiple = "all")
library_ID <- right_join(labels, library_ID, by = "row.ID", multiple = "all")

result_list_iin <- library_ID %>%
  filter(!is.na(correlation.group.ID)) %>%
  group_by(correlation.group.ID) %>%
  summarise(
    mn_mz = mean(row.m.z, na.rm = TRUE),
    mn_rt = mean(row.retention.time, na.rm = TRUE),
    across(-row.m.z, ~paste(., collapse = "|")))%>%
  dplyr::rename(FeatureID=correlation.group.ID)%>%
  dplyr::select(c(1,4:2,19,36,44:49))

result_list_iin$FeatureID <- paste0(result_list_iin$FeatureID, "_i")

result_list <- library_ID %>%
  filter(is.na(correlation.group.ID)) %>%
  dplyr::select(c(1,4:2,18,35,43:48))%>%
  dplyr::rename(FeatureID=row.ID,x=correlation.group.ID, mn_mz=row.m.z,mn_rt=row.retention.time)%>%
  dplyr::rename(row.ID=x)

result_list<-rbind(result_list_iin,result_list)
write_csv(result_list, "data_files/human_BEAM_GNPS_BA/nf_output/library/BEAM_library_matches_iin_corrgroup.csv")
