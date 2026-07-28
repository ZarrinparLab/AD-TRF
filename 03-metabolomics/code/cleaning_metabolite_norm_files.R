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

norm <- read.csv("data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")
norm <- norm %>% rename_with(~str_replace_all(., ".mzML.Peak.area", ".mzML Peak area"))
norm <- norm %>% dplyr::select(-(row.m.z:ccs), -(correlation.group.ID:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID, as.character))

write_csv(norm, "data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean.csv")

id_col <- 1
sample_cols <- if (!is.null(id_col)) colnames(norm)[-id_col] else colnames(norm)

norm_10   <- norm
norm_100  <- norm
norm_1000 <- norm

norm_10[sample_cols]   <- norm[sample_cols]   / 10
norm_100[sample_cols]  <- norm[sample_cols]   / 100
norm_1000[sample_cols] <- norm[sample_cols]   / 1000
write_csv(norm_10, "data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc10.csv")
write_csv(norm_100, "data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc100.csv")
write_csv(norm_1000, "data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_sc1000.csv")


norm_log<-norm%>%column_to_rownames("row.ID")%>%as.matrix()

norm_log <- norm_log %>% decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_log, "data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_log.csv")

## Read in FBMN library IDs from GNPS
mzs <- read.csv("data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")%>%
  dplyr::select(1:3)%>%
  mutate(row.ID=as.character(row.ID))

library_ID <- read.delim("data_files/AD_GNPS_BA/nf_output/library/merged_results_with_gnps.tsv")%>%
  dplyr::rename(row.ID=`X.Scan.`)
library_ID$row.ID <- as.character(library_ID$row.ID)
library_ID <- left_join(mzs, library_ID, by = "row.ID", multiple = "all")


result_list <- library_ID %>%
  dplyr::select(c(1:3,17,34,42:47))%>%
  dplyr::rename(FeatureID=row.ID, mn_mz=row.m.z,mn_rt=row.retention.time)

write_csv(result_list, "data_files/AD_GNPS_BA/nf_output/library/library_matches_cleanwmz.csv")

## Cleanup feature quant table from MZmine ###########################################################################
#human

norm <- read.csv("data_files/human_BEAM_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")
norm <- norm %>% rename_with(~str_replace_all(., ".mzXML.Peak.area", ".mzXML"))
norm <- norm %>% dplyr::select(-(row.m.z:ccs), -(correlation.group.ID:neutral.M.mass))
norm <- norm %>% mutate(across(row.ID, as.character))
write_csv(norm, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean.csv")

id_col <- 1
sample_cols <- if (!is.null(id_col)) colnames(norm)[-id_col] else colnames(norm)

norm_10   <- norm
norm_100  <- norm
norm_1000 <- norm

norm_10[sample_cols]   <- norm[sample_cols]   / 10
norm_100[sample_cols]  <- norm[sample_cols]   / 100
norm_1000[sample_cols] <- norm[sample_cols]   / 1000
write_csv(norm_10, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc10.csv")
write_csv(norm_100, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc100.csv")
write_csv(norm_1000, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_sc1000.csv")

norm_log<-norm%>%
  column_to_rownames("row.ID")%>%as.matrix()

norm_log <- norm %>% column_to_rownames("row.ID")%>%decostand(method = "log")%>%
  as.data.frame()%>%rownames_to_column("row.ID")
write_csv(norm_log, "data_files/human_BEAM_GNPS_BA/nf_output/clustering/BEAM_gnps_quant_clean_log.csv")

## Read in FBMN library IDs from GNPS
mzs <- read.csv("data_files/human_BEAM_GNPS_BA/nf_output/clustering/featuretable_reformated.csv")%>%
  dplyr::select(1:3)%>%
  mutate(row.ID=as.character(row.ID))

library_ID <- read.delim("data_files/human_BEAM_GNPS_BA/nf_output/library/merged_results_with_gnps.tsv")%>%
  dplyr::rename(row.ID=`X.Scan.`)
library_ID$row.ID <- as.character(library_ID$row.ID)
library_ID <- left_join(mzs, library_ID, by = "row.ID", multiple = "all")

result_list <- library_ID %>%
  dplyr::select(c(1:3,17,34,42:47))%>%
  dplyr::rename(FeatureID=row.ID, mn_mz=row.m.z,mn_rt=row.retention.time)

write_csv(result_list, "data_files/human_BEAM_GNPS_BA/nf_output/library/BEAM_library_matches_clean.csv")
