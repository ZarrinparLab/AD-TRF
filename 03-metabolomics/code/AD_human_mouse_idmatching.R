setwd("~/scratch/AD_EM/metabolomics/data_files/ID_matching_human_mouse/")
library(tidyverse)
library(ggpubr)
library(UpSetR)
library(Hmisc)
# Function - flatten correlation matrix (from sthda website)
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
# Read clustering files from Molecular Network of tissues .mgfs
clustering_mn <- read.delim("nf_output/clustering/clusterinfo.tsv") %>% dplyr::select(-8)
colnames(clustering_mn) <- c("clusterIDX", "filename", "specIDX", "Scan", "mz", "charge", "rt")
annottions <- read.delim("nf_output/networking/clustersummary_with_network.tsv") %>% dplyr::select(1, 3, 23)
clustering_mn_match <- clustering_mn %>%
  dplyr::mutate(Feature_organ = case_when(str_detect(pattern = "beam", filename) ~ paste("human", Scan, sep = "_"),
                                          str_detect(pattern = "trf", filename) ~ paste("mouse", Scan, sep = "_")))
network_pairs <- read.delim("nf_output/networking/merged_pairs.tsv")
network_components_info <- read.delim("nf_output/networking/clustersummary_with_network.tsv") %>% dplyr::select(1,3,7,8)
##########################################
# Track features across multiple tissues #
##########################################
network_pairs_info <- network_pairs %>%
  dplyr::filter(abs(DeltaMZ) < 0.02) %>% # because run with orbitrap
  left_join(clustering_mn_match %>% dplyr::select(clusterIDX, Feature_organ, mz, rt, charge), by = c("CLUSTERID1" = "clusterIDX")) %>%
  left_join(clustering_mn_match %>% dplyr::select(clusterIDX, Feature_organ, mz, rt, charge), by = c("CLUSTERID2" = "clusterIDX")) %>%
  dplyr::mutate(DeltaRT = rt.x - rt.y) %>%
  #dplyr::filter(abs(DeltaRT) < 0.3) %>% # same method was used so RT should be very similar
  dplyr::filter(Cosine > 0.7) %>% # make sure molecules are the same by high cosine score
  dplyr::mutate(SameTissue = (gsub("_.*$", "", Feature_organ.x) == gsub("_.*$", "", Feature_organ.y))) %>%
  dplyr::filter(SameTissue == FALSE) %>% # I want to remove matches within the same tissue
  dplyr::select(CLUSTERID1, CLUSTERID2, mz.x, mz.y, DeltaMZ, rt.x, rt.y, DeltaRT, Cosine, Feature_organ.x, Feature_organ.y) %>%
  dplyr::mutate(Matching_Organs = paste(gsub("_.*$", "", Feature_organ.x), gsub("_.*$", "", Feature_organ.y), sep = "_")) %>%
  left_join(network_components_info %>% dplyr::select(`cluster.index`, component), by = c("CLUSTERID1" = "cluster.index")) %>%
  left_join(annottions %>% dplyr::select(`cluster.index`, Compound_Name), by = c("CLUSTERID1" = "cluster.index"))

write.table(network_pairs_info, file = "nf_output/clustering/IDmatching_human_mouse.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(network_pairs_info, file = "nf_output/clustering/IDmatching_all.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
