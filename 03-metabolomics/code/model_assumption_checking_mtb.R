setwd("~/scratch/AD_EM/metabolomics//")

library(data.table)
library(tidyverse)
library(ggpubr)

#############################################################
#inputs

dt<-"data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_feces.txt"
dt2<-"data_files/AD_GNPS_BA/nf_output/clustering/featuretable_reformated_clean_log.csv"
metat<-"data_files/FBMN/metadata_table/metadata_table-noblanks.txt"
#############################################################

mtb<-fread(dt)%>%
  column_to_rownames("FeatureID")


md<-fread(metat)%>%
  filter(ATTRIBUTE_BodyPartName=="feces")%>%
  mutate(ATTRIBUTE_sex=ifelse(ATTRIBUTE_sex=="F","female",ATTRIBUTE_sex),
         genofedZT=paste(ATTRIBUTE_ZT_collection_time, ATTRIBUTE_geno_n_feedingcond, sep = "_"))%>%
  mutate(ATTRIBUTE_timepoint = as.numeric(sub("^ZT", "", ATTRIBUTE_ZT_collection_time)))


prevalence_df <- md %>%
  filter(sampleid %in% colnames(mtb)) %>%
  group_by(ATTRIBUTE_geno_n_feedingcond) %>%
  group_modify(~ {
    samps <- .x$sampleid
    mat <- mtb[, samps]
    tibble(
      feature = rownames(mtb),
      prevalence = rowMeans(mat > 0),
      missingness = rowMeans(mat == 0)
    )
  })

mtb2<-fread(dt2)%>%as.data.frame()%>%
  dplyr::rename(FeatureID=row.ID)

mean_df <- md %>%
  filter(sampleid %in% colnames(mtb2)) %>%
  group_by(ATTRIBUTE_geno_n_feedingcond) %>%
  group_modify(~ {
    samps <- .x$sampleid
    mat <- mtb2[, samps]
    tibble(
      feature = rownames(mtb2),
      mean_clr    = rowMeans(mat)
    )
  })

temporal_variance <- md %>%
  filter(sampleid %in% colnames(mtb2)) %>%
  group_by(ATTRIBUTE_geno_n_feedingcond, ATTRIBUTE_ZT_collection_time) %>%
  group_modify(~ {
    samps <- .x$sampleid
    mat <- mtb2[, samps]
    tibble(
      feature = rownames(mtb2),
      mean_zt = rowMeans(mat)
    )
  }) %>%
  group_by(ATTRIBUTE_geno_n_feedingcond, feature) %>%
  summarise(temporal_variance = var(mean_zt), .groups = "drop")

#p<-ggplot(prevalence_df, aes(x = missingness, color=ATTRIBUTE_geno_n_feedingcond)) +
#p<-ggplot(prevalence_df, aes(x = prevalence, color=ATTRIBUTE_geno_n_feedingcond)) +
#p<-ggplot(mean_df, aes(x = mean_clr, color=ATTRIBUTE_geno_n_feedingcond)) +
p<-ggplot(temporal_variance, aes(x = temporal_variance, color=ATTRIBUTE_geno_n_feedingcond)) +
  scale_x_log10()+
  stat_ecdf() +
  #scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr() +
  #labs(x = "missingness",
  #labs(x = "prevalence",
  #labs(x = "mean peak abundance (log10)", 
  labs(x = "Temporal Variance (log10)", 
       y = "Cumulative Proportion",
       color = "Group")
ggsave("metacycle/SFR26_1013_tempvar_ecdf.pdf", plot=p,height=4, width=4)
#ggsave("metacycle/SFR26_1013_meanabunlog_ecdf.pdf", plot=p,height=4, width=4)
#ggsave("metacycle/SFR26_1013_missingness_ecdf.pdf", plot=p,height=4, width=4)
#ggsave("metacycle/SFR26_1013_prevalence_ecdf.pdf", plot=p,height=4, width=4)


groups <- c("NTg-ALF", "Tg-ALF", "Tg-TRF")

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      prevalence_df$prevalence[prevalence_df$ATTRIBUTE_geno_n_feedingcond == pr[1]],
      prevalence_df$prevalence[prevalence_df$ATTRIBUTE_geno_n_feedingcond == pr[2]]
    )
    data.frame(
      group1 = pr[1],
      group2 = pr[2],
      D = round(ks$statistic, 3),
      p_value = round(ks$p.value, 4)
    )
  }) %>%
  bind_rows()

# group1 group2     D p_value
# D...1 NTg-ALF Tg-ALF 0.058       0
# D...2 NTg-ALF Tg-TRF 0.121       0
# D...3  Tg-ALF Tg-TRF 0.098       0

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      prevalence_df$missingness[prevalence_df$ATTRIBUTE_geno_n_feedingcond == pr[1]],
      prevalence_df$missingness[prevalence_df$ATTRIBUTE_geno_n_feedingcond == pr[2]]
    )
    data.frame(
      group1 = pr[1],
      group2 = pr[2],
      D = round(ks$statistic, 3),
      p_value = round(ks$p.value, 4)
    )
  }) %>%
  bind_rows()

# group1 group2     D p_value
# D...1 NTg-ALF Tg-ALF 0.058       0
# D...2 NTg-ALF Tg-TRF 0.121       0
# D...3  Tg-ALF Tg-TRF 0.098       0

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      mean_df$mean_clr[mean_df$ATTRIBUTE_geno_n_feedingcond == pr[1]],
      mean_df$mean_clr[mean_df$ATTRIBUTE_geno_n_feedingcond == pr[2]]
    )
    data.frame(
      group1 = pr[1],
      group2 = pr[2],
      D = round(ks$statistic, 3),
      p_value = round(ks$p.value, 4)
    )
  }) %>%
  bind_rows()

# group1 group2     D p_value
# D...1 NTg-ALF Tg-ALF 0.044   1e-04
# D...2 NTg-ALF Tg-TRF 0.109   0e+00
# D...3  Tg-ALF Tg-TRF 0.082   0e+00

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      temporal_variance$temporal_variance[temporal_variance$ATTRIBUTE_geno_n_feedingcond == pr[1]],
      temporal_variance$temporal_variance[temporal_variance$ATTRIBUTE_geno_n_feedingcond == pr[2]]
    )
    data.frame(
      group1 = pr[1],
      group2 = pr[2],
      D = round(ks$statistic, 3),
      p_value = round(ks$p.value, 4)
    )
  }) %>%
  bind_rows()

# group1 group2     D p_value
# D...1 NTg-ALF Tg-ALF 0.026  0.0718
# D...2 NTg-ALF Tg-TRF 0.252  0.0000
# D...3  Tg-ALF Tg-TRF 0.234  0.0000
