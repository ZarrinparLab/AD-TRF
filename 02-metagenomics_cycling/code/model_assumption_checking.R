setwd("~/scratch/AD_EM/metagenomics/")

library(data.table)
library(tidyverse)
library(ggpubr)

#############################################################
#inputs

dt<-"data_files/159013_wolR2_per_genome.feature-table.species.tsv"
dt2<-"data_files/159013_wolR2_per_genome.feature-table.species.clr.csv"
metat<-"data_files/stool.metagenomic.metadata.with_amyloid_saczt.tsv"
#############################################################
df_species <- fread(dt) %>%
  dplyr::rename(taxonomy = `OTU ID`) %>%
  column_to_rownames("taxonomy") %>%
  rename_with(~ gsub("^14708\\.", "", .x))

md <- fread(metat)%>%
  mutate(fasted=ifelse(zt_time<20 & zt_time>4, "fasted","fed"),
         sample_name =sub("^[^.]*\\.", "", sample_name),
         mouseid = sub("\\..*", "",sample_name))

prevalence_df <- md %>%
  filter(sample_name %in% colnames(df_species)) %>%
  group_by(genotype_diet) %>%
  group_modify(~ {
    samps <- .x$sample_name
    mat <- df_species[, samps]
    tibble(
      feature = rownames(df_species),
      prevalence = rowMeans(mat > 0),
      missingness = rowMeans(mat == 0)
    )
  })

df_species2 <- fread(dt2) %>%
  dplyr::rename(taxonomy = `OTU ID`) %>%
  column_to_rownames("taxonomy")

mean_df <- md %>%
  filter(sample_name %in% colnames(df_species2)) %>%
  group_by(genotype_diet) %>%
  group_modify(~ {
    samps <- .x$sample_name
    mat <- df_species2[, samps]
    tibble(
      feature = rownames(df_species2),
      mean_clr    = rowMeans(mat)
    )
  })

temporal_variance <- md %>%
  filter(sample_name %in% colnames(df_species2)) %>%
  group_by(genotype_diet, zt_collection_time) %>%
  group_modify(~ {
    samps <- .x$sample_name
    mat <- df_species2[, samps]
    tibble(
      feature = rownames(df_species2),
      mean_zt = rowMeans(mat)
    )
  }) %>%
  group_by(genotype_diet, feature) %>%
  summarise(temporal_variance = var(mean_zt), .groups = "drop")

#p<-ggplot(prevalence_df, aes(x = missingness, color=genotype_diet)) +
p<-ggplot(prevalence_df, aes(x = prevalence, color=genotype_diet)) +
#p<-ggplot(mean_df, aes(x = mean_clr, color=genotype_diet)) +
#p<-ggplot(temporal_variance, aes(x = temporal_variance, color = genotype_diet)) +
  stat_ecdf() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values=c("darkorange","purple4","forestgreen"))+
  theme_pubr() +
  #labs(x = "missingness",
  labs(x = "prevalence",
  #labs(x = "mean abundance (CLR)", 
  #labs(x = "Temporal Variance (CLR, log10)", 
       y = "Cumulative Proportion",
       color = "Group")
#ggsave("cycling_metacycle/SFR26_1013_tempvar_ecdf.pdf", plot=p,height=4, width=4)
#ggsave("cycling_metacycle/SFR26_1013_meanCLRabun_ecdf.pdf", plot=p,height=4, width=4)
#ggsave("cycling_metacycle/SFR26_1013_missingness_ecdf.pdf", plot=p,height=4, width=4)
ggsave("cycling_metacycle/SFR26_1013_prevalence_ecdf.pdf", plot=p,height=4, width=4)

groups <- c("Ntg.ALF", "Tg.ALF", "Tg.TRF")

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      prevalence_df$prevalence[prevalence_df$genotype_diet == pr[1]],
      prevalence_df$prevalence[prevalence_df$genotype_diet == pr[2]]
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
# D...1 Ntg.ALF Tg.ALF 0.082       0
# D...2 Ntg.ALF Tg.TRF 0.064       0
# D...3  Tg.ALF Tg.TRF 0.112       0

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      prevalence_df$missingness[prevalence_df$genotype_diet  == pr[1]],
      prevalence_df$missingness[prevalence_df$genotype_diet  == pr[2]]
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
# D...1 Ntg.ALF Tg.ALF 0.082       0
# D...2 Ntg.ALF Tg.TRF 0.064       0
# D...3  Tg.ALF Tg.TRF 0.112       0

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      mean_df$mean_clr[mean_df$genotype_diet == pr[1]],
      mean_df$mean_clr[mean_df$genotype_diet == pr[2]]
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
# D...1 Ntg.ALF Tg.ALF 0.465       0
# D...2 Ntg.ALF Tg.TRF 0.557       0
# D...3  Tg.ALF Tg.TRF 0.559       0

combn(groups, 2, simplify = FALSE) %>%
  lapply(function(pr) {
    ks <- ks.test(
      temporal_variance$temporal_variance[temporal_variance$genotype_diet == pr[1]],
      temporal_variance$temporal_variance[temporal_variance$genotype_diet == pr[2]]
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
# D...1 Ntg.ALF Tg.ALF 0.431       0
# D...2 Ntg.ALF Tg.TRF 0.464       0
# D...3  Tg.ALF Tg.TRF 0.437       0

