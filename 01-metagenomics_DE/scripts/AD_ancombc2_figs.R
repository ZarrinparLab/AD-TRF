setwd("~/scratch/AD_EM/metagenomics/ancombc")

library(data.table)
library(tidyverse)

#####################################################################################

#get zt_time results

ZT_res<-function(primres){
  df_zt_time = primres %>%
    dplyr::select(taxon, ends_with("zt_time")) 
  df_fig_zt_time = df_zt_time %>%
    #dplyr::filter(diff_zt_time == 1)%>% 
    dplyr::arrange(desc(lfc_zt_time)) %>%
    dplyr::filter(p_zt_time < 0.05 & passed_ss_zt_time == 1) %>%
    #dplyr::filter(q_zt_time < 0.05 & passed_ss_zt_time == 1) %>%
    dplyr::mutate(direct = ifelse(lfc_zt_time > 0, "Positive LFC", "Negative LFC"))
  df_fig_zt_time$taxon = factor(df_fig_zt_time$taxon, levels = df_fig_zt_time$taxon)
  df_fig_zt_time$direct = factor(df_fig_zt_time$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_zt_time)
}

ZT_fig<-function(df_fig_zt_time){
  fig_zt_time = df_fig_zt_time %>%
    ggplot(aes(x = taxon, y = lfc_zt_time, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_zt_time - se_zt_time, ymax = lfc_zt_time + se_zt_time), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = "Log fold changes as one unit increase of ZT") + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_zt_time)
}

fasted_res<-function(primres){
  df_time = primres %>%
    dplyr::select(taxon, ends_with("fastedfed")) 
  df_fig_fastedfed = df_time %>%
    #dplyr::filter(diff_zt_time == 1)%>% 
    dplyr::arrange(desc(lfc_fastedfed)) %>%
    dplyr::filter(p_fastedfed < 0.05 & passed_ss_fastedfed == 1) %>%
    #dplyr::filter(q_fastedfed < 0.05 & passed_ss_fastedfed == 1) %>%
    dplyr::mutate(direct = ifelse(lfc_fastedfed > 0, "Positive LFC", "Negative LFC"))
  df_fig_fastedfed$taxon = factor(df_fig_fastedfed$taxon, levels = df_fig_fastedfed$taxon)
  df_fig_fastedfed$direct = factor(df_fig_fastedfed$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_fastedfed)
}

fastedfed_fig<-function(df_fig_fastedfed){
  fig_time = df_fig_fastedfed %>%
    ggplot(aes(x = taxon, y = lfc_fastedfed, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_fastedfed - se_fastedfed, ymax = lfc_fastedfed+ se_fastedfed), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = "Log fold changes fasted vs. fed") + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_time)
}

light_res<-function(primres){
  df_time = primres %>%
    dplyr::select(taxon, ends_with("diurnallight")) 
  df_fig_diurnallight = df_time %>%
    #dplyr::filter(diff_diurnallight == 1)%>% 
    dplyr::mutate(lfc_diurnallight=-1*lfc_diurnallight)%>%
    dplyr::arrange(desc(lfc_diurnallight)) %>%
    dplyr::filter(p_diurnallight < 0.05 & passed_ss_diurnallight == 1) %>%
    #dplyr::filter(q_diurnallight < 0.05 & passed_ss_diurnallight == 1) %>%
    dplyr::mutate(direct = ifelse(lfc_diurnallight > 0, "Positive LFC", "Negative LFC"))
  df_fig_diurnallight$taxon = factor(df_fig_diurnallight$taxon, levels = df_fig_diurnallight$taxon)
  df_fig_diurnallight$direct = factor(df_fig_diurnallight$direct, 
                                   levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_diurnallight)
}

diurnallight_fig<-function(df_fig_diurnallight){
  fig_time = df_fig_diurnallight %>%
    ggplot(aes(x = taxon, y = lfc_diurnallight, fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lfc_diurnallight - se_diurnallight, ymax = lfc_diurnallight+ se_diurnallight), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = "Log fold changes light vs. dark") + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_time)
}

cond_res<-function(primres,group){
  if(group=="ALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.ALF")) 
    df_fig_time = df_time %>%
      #dplyr::filter(diff_genotype_dietTg.ALF == 1)%>% 
      dplyr::arrange(desc(lfc_genotype_dietTg.ALF)) %>%
      dplyr::filter(p_genotype_dietTg.ALF < 0.05 & passed_ss_genotype_dietTg.ALF == 1) %>%
      dplyr::mutate(direct = ifelse(lfc_genotype_dietTg.ALF > 0, "Positive LFC", "Negative LFC"))
  }
  else if(group=="NtgALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Ntg.ALF")) 
    df_fig_time = df_time %>%
      #dplyr::filter(diff_genotype_dietNtg.ALF == 1)%>% 
      dplyr::mutate(lfc_genotype_dietNtg.ALF=-1*lfc_genotype_dietNtg.ALF)%>%
      dplyr::arrange(desc(lfc_genotype_dietNtg.ALF)) %>%
      dplyr::filter(p_genotype_dietNtg.ALF < 0.05 & passed_ss_genotype_dietNtg.ALF == 1) %>%
      dplyr::mutate(direct = ifelse(lfc_genotype_dietNtg.ALF > 0, "Positive LFC", "Negative LFC"))
  }
  else{
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.TRF")) 
    df_fig_time = df_time %>%
      #dplyr::filter(diff_genotype_dietTg.TRF == 1)%>% 
      dplyr::arrange(desc(lfc_genotype_dietTg.TRF)) %>%
      dplyr::filter(p_genotype_dietTg.TRF < 0.05 & passed_ss_genotype_dietTg.TRF == 1) %>%
      dplyr::mutate(direct = ifelse(lfc_genotype_dietTg.TRF > 0, "Positive LFC", "Negative LFC"))
  }
  df_fig_time$taxon = factor(df_fig_time$taxon, levels = df_fig_time$taxon)
  df_fig_time$direct = factor(df_fig_time$direct, 
                                  levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_time)
}

condprim_fig<-function(df_fig_time,group,title_name){
  if(group=="ALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = lfc_genotype_dietTg.ALF, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = lfc_genotype_dietTg.ALF - se_genotype_dietTg.ALF, ymax = lfc_genotype_dietTg.ALF + se_genotype_dietTg.ALF), 
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("purple4","darkorange"))
  }
  else if(group=="NtgALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = lfc_genotype_dietNtg.ALF, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = lfc_genotype_dietNtg.ALF - se_genotype_dietNtg.ALF, ymax = lfc_genotype_dietNtg.ALF + se_genotype_dietNtg.ALF), 
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("purple4","darkorange"))
    
  }
  else{
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = lfc_genotype_dietTg.TRF, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = lfc_genotype_dietTg.TRF - se_genotype_dietTg.TRF, ymax = lfc_genotype_dietTg.TRF + se_genotype_dietTg.TRF), 
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("forestgreen","purple4"))
  }
  fig_time<-fig_time+ labs(x = NULL, y = "Log fold change", 
         title = title_name) + 
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_time)
}

condzt<-function(primres,group){
  if(group=="ALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.ALF:zt_time")) 
    df_fig_time = df_time %>%
      dplyr::arrange(desc(`lfc_genotype_dietTg.ALF:zt_time`)) %>%
      dplyr::filter(`p_genotype_dietTg.ALF:zt_time` < 0.05 & `passed_ss_genotype_dietTg.ALF:zt_time` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietTg.ALF:zt_time` > 0, "Positive LFC", "Negative LFC"))
  }
  else if(group=="NtgALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Ntg.ALF:zt_time")) 
    df_fig_time = df_time %>%
      dplyr::mutate(`lfc_genotype_dietNtg.ALF:zt_time`=-1*`lfc_genotype_dietNtg.ALF:zt_time`)%>%
      dplyr::arrange(desc(`lfc_genotype_dietNtg.ALF:zt_time`)) %>%
      dplyr::filter(`p_genotype_dietNtg.ALF:zt_time` < 0.05 & `passed_ss_genotype_dietNtg.ALF:zt_time` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietNtg.ALF:zt_time` > 0, "Positive LFC", "Negative LFC"))
  }
  else{
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.TRF:zt_time")) 
    df_fig_time = df_time %>%
      dplyr::arrange(desc(`lfc_genotype_dietTg.TRF:zt_time`)) %>%
      dplyr::filter(`p_genotype_dietTg.TRF:zt_time` < 0.05 & `passed_ss_genotype_dietTg.TRF:zt_time` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietTg.TRF:zt_time` > 0, "Positive LFC", "Negative LFC"))
  }
  df_fig_time$taxon = factor(df_fig_time$taxon, levels = df_fig_time$taxon)
  df_fig_time$direct = factor(df_fig_time$direct, 
                                  levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_time)
}

fig_condzt<-function(df_fig_time,group,title_name){
  if(group=="ALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietTg.ALF:zt_time`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietTg.ALF:zt_time` - `se_genotype_dietTg.ALF:zt_time`, ymax = `lfc_genotype_dietTg.ALF:zt_time`+ `se_genotype_dietTg.ALF:zt_time`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("darkorange")) #scale_fill_manual(values=c("purple4","darkorange"))
  }
  else if(group=="NtgALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietNtg.ALF:zt_time`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietNtg.ALF:zt_time` - `se_genotype_dietNtg.ALF:zt_time`, ymax = `lfc_genotype_dietNtg.ALF:zt_time`+ `se_genotype_dietNtg.ALF:zt_time`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("darkorange")) #scale_fill_manual(values=c("purple4","darkorange"))
  }
  else{
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietTg.TRF:zt_time`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietTg.TRF:zt_time` - `se_genotype_dietTg.TRF:zt_time`, ymax = `lfc_genotype_dietTg.TRF:zt_time`+ `se_genotype_dietTg.TRF:zt_time`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("forestgreen","purple4"))
  }
  fig_time<-fig_time+ labs(x = NULL, y = "Log fold change", 
         title = title_name) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_time)
}

condfasted<-function(primres,group){
  if(group=="ALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.ALF:fastedfed")) 
    df_fig_time = df_time %>%
      dplyr::arrange(desc(`lfc_genotype_dietTg.ALF:fastedfed`)) %>%
      dplyr::filter(`p_genotype_dietTg.ALF:fastedfed` < 0.05 & `passed_ss_genotype_dietTg.ALF:fastedfed` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietTg.ALF:fastedfed` > 0, "Positive LFC", "Negative LFC"))
  }
  else if(group=="NtgALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Ntg.ALF:fastedfed")) 
    df_fig_time = df_time %>%
      dplyr::mutate(`lfc_genotype_dietNtg.ALF:fastedfed`=-1*`lfc_genotype_dietNtg.ALF:fastedfed`)%>%
      dplyr::arrange(desc(`lfc_genotype_dietNtg.ALF:fastedfed`)) %>%
      dplyr::filter(`p_genotype_dietNtg.ALF:fastedfed` < 0.05 & `passed_ss_genotype_dietNtg.ALF:fastedfed` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietNtg.ALF:fastedfed` > 0, "Positive LFC", "Negative LFC"))
  }
  else{
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.TRF:fastedfed")) 
    df_fig_time = df_time %>%
      dplyr::arrange(desc(`lfc_genotype_dietTg.TRF:fastedfed`)) %>%
      dplyr::filter(`p_genotype_dietTg.TRF:fastedfed` < 0.05 & `passed_ss_genotype_dietTg.TRF:fastedfed` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietTg.TRF:fastedfed` > 0, "Positive LFC", "Negative LFC"))
  }
  df_fig_time$taxon = factor(df_fig_time$taxon, levels = df_fig_time$taxon)
  df_fig_time$direct = factor(df_fig_time$direct, 
                              levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_time)
}

fig_condfasted<-function(df_fig_time,group,title_name){
  if(group=="ALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietTg.ALF:fastedfed`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietTg.ALF:fastedfed` - `se_genotype_dietTg.ALF:fastedfed`, ymax = `lfc_genotype_dietTg.ALF:fastedfed`+ `se_genotype_dietTg.ALF:fastedfed`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("darkorange")) #scale_fill_manual(values=c("purple4","darkorange"))
  }
  else if(group=="NtgALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietNtg.ALF:fastedfed`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietNtg.ALF:fastedfed` - `se_genotype_dietNtg.ALF:fastedfed`, ymax = `lfc_genotype_dietNtg.ALF:fastedfed`+ `se_genotype_dietNtg.ALF:fastedfed`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("purple4","darkorange"))
  }
  else{
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietTg.TRF:fastedfed`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietTg.TRF:fastedfed` - `se_genotype_dietTg.TRF:fastedfed`, ymax = `lfc_genotype_dietTg.TRF:fastedfed`+ `se_genotype_dietTg.TRF:fastedfed`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("forestgreen","purple4"))
  }
  fig_time<-fig_time+ labs(x = NULL, y = "Log fold change", 
                           title = title_name) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_time)
}

condlight<-function(primres,group){
  if(group=="ALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.ALF:diurnallight")) 
    df_fig_time = df_time %>%
      dplyr::mutate(`lfc_genotype_dietTg.ALF:diurnallight`=-1*`lfc_genotype_dietTg.ALF:diurnallight`)%>%
      dplyr::arrange(desc(`lfc_genotype_dietTg.ALF:diurnallight`)) %>%
      dplyr::filter(`p_genotype_dietTg.ALF:diurnallight` < 0.05 & `passed_ss_genotype_dietTg.ALF:diurnallight` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietTg.ALF:diurnallight` > 0, "Positive LFC", "Negative LFC"))
  }
  else if(group=="NtgALF"){
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Ntg.ALF:diurnallight")) 
    df_fig_time = df_time %>%
      dplyr::arrange(desc(`lfc_genotype_dietNtg.ALF:diurnallight`)) %>%
      dplyr::filter(`p_genotype_dietNtg.ALF:diurnallight` < 0.05 & `passed_ss_genotype_dietNtg.ALF:diurnallight` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietNtg.ALF:diurnallight` > 0, "Positive LFC", "Negative LFC"))
  }
  else{
    df_time = primres %>%
      dplyr::select(taxon, ends_with("Tg.TRF:diurnallight")) 
    df_fig_time = df_time %>%
      dplyr::mutate(`lfc_genotype_dietTg.TRF:diurnallight`=-1*`lfc_genotype_dietTg.TRF:diurnallight`)%>%
      dplyr::arrange(desc(`lfc_genotype_dietTg.TRF:diurnallight`)) %>%
      dplyr::filter(`p_genotype_dietTg.TRF:diurnallight` < 0.05 & `passed_ss_genotype_dietTg.TRF:diurnallight` == 1) %>%
      dplyr::mutate(direct = ifelse(`lfc_genotype_dietTg.TRF:diurnallight` > 0, "Positive LFC", "Negative LFC"))
  }
  df_fig_time$taxon = factor(df_fig_time$taxon, levels = df_fig_time$taxon)
  df_fig_time$direct = factor(df_fig_time$direct, 
                              levels = c("Positive LFC", "Negative LFC"))
  return(df_fig_time)
}

fig_condlight<-function(df_fig_time,group,title_name){
  if(group=="ALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietTg.ALF:diurnallight`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietTg.ALF:diurnallight` - `se_genotype_dietTg.ALF:diurnallight`, ymax = `lfc_genotype_dietTg.ALF:diurnallight`+ `se_genotype_dietTg.ALF:diurnallight`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("darkorange")) #scale_fill_manual(values=c("purple4","darkorange"))
  }
  else if(group=="NtgALF"){
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietNtg.ALF:diurnallight`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietNtg.ALF:diurnallight` - `se_genotype_dietNtg.ALF:diurnallight`, ymax = `lfc_genotype_dietNtg.ALF:diurnallight`+ `se_genotype_dietNtg.ALF:diurnallight`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("darkorange")) #scale_fill_manual(values=c("purple4","darkorange"))
  }
  else{
    fig_time = df_fig_time %>%
      ggplot(aes(x = taxon, y = `lfc_genotype_dietTg.TRF:diurnallight`, fill = direct)) + 
      geom_bar(stat = "identity", width = 0.7, color = "black", 
               position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = `lfc_genotype_dietTg.TRF:diurnallight` - `se_genotype_dietTg.TRF:diurnallight`, ymax = `lfc_genotype_dietTg.TRF:diurnallight`+ `se_genotype_dietTg.TRF:diurnallight`),
                    width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("forestgreen","purple4"))
  }
  fig_time<-fig_time+ labs(x = NULL, y = "Log fold change", 
                           title = title_name) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
  return(fig_time)
}

global_res<-function(primres,globalres,group){
  df_geno = primres %>%
    dplyr::select(taxon, contains("genotype_diet")) 
  df_fig_global = df_geno %>%
    dplyr::left_join(globalres %>%
                       dplyr::transmute(taxon, 
                                        p_val=p_val,
                                        diff_geno = diff_abn, 
                                        passed_ss = passed_ss))%>%
    dplyr::filter(p_val<0.05 & passed_ss==1) #%>% dplyr::filter(diff_geno == 1) %>%
  if(group=="NtgALF"){
    df_fig_global = df_fig_global%>%
      dplyr::mutate(lfc_ALF = lfc_genotype_dietTg.ALF,
                    lfc_TRF = lfc_genotype_dietTg.TRF) %>%
      dplyr::transmute(taxon,
                       `TgALF - NTgALF` = round(lfc_ALF, 2),
                       `TgTRF - NTgALF` = round(lfc_TRF, 2)) %>%
      tidyr::pivot_longer(cols = `TgALF - NTgALF`:`TgTRF - NTgALF`, 
                          names_to = "group", values_to = "value") %>%
      dplyr::arrange(taxon)
    
    df_fig_global$group = factor(df_fig_global$group, 
                                 levels = c("TgALF - NTgALF",
                                            "TgTRF - NTgALF"))
  }
  else{
    df_fig_global = df_fig_global%>%
      dplyr::mutate(lfc_ALF = -1*lfc_genotype_dietNtg.ALF,
                                  lfc_TRF = lfc_genotype_dietTg.TRF) %>%
      dplyr::transmute(taxon,
                       `TgALF - NTgALF` = round(lfc_ALF, 2),
                       `TgTRF - TgALF` = round(lfc_TRF, 2)) %>%
      tidyr::pivot_longer(cols = `TgALF - NTgALF`:`TgTRF - TgALF`, 
                          names_to = "group", values_to = "value") %>%
      dplyr::arrange(taxon)
    
    df_fig_global$group = factor(df_fig_global$group, 
                                 levels = c("TgALF - NTgALF",
                                            "TgTRF - TgALF"))
  }

  
  return(df_fig_global)
}

global_fig<-function(df_fig_global){
  lo = floor(min(df_fig_global$value))
  up = ceiling(max(df_fig_global$value))
  mid = (lo + up)/2
  fig_global = df_fig_global %>%
    ggplot(aes(x = group, y = taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = NULL) +
    geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
    labs(x = NULL, y = NULL, title = "Log fold changes for globally significant taxa") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  return(fig_global)
}

dunnett_res<-function(dunnres,group){
  if(group=="NtgALF"){
    df_fig_dunn1 = dunnres %>%
      # dplyr::filter(diff_genotype_dietTg.ALF == 1 | 
      #                 diff_genotype_dietTg.TRF == 1) %>%
      dplyr::filter(p_genotype_dietTg.ALF < 0.05 | 
                      p_genotype_dietTg.TRF < 0.05) %>%
      dplyr::mutate(lfc1 = ifelse(p_genotype_dietTg.ALF < 0.05, 
                                  round(lfc_genotype_dietTg.ALF, 2), 0),
                    lfc2 = ifelse(p_genotype_dietTg.TRF < 0.05, 
                                  round(lfc_genotype_dietTg.TRF, 2), 0)) %>%
      tidyr::pivot_longer(cols = lfc1:lfc2, 
                          names_to = "group", values_to = "value") %>%
      dplyr::arrange(taxon)
    
    df_fig_dunn2 = dunnres %>%
      dplyr::filter(p_genotype_dietTg.ALF < 0.05 | 
                      p_genotype_dietTg.TRF < 0.05) %>%
      dplyr::mutate(lfc1 = ifelse(passed_ss_genotype_dietTg.ALF == 1 & p_genotype_dietTg.ALF < 0.05, 
                                  "aquamarine3", "black"),
                    lfc2 = ifelse(passed_ss_genotype_dietTg.TRF == 1 & p_genotype_dietTg.TRF < 0.05, 
                                  "aquamarine3", "black")) %>%
      tidyr::pivot_longer(cols = lfc1:lfc2, 
                          names_to = "group", values_to = "color") %>%
      dplyr::arrange(taxon)
    
    df_fig_dunn = df_fig_dunn1 %>%
      dplyr::left_join(df_fig_dunn2, by = c("taxon", "group"))
    
    df_fig_dunn$group = recode(df_fig_dunn$group, 
                               `lfc1` = "NTg-ALF - Tg-ALF",
                               `lfc2` = "NTg-ALF - Tg-TRF")
    df_fig_dunn$group = factor(df_fig_dunn$group, 
                               levels = c("NTg-ALF - Tg-ALF",
                                          "NTg-ALF - Tg-TRF"))
  }
  else{
    df_fig_dunn1 = dunnres %>%
      # dplyr::filter(diff_genotype_dietNtg.ALF == 1 | 
      #                 diff_genotype_dietTg.TRF == 1) %>%
      dplyr::filter(p_genotype_dietNtg.ALF < 0.05 | 
                      p_genotype_dietTg.TRF < 0.05) %>%
      dplyr::mutate(lfc_genotype_dietNtg.ALF=-1*lfc_genotype_dietNtg.ALF)%>%
      dplyr::mutate(lfc1 = ifelse(p_genotype_dietNtg.ALF < 0.05, 
                                  round(lfc_genotype_dietNtg.ALF, 2), 0),
                    lfc2 = ifelse(p_genotype_dietTg.TRF < 0.05, 
                                  round(lfc_genotype_dietTg.TRF, 2), 0)) %>%
      tidyr::pivot_longer(cols = lfc1:lfc2, 
                          names_to = "group", values_to = "value") %>%
      dplyr::arrange(taxon)
    
    df_fig_dunn2 = dunnres %>%
      dplyr::filter(p_genotype_dietNtg.ALF < 0.05 | 
                      p_genotype_dietTg.TRF < 0.05) %>%
      dplyr::mutate(lfc1 = ifelse(passed_ss_genotype_dietNtg.ALF == 1 & p_genotype_dietNtg.ALF < 0.05, 
                                  "aquamarine3", "black"),
                    lfc2 = ifelse(passed_ss_genotype_dietTg.TRF == 1 & p_genotype_dietTg.TRF < 0.05, 
                                  "aquamarine3", "black")) %>%
      tidyr::pivot_longer(cols = lfc1:lfc2, 
                          names_to = "group", values_to = "color") %>%
      dplyr::arrange(taxon)
    
    df_fig_dunn = df_fig_dunn1 %>%
      dplyr::left_join(df_fig_dunn2, by = c("taxon", "group"))
    
    df_fig_dunn$group = recode(df_fig_dunn$group, 
                               `lfc1` = "NTg-ALF - Tg-ALF",
                               `lfc2` = "Tg-ALF - Tg-TRF")
    df_fig_dunn$group = factor(df_fig_dunn$group, 
                               levels = c("NTg-ALF - Tg-ALF",
                                          "Tg-ALF - Tg-TRF"))
  }
  return(df_fig_dunn)
}

dunett_fig<-function(df_fig_dunn,title_name){
  lo = floor(min(df_fig_dunn$value))
  up = ceiling(max(df_fig_dunn$value))
  mid = (lo + up)/2
  fig_dunn = df_fig_dunn %>%
    ggplot(aes(x = group, y = taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = NULL) +
    geom_text(aes(group, taxon, label = value, color = color), size = 4) +
    scale_color_identity(guide = FALSE) +
    labs(x = NULL, y = NULL, title = title_name) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  return(fig_dunn)
}

#####################################################################################

#ZT results
res_prim<-fread("ancombbc_primaryres_B10.txt")
df_fig_zt_time<-ZT_res(res_prim)
ZT_fig(df_fig_zt_time)
ggsave("lfc_zt_time_ancombc_B10.pdf", height=10, width=30)

res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
df_fig_zt_time<-ZT_res(res_prim)
ZT_fig(df_fig_zt_time)
ggsave("lfc_zt_timetrf_ancombc_B10.pdf", height=11, width=15)

#fasted/fed results
res_prim<-fread("ancombbc_primaryresfasted_B10.txt")
df_fig_fastedfed<-fasted_res(res_prim)
fastedfed_fig(df_fig_fastedfed)
ggsave("lfc_fastedfed_ancombc_B10.pdf", height=11, width=30)

res_prim<-fread("ancombbc_primaryresfastedtrf_B10.txt")
df_fig_fastedfed<-fasted_res(res_prim)
fastedfed_fig(df_fig_fastedfed)
ggsave("lfc_fastedfedtrf_ancombc_B10.pdf", height=11, width=6)

#light/dark results
res_prim<-fread("ancombbc_primaryreslighttrf_B10.txt")
df_fig_lightdark<-light_res(res_prim)
diurnallight_fig(df_fig_lightdark)
ggsave("lfc_lightdarktrf_ancombc_B10.pdf", height=11, width=20)

#results for NtgALFtgALF
res_prim<-fread("ancombbc_primaryres_B10.txt")
df_fig_alf_time<-cond_res(res_prim,"ALF")
condprim_fig(df_fig_alf_time,"ALF","Log fold changes NTg-ALF vs. Tg-ALF")
ggsave("lfc_zt_tgALFntgALF_ancombc_B10.pdf", height=10, width=5)

res_prim<-fread("ancombbc_primaryresfasted_B10.txt")
df_fig_alf_time<-cond_res(res_prim,"ALF")
condprim_fig(df_fig_alf_time,"ALF","Log fold changes NTg-ALF vs. Tg-ALF")
ggsave("lfc_zt_tgALFntgALFfasted_ancombc_B10.pdf", height=11, width=7)

res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
df_fig_alf_time<-cond_res(res_prim,"NtgALF")
condprim_fig(df_fig_alf_time,"NtgALF","Log fold changes NTg-ALF vs. Tg-ALF")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFntgALFtrf_ancombc_B10_cln.pdf", height=4.5, width=10)

res_prim<-fread("ancombbc_primaryresfastedtrf_B10.txt")
df_fig_alf_time<-cond_res(res_prim,"NtgALF")
condprim_fig(df_fig_alf_time,"NtgALF","Log fold changes NTg-ALF vs. Tg-ALF")
ggsave("lfc_zt_tgALFntgALFfastedtrf_ancombc_B10.pdf", height=11, width=9)

res_prim<-fread("ancombbc_primaryreslighttrf_B10.txt")
df_fig_alf_time<-cond_res(res_prim,"NtgALF")
condprim_fig(df_fig_alf_time,"NtgALF","Log fold changes NTg-ALF vs. Tg-ALF")
ggsave("lfc_zt_tgALFntgALFlighttrf_ancombc_B10.pdf", height=10.5, width=5)

#results for NtgALFtgALF w/ interaction
res_prim<-fread("ancombbc_primaryres_B10.txt")
df_fig_alf_time<-condzt(res_prim,"ALF")
fig_condzt(df_fig_alf_time,"ALF","NTg-ALF vs. Tg-ALF w/ time interaction")
ggsave("lfc_zt_tgALFntgALFwzt_ancombc_B10.pdf", height=9, width=3.5)

res_prim<-fread("ancombbc_primaryresfasted_B10.txt")
df_fig_alf_time<-condfasted(res_prim,"ALF")
fig_condfasted(df_fig_alf_time,"ALF","NTg-ALF vs. Tg-ALF w/ time interaction")
ggsave("lfc_zt_tgALFntgALFwztfasted_ancombc_B10.pdf", height=10, width=20)

res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
df_fig_alf_time<-condzt(res_prim,"NtgALF")
fig_condzt(df_fig_alf_time,"NtgALF","NTg-ALF vs. Tg-ALF w/ time interaction")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFntgALFwzttrf_ancombc_B10_cln.pdf", height=3, width=9.5)

res_prim<-fread("ancombbc_primaryresfastedtrf_B10.txt")
df_fig_alf_time<-condfasted(res_prim,"NtgALF")%>%
  filter(`lfc_genotype_dietNtg.ALF:fastedfed`>0.75|`lfc_genotype_dietNtg.ALF:fastedfed`< -0.75)
fig_condfasted(df_fig_alf_time,"NtgALF","NTg-ALF vs. Tg-ALF w/ time interaction")
ggsave("lfc_zt_tgALFntgALFwztfastedtrf_ancombc_B10.pdf", height=12, width=20)
ggsave("lfc_zt_tgALFntgALFwztfastedtrffold0.75_ancombc_B10.pdf", height=11, width=7)

df_fig_alf_time<-condfasted(res_prim,"NtgALF")%>%
  filter(str_detect(taxon,"Prevotella"))
fig_condfasted(df_fig_alf_time,"NtgALF","NTg-ALF vs. Tg-ALF w/ time interaction")
ggsave("lfc_zt_tgALFntgALFwztfastedtrfprevotella_ancombc_B10.pdf", height=9, width=7)

res_prim<-fread("ancombbc_primaryreslighttrf_B10.txt")
df_fig_alf_time<-condlight(res_prim,"NtgALF")
fig_condlight(df_fig_alf_time,"NtgALF","NTg-ALF vs. Tg-ALF w/ time interaction")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFntgALFwztlighttrf_ancombc_B10_cln.pdf", height=2, width=8.5)

res_prim<-fread("ancombbc_primaryreslighttrf2_B10.txt")
df_fig_alf_time<-condlight(res_prim,"NtgALF")
fig_condlight(df_fig_alf_time,"NtgALF","NTg-ALF vs. Tg-ALF w/ time interaction")+
  scale_fill_manual(values=c("purple4"))+ coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFntgALFwztlighttrf2_ancombc_B10_cln.pdf", height=2, width=8)

#results for tgALFtgTRF
res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
df_fig_trf_time<-cond_res(res_prim,"TRF")
condprim_fig(df_fig_trf_time,"TRF","Log fold changes Tg-ALF vs. Tg-TRF")
ggsave("lfc_zt_tgALFtgTRF_ancombc_B10.pdf", height=10, width=20)

res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
df_fig_trf_time<-cond_res(res_prim,"TRF")%>%
  filter(lfc_genotype_dietTg.TRF>2.25|lfc_genotype_dietTg.TRF<0)
condprim_fig(df_fig_trf_time,"TRF","Log fold changes Tg-ALF vs. Tg-TRF")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFtgTRFfold2_ancombc_B10_cln.pdf", height=4.5, width=10.5)

res_prim<-fread("ancombbc_primaryresfastedtrf_B10.txt")
df_fig_trf_time<-cond_res(res_prim,"TRF")%>%
  filter(lfc_genotype_dietTg.TRF>1.5|lfc_genotype_dietTg.TRF<0)
condprim_fig(df_fig_trf_time,"TRF","Log fold changes Tg-ALF vs. Tg-TRF")
ggsave("lfc_zt_tgALFtgTRFfasted_ancombc_B10.pdf", height=11, width=20)
ggsave("lfc_zt_tgALFtgTRFfastedfold1.5_ancombc_B10.pdf", height=11, width=10)

res_prim<-fread("ancombbc_primaryresfastedtrf_B10.txt")
df_fig_trf_time<-cond_res(res_prim,"TRF")%>%
  #filter(str_detect(taxon,"Limosilactobacillus"))
  filter(str_detect(taxon,"Ligilactobacillus"))
condprim_fig(df_fig_trf_time,"TRF","Log fold changes Tg-ALF vs. Tg-TRF")
ggsave("lfc_zt_tgALFtgTRFfastedligi_ancombc_B10.pdf", height=9, width=5)

res_prim<-fread("ancombbc_primaryreslighttrf_B10.txt")
df_fig_trf_time<-cond_res(res_prim,"TRF")
condprim_fig(df_fig_trf_time,"TRF","Log fold changes Tg-ALF vs. Tg-TRF")
ggsave("lfc_zt_tgALFtgTRFlight_ancombc_B10.pdf", height=11, width=20)

#results for tgALFtgTRF w/ interaction
res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
df_fig_trf_time<-condzt(res_prim,"TRF")
fig_condzt(df_fig_trf_time,"TRF","Tg-ALF vs. Tg-TRF w/ time interaction")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFtgTRFwzt_ancombc_B10_cln.pdf", height=4, width=10)

res_prim<-fread("ancombbc_primaryresfastedtrf_B10.txt")
df_fig_trf_time<-condfasted(res_prim,"TRF") %>%
  filter(`lfc_genotype_dietTg.TRF:fastedfed`>0.75|`lfc_genotype_dietTg.TRF:fastedfed`< -0.75)
  #filter(str_detect(taxon,"Limosilactobacillus"))
  #filter(str_detect(taxon,"Ligilactobacillus"))
fig_condfasted(df_fig_trf_time,"TRF","Tg-ALF vs. Tg-TRF w/ time interaction")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
#ggsave("lfc_zt_tgALFtgTRFwfasted_ancombc_B10.pdf", height=11, width=20)
ggsave("lfc_zt_tgALFtgTRFwfasted_ancombc_B10_cln.pdf", height=7, width=10)

res_prim<-fread("ancombbc_primaryreslighttrf_B10.txt")
df_fig_trf_time<-condlight(res_prim,"TRF")
fig_condlight(df_fig_trf_time,"TRF","Tg-ALF vs. Tg-TRF w/ time interaction")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("lfc_zt_tgALFtgTRFwlight_ancombc_B10_cln.pdf", height=3.5, width=9.5)

res_prim<-fread("ancombbc_primaryreslighttrf2_B10.txt")#226
df_fig_trf_time<-condlight(res_prim,"TRF")%>%
  filter(`lfc_genotype_dietTg.TRF:diurnallight`>1|`lfc_genotype_dietTg.TRF:diurnallight`< -1)
fig_condlight(df_fig_trf_time,"TRF","Tg-ALF vs. Tg-TRF w/ time interaction")+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
#ggsave("lfc_zt_tgALFtgTRFwlight2_ancombc_B10.pdf", height=10.5, width=20)
ggsave("lfc_zt_tgALFtgTRFwlight2_ancombc_B10_cln.pdf", height=6, width=11)

#global test
res_prim<-fread("ancombbc_primaryres_B10.txt")
res_global<-fread("ancombbc_globalres_B10.txt")
df_fig_global<-global_res(res_prim,res_global,"NtgALF")
global_fig(df_fig_global)
ggsave("lfc_globalsig_ancombc_B10.pdf", height=15, width=15)

res_prim<-fread("ancombbc_primaryresfasted_B10.txt")
res_global<-fread("ancombbc_globalresfasted_B10.txt")
df_fig_global<-global_res(res_prim,res_global,"NtgALF")
global_fig(df_fig_global)
ggsave("lfc_globalsigfasted_ancombc_B10.pdf", height=20, width=15)

res_prim<-fread("ancombbc_primaryrestrf_B10.txt")
res_global<-fread("ancombbc_globalrestrf_B10.txt")
df_fig_global<-global_res(res_prim,res_global,"tgALF")
global_fig(df_fig_global)
ggsave("lfc_globalsigtrf_ancombc_B10.pdf", height=15, width=15)

res_prim<-fread("ancombbc_primaryreslighttrf_B10.txt")
res_global<-fread("ancombbc_globalreslighttrf_B10.txt")
df_fig_global<-global_res(res_prim,res_global,"tgALF")
global_fig(df_fig_global)
ggsave("lfc_globalsiglighttrf_ancombc_B10.pdf", height=20, width=15)

#dunnett pairwise
res_dunn<-fread("ancombbc_dunnett_B10.txt")
df_fig_dunn<-dunnett_res(res_dunn,"NtgALF")
dunett_fig(df_fig_dunn,"Log fold changes as compared to NTg-ALF")
ggsave("lfc_zt_tgALFntgALF_ancombc_dunnett_B10.pdf", height=4, width=11)

res_dunn<-fread("ancombbc_dunnettfasted_B10.txt")
df_fig_dunn<-dunnett_res(res_dunn,"NtgALF")
dunett_fig(df_fig_dunn,"Log fold changes as compared to NTg-ALF")
ggsave("lfc_zt_tgALFntgALFfasted_ancombc_dunnett_B10.pdf", height=6, width=12)

res_dunn<-fread("ancombbc_dunnetttrf_B10.txt")
df_fig_dunn<-dunnett_res(res_dunn,"tgALF")
dunett_fig(df_fig_dunn,"Log fold changes as compared to Tg-ALF")
ggsave("lfc_zt_tgALFtgTRF_ancombc_dunnett_B10.pdf", height=15, width=12)

res_dunn<-fread("ancombbc_dunnettlighttrf_B10.txt")
df_fig_dunn<-dunnett_res(res_dunn,"tgALF")
dunett_fig(df_fig_dunn,"Log fold changes as compared to Tg-ALF")
ggsave("lfc_zt_tgALFtgTRFlighttrf_ancombc_dunnett_B10.pdf", height=8, width=12)

# pairwise comparison
res_pair<-fread("ancombbc_pairwisecomparres_B10.txt")
res_trend<-fread("ancombbc_pattern_B10.txt")
res_pair<-fread("ancombbc_pairwisecomparrestrf_B10.txt")
res_trend<-fread("ancombbc_patterntrf_B10.txt")
res_pair<-fread("ancombbc_pairwisecomparresfasted_B10.txt")
res_trend<-fread("ancombbc_patternfasted_B10.txt")
res_pair<-fread("ancombbc_pairwisecomparreslighttrf_B10.txt")
res_trend<-fread("ancombbc_patternlighttrf_B10.txt")

taxa<-fread("/mnt/zarrinpar/Pynchon/Databases/wol2/taxonomy/lineages.split.tsv")%>%
  dplyr::rename(taxon=genome)
prim_human<-fread("~/scratch/AD_EM/metagenomics/human_BEAM/ancombbc_primaryreshuman_B10.txt")%>%
  filter(p_memdxmci<0.05)%>%
  left_join(.,taxa,by="taxon")%>%
  mutate(label=paste(taxon,Species, sep=" "))%>%
  dplyr::mutate(direct = ifelse(lfc_memdxmci > 0, "Positive LFC", "Negative LFC"))%>%
  arrange(desc(lfc_memdxmci))%>%
  filter(lfc_memdxmci > 1.75 | lfc_memdxmci < -1.75)

prim_human$label = factor(prim_human$label, levels = prim_human$label)
prim_human$direct = factor(prim_human$direct, 
                            levels = c("Positive LFC", "Negative LFC"))


fig_time = prim_human %>%
  ggplot(aes(x = label, y = lfc_memdxmci, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_memdxmci - se_memdxmci, ymax = lfc_memdxmci + se_memdxmci), 
                width = 0.2, position = position_dodge(0.05), color = "black") + scale_fill_manual(values=c("purple4","darkorange"))+
  labs(x = NULL, y = "Log fold change") +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1))+
  coord_flip()+theme(legend.position = "top", axis.text.x = element_text(angle = 0, hjust = 1))
ggsave("~/scratch/AD_EM/metagenomics/human_BEAM/ancombbc_human1.75fold_ncmci.pdf", height=6, width=5)



# 
# df_fig_pair1 = res_pair %>%
#   dplyr::filter(diff_genotype_dietTg.ALF == 1 |
#                 diff_genotype_dietTg.TRF == 1 |
#                 diff_genotype_dietTg.TRF_genotype_dietTg.ALF == 1) %>%
#   dplyr::mutate(lfc1 = ifelse(diff_genotype_dietTg.ALF == 1,
#                               round(lfc_genotype_dietTg.ALF, 2), 0),
#                 lfc2 = ifelse(diff_genotype_dietTg.TRF == 1,
#                               round(lfc_genotype_dietTg.TRF, 2), 0),
#                 lfc3 = ifelse(diff_genotype_dietTg.TRF_genotype_dietTg.ALF == 1,
#                               round(lfc_genotype_dietTg.TRF_genotype_dietTg.ALF, 2), 0)) %>%
#   tidyr::pivot_longer(cols = lfc1:lfc3,
#                       names_to = "group", values_to = "value") %>%
#   dplyr::arrange(taxon)
# 
# df_fig_pair2 = res_pair %>%
#   dplyr::filter(diff_bmioverweight == 1 |
#                   diff_bmilean == 1 |
#                   diff_bmilean_bmioverweight == 1) %>%
#   dplyr::mutate(lfc1 = ifelse(passed_ss_bmioverweight == 1 & diff_bmioverweight == 1,
#                               "aquamarine3", "black"),
#                 lfc2 = ifelse(passed_ss_bmilean == 1 & diff_bmilean == 1,
#                               "aquamarine3", "black"),
#                 lfc3 = ifelse(passed_ss_bmilean_bmioverweight == 1 & diff_bmilean_bmioverweight == 1,
#                               "aquamarine3", "black")) %>%
#   tidyr::pivot_longer(cols = lfc1:lfc3,
#                       names_to = "group", values_to = "color") %>%
#   dplyr::arrange(taxon)
# 
# df_fig_pair = df_fig_pair1 %>%
#   dplyr::left_join(df_fig_pair2, by = c("taxon", "group"))
# 
# df_fig_pair$group = recode(df_fig_pair$group,
#                            `lfc1` = "Overweight - Obese",
#                            `lfc2` = "Lean - Obese",
#                            `lfc3` = "Lean - Overweight")
# df_fig_pair$group = factor(df_fig_pair$group,
#                            levels = c("Overweight - Obese",
#                                       "Lean - Obese",
#                                       "Lean - Overweight"))
# 
# lo = floor(min(df_fig_pair$value))
# up = ceiling(max(df_fig_pair$value))
# mid = (lo + up)/2
# fig_pair = df_fig_pair %>%
#   ggplot(aes(x = group, y = taxon, fill = value)) +
#   geom_tile(color = "black") +
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
#                        na.value = "white", midpoint = mid, limit = c(lo, up),
#                        name = NULL) +
#   geom_text(aes(group, taxon, label = value, color = color), size = 4) +
#   scale_color_identity(guide = FALSE) +
#   labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# fig_pair
