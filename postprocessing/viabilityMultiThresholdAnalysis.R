library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# library(PharmacoGx)
# library(caret)
source("viabilityModeling.R")


# Initialize
set <- 'NCI60_2021'
no_models <- 50


# Load data---------

# Load cmap data
drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]
data_cmap <- readRDS("../preprocessing/preprocessed_data/all_cmap_sigs_with_pert_info.rds")
TF_activities_2 <- read.delim('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv',row.names = 1)

data_cmap <- data_cmap %>% filter(sig_id %in% rownames(TF_activities_2))
data_cmap <- data_cmap %>% filter(canonical_smiles %in% rownames(drug_targets) | pert_iname %in% rownames(drug_targets))
gc()

# Load lethality/viability----------------------
data <- readRDS(paste0('./',set,'.rds'))

# Create initial linear models using different number of models----------------
# Do it for multiple thresholds to check. And then choose one
initial_models <- c('elasticnet')
performance_df <- data.frame()
for (md in initial_models){
  test_corr <- NULL
  test_pval <- NULL
  for (th in 1:no_models){
    drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
    merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
    consensus_info <- readRDS('../results/consensus_interaction_info.rds')
    colnames(consensus_info)[2] <- 'variable'
    merged_interactions_all <- left_join(consensus_info,merged_interactions_all)
    merged_interactions_all <- merged_interactions_all %>% filter(!is.na(`Prior knowledge`))
    merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = 1*(Inferred=='Interaction'))
    merged_interactions_all <- merged_interactions_all %>% mutate(`Prior knowledge` = 1*(`Prior knowledge`=='Interaction'))
    merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred",
                                                                  "model_no","consensus_inferrence","mean_frequency")
    merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = ifelse((consensus_inferrence==3) & (mean_frequency>=th/no_models),1,0))
    merged_interactions_all <- merged_interactions_all %>% select(-mean_frequency,-consensus_inferrence)
    merged_interactions_all <- distinct(merged_interactions_all)
    
    inferred_interactions <- merged_interactions_all %>% select(drug,variable,Inferred)
    inferred_interactions <- distinct(inferred_interactions)
    inferred_interactions <- inferred_interactions %>% group_by(drug) %>% mutate(interactions = sum(Inferred)) %>% ungroup()
    inferred_interactions <- inferred_interactions %>% filter(interactions>0)
    inferred_interactions <- inferred_interactions %>% select(-interactions)
    inferred_interactions <- inferred_interactions %>% spread('variable','Inferred')
    gc()
    drug_targets_inferred <- as.matrix(distinct(inferred_interactions %>% column_to_rownames('drug')))
    drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions_all$drug),
                                 which(colnames(drug_targets) %in% merged_interactions_all$variable)]
    
    data_new <- data
    data_new <- data_new %>% filter(smiles %in% merged_interactions_all$drug)
    data_new[,3:(ncol(data_new)-9)] <- apply(data_new[,3:(ncol(data_new)-9)],c(1,2),as.numeric)
    changed <- 0
    for (i in 1:nrow(data_new)){
      drug <- data_new$smiles[i]
      if (drug %in% rownames(drug_targets_inferred)){
        for (j in 3:(ncol(data_new)-9)){
          target <-   colnames(data_new)[j]
          val <- drug_targets_inferred[drug,target]
          if (val!=data_new[i,j]){
            changed <- changed+1
          }
          data_new[i,j] <- val
        }
      }
    }
    ### Call function ###
    total_results <- viability_model(data_new,
                                     colnames(drug_targets),
                                     dataset=set,
                                     lethality_data_path = "./",
                                     model = md,
                                     no_models=50)
    
    res_test_all <- total_results$res_test_all
    res_test_all <- res_test_all %>% group_by(data) %>% mutate(correlation = cor(predicted,EC50)) %>% ungroup()
    res_test_all <- res_test_all %>% mutate(frequency=th/no_models)
    res_test_all <- res_test_all %>% mutate(model=md)
    performance_df <- rbind(performance_df,
                            res_test_all)
    saveRDS(performance_df,'../results/viability_analysis_results/viability_results_LOOCV_multi_threshold_analysis.rds')
    print(paste0('Using ',md,' finished: ',th,'/',no_models))
  }
}
saveRDS(performance_df,'../results/viability_analysis_results/viability_results_LOOCV_multi_threshold_analysis.rds')

