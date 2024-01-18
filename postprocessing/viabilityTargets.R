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
md <- 'lasso'
type <- 'prior'
no_models <- 50
th <- 44


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
# # PRISM has sensitivity data
# available_datasets <- availablePSets()
# # Load data through PharmacoDB
# dataset <- PharmacoGx::downloadPSet('NCI60_2021')
# ## Plot Drug Dose response curves, using the same names for compounds and cell lines as PharmacoDB
# drugDoseResponseCurve(dataset, drug="Lestaurtinib", cell="A-549") #CCS-1477 IS A KNOWN MYC INHIBITOR
# drugDoseResponseCurve(dataset, drug="Cbp-IN-1", cell="A-549")
# cell_line_sensitivity <- dataset@sensitivity[["profiles"]]
# pert_info <-dataset@sensitivity[["info"]]
# drug_info <- dataset@drug
# # get cmap drugs in this dataset
# drug_info <- drug_info %>% filter((smiles %in% data_cmap$canonical_smiles) |
#                              (inchikey %in% data_cmap$inchi_key) | 
#                              (cid %in% data_cmap$pubchem_cid) | 
#                              (tolower(drugid) %in% tolower(data_cmap$pert_iname)))
# smile_merged <- left_join(drug_info,data_cmap %>% dplyr::select(canonical_smiles,pert_iname,cell_id,sig_id),
#                           by=c('smiles'='canonical_smiles')) %>% filter(!is.na(sig_id))
# inchi_merged <- left_join(drug_info,data_cmap %>% dplyr::select(inchi_key,pert_iname,cell_id,sig_id),
#                           by=c('inchikey'='inchi_key')) %>% filter(!is.na(sig_id))
# cid_merged <- left_join(drug_info,data_cmap %>% dplyr::select(pubchem_cid,pert_iname,cell_id,sig_id),
#                           by=c('cid'='pubchem_cid')) %>% filter(!is.na(sig_id))
# drugid_merged <- left_join(drug_info %>%  mutate(drugid=tolower(drugid)),
#                            data_cmap %>% dplyr::select(pert_iname,cell_id,sig_id) %>%
#                              mutate(pert_iname=tolower(pert_iname)),
#                         by=c('drugid'='pert_iname')) %>% filter(!is.na(sig_id)) %>% dplyr::select(-drugid)
# drugid_merged <- left_join(drug_info,drugid_merged) %>% filter(!is.na(sig_id))
# drugid_merged <- left_join(drugid_merged,data_cmap %>% dplyr::select(pert_iname,sig_id)) %>% filter(!is.na(sig_id))
# drugid_merged <- drugid_merged %>% dplyr::select(all_of(colnames(smile_merged)))
# drug_info <- rbind(smile_merged,inchi_merged,cid_merged,drugid_merged) %>% unique()
# # # drug original with drugbank
# # drug_info_original <- dataset@drug
# # drug_info_original <- drug_info_original %>% filter(smiles %in% rownames(drug_targets_drugbank) | toupper(drugid) %in% toupper(rownames(drug_targets_drugbank)))
# 
# ## Keep perturbations only for those drugs
# # pert_info <- pert_info %>% filter(drugid %in% c(drug_info$drugid,drug_info_original$drugid))
# pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
# pert_info <- pert_info %>% mutate(cellid=gsub('-','',cellid))
# pert_info <- pert_info %>% mutate(cellid=toupper(cellid))
# pert_info <- pert_info %>% filter(cellid %in% data_cmap$cell_id)
# drug_info <- drug_info %>% filter(drugid %in% pert_info$drugid)
# drug_info <- drug_info %>% filter(cell_id %in% pert_info$cellid)
# # drug_info_original <- drug_info_original  %>% filter(drugid %in% pert_info$drugid)
# # pert_info <- pert_info %>% filter(drugid %in% c(drug_info$drugid,drug_info_original$drugid))
# pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
# # filter sensitivity data
# cell_line_sensitivity <- cell_line_sensitivity[rownames(pert_info),]
# gc()
# # filter to keep drugs that are in our data
# cell_line_sensitivity$EC50 <- log10(cell_line_sensitivity$EC50)
# drug_info <- drug_info %>% filter(smiles %in% rownames(drug_targets))
# # drug_info_original <- drug_info_original %>% filter((smiles %in% rownames(drug_targets)) | (drugid %in% rownames(drug_targets)))
# pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
# # pert_info <- pert_info %>% filter(drugid %in% c(drug_info$drugid,drug_info_original$drugid))
# cell_line_sensitivity <- cell_line_sensitivity[rownames(pert_info),]
# cell_line_sensitivity <- cell_line_sensitivity %>% rownames_to_column('exp_id') %>%
#   select(exp_id,EC50) %>% unique()
# # create first cell-line one-hot encoding matrix
# cell_info <- pert_info %>% rownames_to_column('exp_id') %>% select(exp_id,cellid) %>% unique() %>% mutate(value=1)
# cell_info <- cell_info %>% spread('cellid','value')
# cell_info[is.na(cell_info)] <- 0
# 
# pert_info <- pert_info %>% rownames_to_column('exp_id') %>% select(exp_id,drugid,cellid)
# # pert_info_original <- left_join(pert_info,drug_info_original %>% select(drugid,smiles) %>% unique())
# # pert_info_original <- pert_info_original %>% filter(drugid %in% drug_info_original$drugid)
# # pert_info_original <- distinct(pert_info_original)
# # pert_info_original <- pert_info_original %>% mutate(smiles = ifelse(!(smiles %in% rownames(drug_targets)),drugid,smiles))
# # pert_info_original <- distinct(pert_info_original)
# pert_info <- left_join(pert_info,drug_info %>% select(drugid,smiles) %>% unique())
# # pert_info <- rbind(pert_info,pert_info_original)
# pert_info <- left_join(pert_info,as.data.frame(drug_targets) %>% rownames_to_column('smiles'))
# pert_info <- left_join(pert_info,cell_info)
# pert_info <-  pert_info %>% filter(!is.na(smiles)) 
# # get all data
# data <- left_join(pert_info,cell_line_sensitivity)
# data <- data %>% group_by(drugid,cellid) %>% mutate(EC50=mean(EC50)) %>% ungroup()
# data <- data %>% filter(!is.na(EC50))
# data <- data %>% select(-cellid,-exp_id) %>% unique()
# saveRDS(data,'NCI60_2021.rds')

# Skip above and run this
data <- readRDS(paste0('./',set,'.rds'))
  
### Use DTI and cell line info as inputs---------------------------------------
drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')

merged_interactions_all <- merged_interactions_all %>% filter(!is.na(`Prior knowledge`))
merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = 1*(Inferred=='Interaction'))
merged_interactions_all <- merged_interactions_all %>% mutate(`Prior knowledge` = 1*(`Prior knowledge`=='Interaction'))
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
  mutate(model_counts = sum(Inferred==1)) %>% ungroup()

if (type=='threshold'){
  merged_interactions_all <- merged_interactions_all %>% mutate(model_counts=model_counts/no_models)
  # print(all(merged_interactions_all %>% filter(model_counts==0) %>% mutate(logic = (Inferred==`Prior knowledge`)) %>% select(logic)))
  merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = ifelse(model_counts>=th/no_models,1,0))
}else if (type=='frequency'){
  merged_interactions_all <- merged_interactions_all %>% mutate(Inferred=model_counts/no_models)
} else{
}

if (type!='prior') {
  inferred_interactions <- merged_interactions_all %>% select(drug,variable,Inferred)
  inferred_interactions <- distinct(inferred_interactions)
  inferred_interactions <- inferred_interactions %>% group_by(drug) %>% mutate(interactions = sum(Inferred)) %>% ungroup()
  inferred_interactions <- inferred_interactions %>% filter(interactions>0)
  inferred_interactions <- inferred_interactions %>% select(-interactions)
  inferred_interactions <- inferred_interactions %>% spread('variable','Inferred')
  gc()
  drug_targets_inferred <- as.matrix(distinct(inferred_interactions %>% column_to_rownames('drug')))
}
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions_all$drug),
                             which(colnames(drug_targets) %in% merged_interactions_all$variable)]

data_new <- data
data_new <- data_new %>% filter(smiles %in% merged_interactions_all$drug)
# data_new <- data_new %>% filter(smiles %in% rownames(drug_targets))
data_new[,3:(ncol(data_new)-9)] <- apply(data_new[,3:(ncol(data_new)-9)],c(1,2),as.numeric)
changed <- 0
if(type!='prior'){
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
  # data_new[,3:(ncol(data_new)-9)] <- scale(data_new[,3:(ncol(data_new)-9)],scale=F)
}
# hist(as.matrix(data_new[,3:(ncol(data_new)-9)]),10)

### Train ML model with LOOCV procedure----------------------------------
### Call function ###
total_results <- viability_model(data_new,
                                 dataset=set,
                                 lethality_data_path = "./",
                                 model = md,
                                 no_models=50)

mdls <- total_results$mdls
frequency_results <- total_results$frequency_results
res_test_frequency <- total_results$res_test_frequency
res_test_frequency_all <- total_results$res_test_frequency_all
  
saveRDS(frequency_results,paste0('../results/viability_analysis_results/',type,'_results_',toupper(md),'.rds'))
saveRDS(mdls,paste0('../results/viability_analysis_results/','models_',type,'_',toupper(md),'_a549.rds'))
saveRDS(res_test_frequency,paste0('../results/viability_analysis_results/','predictions_A549_',type,'_',toupper(md),'.rds'))
saveRDS(res_test_frequency_all,paste0('../results/viability_analysis_results/','predictions_all_',type,'_',toupper(md),'.rds'))

colnames(frequency_results)[c(3,5)] <- c('LOOCV A549','LOOCV A549 shuffled')
frequency_results <- frequency_results %>% gather('data','r')
frequency_results <- frequency_results %>% filter(!(data %in% c('LOOCV A549','LOOCV A549 shuffled')))
# ggboxplot(frequency_results,x='data',y='r',color='data',add='jitter') + 
#   ylab('pearson`s r')+
#   theme(text = element_text(family = 'Arial',size=24),legend.position = 'none')+
#   stat_compare_means(comparisons = list(c('test','shuffled')),method = 'wilcox.test',
#                      label.y = 0.8,size=6,tip.length = 0.01)

p1 <- ggscatter(res_test_frequency_all%>% filter(data=='test'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 10) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('LOOCV performance of models')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5)) #+
  # geom_point(data = res_test_frequency %>% filter(data=='test'),aes(x=predicted,y=EC50),color='orange',shape=17,size = 4)
p2 <- ggscatter(res_test_frequency_all%>% filter(data=='shuffled'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 10) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('LOOCV performance of shuffled models')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
p <- p1+p2
print(p)
ggsave(paste0('../results/viability_analysis_results/','performance_',type,'_',toupper(md),'_a549.eps'),
       plot=p,
       device = cairo_ps,
       height = 10,
       width = 17,
       units = 'in',
       dpi = 600)
ggsave(paste0('../results/viability_analysis_results/','performance_',type,'_',toupper(md),'_a549.png'),
       plot=p,
       height = 10,
       width = 17,
       units = 'in',
       dpi = 600)


# # Create initial RF model using different number of models----------------
# # Now predict test set with the inferred
# # Do it for multiple thresholds to check. And then choose one
# train_corr <- NULL
# test_corr <- NULL
# train_pval <- NULL
# test_pval <- NULL
# performance_df <- data.frame() 
# no_models <- 50
# for (th in 1:no_models){
#   drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
#   merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
#   #drug_targets_prior <- as.data.frame(drug_targets) %>% rownames_to_column('drug') %>% 
#   #  gather('variable','Prior knowledge',-drug)
#   #merged_interactions_all <- merged_interactions_all %>% select(-`Prior knowledge`)
#   #merged_interactions_all <- left_join(merged_interactions_all,drug_targets_prior)
#   merged_interactions_all <- merged_interactions_all %>% filter(!is.na(`Prior knowledge`))
#   merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = 1*(Inferred=='Interaction'))
#   merged_interactions_all <- merged_interactions_all %>% mutate(`Prior knowledge` = 1*(`Prior knowledge`=='Interaction'))
#   merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
#   merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
#     mutate(model_counts = sum(Inferred==1)) %>% ungroup()
#   merged_interactions_all <- merged_interactions_all %>% mutate(model_counts=model_counts/no_models)
#   # print(all(merged_interactions_all %>% filter(model_counts==0) %>% mutate(logic = (Inferred==`Prior knowledge`)) %>% select(logic)))
#   merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = ifelse(model_counts>=th/no_models,1,0))
#   inferred_interactions <- merged_interactions_all %>% select(drug,variable,Inferred)
#   inferred_interactions <- inferred_interactions %>% unique()
#   inferred_interactions <- inferred_interactions %>% group_by(drug) %>% mutate(interactions = sum(Inferred)) %>% ungroup()
#   inferred_interactions <- inferred_interactions %>% filter(interactions>0)
#   inferred_interactions <- inferred_interactions %>% select(-interactions)
#   inferred_interactions <- inferred_interactions %>% spread('variable','Inferred')
#   gc()
#   drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions_all$drug),
#                                which(colnames(drug_targets) %in% merged_interactions_all$variable)]
#   drug_targets <- drug_targets[which(rownames(drug_targets) %in% drug_info$smiles),]
#   
#   data_new <- data
#   drug_targets_inferred <- as.matrix(distinct(inferred_interactions %>% column_to_rownames('drug')))
#   data_new <- data_new %>% filter(smiles %in% merged_interactions_all$drug)
#   changed <- 0
#   for (i in 1:nrow(data_new)){
#     drug <- data_new$smiles[i]
#     if (drug %in% rownames(drug_targets_inferred)){
#       for (j in 3:(ncol(data_new)-9)){
#         target <-   colnames(data_new)[j]
#         val <- drug_targets_inferred[drug,target]
#         if (val!=data_new[i,j]){
#           changed <- changed+1
#         }
#         data_new[i,j] <- val
#       }
#     }
#   }
#   train_corr_loocv <- NULL
#   test_corr_all_loocv <- NULL
#   test_preds_loocv <- NULL
#   test_true_loocv <-  NULL
#   for (loocv_index in 1:nrow(a549_data)){
#     loocv_smile <- a549_data$smiles[loocv_index]
#     # find test_index
#     test_index <- which(data_new$smiles==loocv_smile)
#     test_index_a549 <- which(data_new$A549==1 & data_new$smiles==loocv_smile)
#     train_exp_ids <- pert_info$exp_id[-test_index]
#     # Subset your data into training and testing sets based on the indices
#     train_data_new <- data_new[-test_index, ]
#     train_drugs <- train_data_new$smiles
#     train_data_new <- train_data_new %>% select(-drugid,-smiles)
#     rownames(train_data_new) <- NULL
#     # train_data_new <- train_data_new %>% column_to_rownames('exp_id')
#     test_data_new <- data_new[test_index, ]
#     test_drugs <- test_data_new$smiles
#     test_data_new <- test_data_new %>% select(-drugid,-smiles)
#     rownames(test_data_new) <- NULL
#     # test_data_new <- test_data_new %>% column_to_rownames('exp_id')
#     ### keep separately only LOOCV A549
#     test_data_a549 <- data_new[test_index_a549, ]
#     test_drugs <- test_data_a549$smiles
#     test_data_a549 <- test_data_a549 %>% select(-drugid,-smiles)
#     rownames(test_data_a549) <- NULL
#     gc()
#     
#     # build new model
#     ctrl <- trainControl(method = "cv", number = 10)
#     mdl <- train(EC50 ~ ., data = train_data_new, method = "rf", trControl = ctrl,trace=F)
#     y_train_new <- predict(mdl,newdata = train_data_new)
#     train_corr_loocv[loocv_index] <- cor(y_train_new,train_data_new$EC50)
#     # first correlation for all test samples
#     y_new <- predict(mdl,newdata = test_data_new)
#     test_corr_all_loocv[loocv_index] <- cor(y_new,test_data_new$EC50)
#     # only A549
#     y_new <- predict(mdl,newdata = test_data_a549)
#     test_preds_loocv <- c(test_preds_loocv,y_new)
#     test_true_loocv <- c(test_true_loocv,test_data_a549$EC50)
#   }
#   test_corr[th] <- cor(test_preds_loocv,test_true_loocv)
#   train_corr[th] <- mean(train_corr_loocv)
#   performance_df <- rbind(performance_df,
#                           data.frame(train=train_corr_loocv,
#                                      test=test_corr_all_loocv,
#                                      'LOOCV A549'=test_corr[th],
#                                      frequency = rep(th,length(train_corr_loocv))/no_models))
#   # saveRDS(performance_df,'../results/viability_results_LOOCV_analysis.rds')
#   message(paste0('Finished model ',th))
# }
# colnames(performance_df)[3] <- "LOOCV A549"
# # saveRDS(performance_df,'../results/viability_results_LOOCV_analysis.rds')
# # plot(seq(1:no_models)/no_models,train_corr)
# # plot(seq(1:no_models)/no_models,test_corr)
# # results_ensembles <- data.frame(frequency = seq(1:no_models)/no_models,
# #                                 train= train_corr,
# #                                 test=test_corr)
# # results_ensembles <- results_ensembles %>% gather('data','r',-frequency)
# # saveRDS(results_ensembles,'../results/A549_viability_results_analysis.rds')
# # results_ensembles <- readRDS('../results/A549_viability_results_analysis.rds')
# results_plot_all <- performance_df %>% gather('data','r',-frequency) %>% 
#   group_by(data,frequency) %>% mutate(mu = mean(r)) %>% mutate(se = sd(r)/sqrt(n())) %>%
#   ungroup()
# ggplot(results_plot_all,
#        aes(x=frequency,y=mu,color=data)) + 
#   geom_point() + geom_line()+
#   geom_errorbar(data = results_plot_all %>% filter(data!='LOOCV A549'),
#                 aes(ymax = mu + se, ymin = mu-se)) +
#   geom_hline(yintercept =0,linetype = 'dashed',color='black',lwd=1)+
#   ggtitle('Performance of Random Forest models in predicting lethality') +
#   ylab('pearson`s r') +
#   theme(title  = element_text(hjust=0.5))+
#   theme_pubr(base_family = 'Arial',
#              base_size = 24)+
#   theme(text=element_text(family = 'Arial',size=24),
#         plot.title  = element_text(hjust=0.5),
#         legend.position = 'right' )
# ggsave('../article_supplementary_info/supple_fig11.png',
#        scale = 1,
#        width = 16,
#        height = 12,
#        units = "in",
#        dpi = 600)