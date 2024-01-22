library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(PharmacoGx)
library(caret)
source("viabilityModeling.R")

# Initialize
set <- 'NCI60_2021'
md <- 'rf'
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
# PRISM has sensitivity data
available_datasets <- availablePSets()
# Load data through PharmacoDB
dataset <- PharmacoGx::downloadPSet('NCI60_2021')
## Plot Drug Dose response curves, using the same names for compounds and cell lines as PharmacoDB
drugDoseResponseCurve(dataset, drug="Lestaurtinib", cell="A-549") #CCS-1477 IS A KNOWN MYC INHIBITOR
drugDoseResponseCurve(dataset, drug="Cbp-IN-1", cell="A-549")
cell_line_sensitivity <- dataset@sensitivity[["profiles"]]
pert_info <-dataset@sensitivity[["info"]]
drug_info <- dataset@drug
# get cmap drugs in this dataset
drug_info <- drug_info %>% filter((smiles %in% data_cmap$canonical_smiles) |
                             (inchikey %in% data_cmap$inchi_key) |
                             (cid %in% data_cmap$pubchem_cid) |
                             (tolower(drugid) %in% tolower(data_cmap$pert_iname)))
smile_merged <- left_join(drug_info,data_cmap %>% dplyr::select(canonical_smiles,pert_iname,cell_id,sig_id),
                          by=c('smiles'='canonical_smiles')) %>% filter(!is.na(sig_id))
inchi_merged <- left_join(drug_info,data_cmap %>% dplyr::select(inchi_key,pert_iname,cell_id,sig_id),
                          by=c('inchikey'='inchi_key')) %>% filter(!is.na(sig_id))
cid_merged <- left_join(drug_info,data_cmap %>% dplyr::select(pubchem_cid,pert_iname,cell_id,sig_id),
                          by=c('cid'='pubchem_cid')) %>% filter(!is.na(sig_id))
drugid_merged <- left_join(drug_info %>%  mutate(drugid=tolower(drugid)),
                           data_cmap %>% dplyr::select(pert_iname,cell_id,sig_id) %>%
                             mutate(pert_iname=tolower(pert_iname)),
                        by=c('drugid'='pert_iname')) %>% filter(!is.na(sig_id)) %>% dplyr::select(-drugid)
drugid_merged <- left_join(drug_info,drugid_merged) %>% filter(!is.na(sig_id))
drugid_merged <- left_join(drugid_merged,data_cmap %>% dplyr::select(pert_iname,sig_id)) %>% filter(!is.na(sig_id))
drugid_merged <- drugid_merged %>% dplyr::select(all_of(colnames(smile_merged)))
drug_info <- rbind(smile_merged,inchi_merged,cid_merged,drugid_merged) %>% unique()
# # drug original with drugbank
# drug_info_original <- dataset@drug
# drug_info_original <- drug_info_original %>% filter(smiles %in% rownames(drug_targets_drugbank) | toupper(drugid) %in% toupper(rownames(drug_targets_drugbank)))

## Keep perturbations only for those drugs
# pert_info <- pert_info %>% filter(drugid %in% c(drug_info$drugid,drug_info_original$drugid))
pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
pert_info <- pert_info %>% mutate(cellid=gsub('-','',cellid))
pert_info <- pert_info %>% mutate(cellid=toupper(cellid))
pert_info <- pert_info %>% filter(cellid %in% data_cmap$cell_id)
drug_info <- drug_info %>% filter(drugid %in% pert_info$drugid)
drug_info <- drug_info %>% filter(cell_id %in% pert_info$cellid)
# drug_info_original <- drug_info_original  %>% filter(drugid %in% pert_info$drugid)
# pert_info <- pert_info %>% filter(drugid %in% c(drug_info$drugid,drug_info_original$drugid))
pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
# filter sensitivity data
cell_line_sensitivity <- cell_line_sensitivity[rownames(pert_info),]
gc()
# filter to keep drugs that are in our data
cell_line_sensitivity$EC50 <- log10(cell_line_sensitivity$EC50)
drug_info <- drug_info %>% filter(smiles %in% rownames(drug_targets))
# drug_info_original <- drug_info_original %>% filter((smiles %in% rownames(drug_targets)) | (drugid %in% rownames(drug_targets)))
pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
# pert_info <- pert_info %>% filter(drugid %in% c(drug_info$drugid,drug_info_original$drugid))
cell_line_sensitivity <- cell_line_sensitivity[rownames(pert_info),]
cell_line_sensitivity <- cell_line_sensitivity %>% rownames_to_column('exp_id') %>%
  select(exp_id,EC50) %>% unique()
# create first cell-line one-hot encoding matrix
cell_info <- pert_info %>% rownames_to_column('exp_id') %>% select(exp_id,cellid) %>% unique() %>% mutate(value=1)
cell_info <- cell_info %>% spread('cellid','value')
cell_info[is.na(cell_info)] <- 0

pert_info <- pert_info %>% rownames_to_column('exp_id') %>% select(exp_id,drugid,cellid)
# pert_info_original <- left_join(pert_info,drug_info_original %>% select(drugid,smiles) %>% unique())
# pert_info_original <- pert_info_original %>% filter(drugid %in% drug_info_original$drugid)
# pert_info_original <- distinct(pert_info_original)
# pert_info_original <- pert_info_original %>% mutate(smiles = ifelse(!(smiles %in% rownames(drug_targets)),drugid,smiles))
# pert_info_original <- distinct(pert_info_original)
pert_info <- left_join(pert_info,drug_info %>% select(drugid,smiles) %>% unique())
# pert_info <- rbind(pert_info,pert_info_original)
pert_info <- left_join(pert_info,as.data.frame(drug_targets) %>% rownames_to_column('smiles'))
pert_info <- left_join(pert_info,cell_info)
pert_info <-  pert_info %>% filter(!is.na(smiles))
# get all data
data <- left_join(pert_info,cell_line_sensitivity)
data <- data %>% group_by(drugid,cellid) %>% mutate(EC50=mean(EC50)) %>% ungroup()
data <- data %>% filter(!is.na(EC50))
data <- data %>% select(-cellid,-exp_id) %>% unique()
# saveRDS(data,'NCI60_2021.rds')

# Skip above and run this
data <- readRDS(paste0('./',set,'.rds'))
  
### Use DTI and cell line info as inputs---------------------------------------
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
# merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
#   mutate(model_counts = sum(Inferred==1)) %>% ungroup()
if (type=='threshold'){
  # merged_interactions_all <- merged_interactions_all %>% mutate(model_counts=model_counts/no_models)
  # print(all(merged_interactions_all %>% filter(model_counts==0) %>% mutate(logic = (Inferred==`Prior knowledge`)) %>% select(logic)))
  merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = ifelse((consensus_inferrence==3) & (mean_frequency>=th/no_models),1,0))
  merged_interactions_all <- merged_interactions_all %>% select(-mean_frequency,-consensus_inferrence)
  merged_interactions_all <- distinct(merged_interactions_all)
}else if (type=='frequency'){
  merged_interactions_all <- merged_interactions_all %>% mutate(Inferred=mean_frequency)
  merged_interactions_all <- merged_interactions_all %>% select(-mean_frequency,-consensus_inferrence)
  merged_interactions_all <- distinct(merged_interactions_all)
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
                                 colnames(drug_targets),
                                 dataset=set,
                                 lethality_data_path = "./",
                                 model = md,
                                 no_models=50)

mdls <- total_results$mdls
# final_results <- total_results$final_results
res_test <- total_results$res_test
res_test_all <- total_results$res_test_all

#saveRDS(final_results,paste0(type,'_results_',toupper(md),'_',toupper(set),'.rds'))
saveRDS(mdls,paste0('models_',type,'_',toupper(md),'_',toupper(set),'.rds'))
saveRDS(res_test,paste0('predictions_A549_',type,'_',toupper(md),'_',toupper(set),'.rds'))
saveRDS(res_test_all,paste0('predictions_all_',type,'_',toupper(md),'_',toupper(set),'.rds'))

# colnames(final_results)[c(3,5)] <- c('LOOCV A549','LOOCV A549 shuffled')
# final_results <- final_results %>% gather('data','r')
# final_results <- final_results %>% filter(!(data %in% c('LOOCV A549','LOOCV A549 shuffled')))
# ggboxplot(final_results,x='data',y='r',color='data',add='jitter') +
#   ylab('pearson`s r')+
#   theme(text = element_text(family = 'Arial',size=24),legend.position = 'none')+
#   stat_compare_means(comparisons = list(c('test','shuffled')),method = 'wilcox.test',
#                      label.y = 0.8,size=6,tip.length = 0.01)

p1 <- ggscatter(res_test_all %>% filter(data=='test'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 10) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('LOOCV performance of models')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5)) #+
# geom_point(data = res_test %>% filter(data=='test'),aes(x=predicted,y=EC50),color='orange',shape=17,size = 4)
p2 <- ggscatter(res_test_all %>% filter(data=='shuffled'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 10) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('LOOCV performance of shuffled models')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=24,family='Arial'),
        plot.title = element_text(hjust = 0.5))
p <- p1+p2
print(p)
ggsave(paste0('../results/viability_analysis_results/','performance_',type,'_',toupper(md),'_',toupper(set),'.eps'),
       plot=p,
       device = cairo_ps,
       height = 10,
       width = 17,
       units = 'in',
       dpi = 600)
ggsave(paste0('../results/viability_analysis_results/','performance_',type,'_',toupper(md),'_',toupper(set),'.png'),
       plot=p,
       height = 10,
       width = 17,
       units = 'in',
       dpi = 600)
