library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(viridis)
library(patchwork)

### Load drug-target interaction data in A375 and performance of models------------------------
cmap_drugs <- readRDS('../preprocessing/preprocessed_data/l1000_drugs_with_targets_all.rds')
no_models <- 50
drug_targets <- readRDS('../preprocessing/preprocessed_data/drug_targets_space.rds')
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]
cmap_drugs <- cmap_drugs %>% filter(canonical_smiles %in% rownames(drug_targets))
pert_info <- read.delim(file = "../data/GSE92742_Broad_LINCS_pert_info.txt")
pert_info <- pert_info %>% filter(pert_type == "trt_cp")
pert_info <- pert_info %>% filter(canonical_smiles %in% cmap_drugs$canonical_smiles)
pert_info <- pert_info %>% select(canonical_smiles,pert_iname) %>% unique()
cmap_drugs <- left_join(pert_info,cmap_drugs)
lestaurtinib_moas <- c('FLT3 inhibitor', 'growth factor receptor inhibitor', 'JAK inhibitor')
cmap_drugs <- cmap_drugs %>% mutate(MOA = str_split(MOA,pattern = ", ")) %>% unnest(MOA) %>% unique()

top_moas <- cmap_drugs %>% group_by(MOA) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_moas <- top_moas %>% select(MOA,counts) %>% unique()
top_moas <- top_moas %>% arrange(-counts)
top_moas <- top_moas[1:10,]

top_diseases <- cmap_drugs %>% group_by(Disease.Area) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_diseases <- cmap_drugs %>% group_by(Disease.Area) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_diseases <- top_diseases %>% select(Disease.Area,counts) %>% unique()
top_diseases <- top_diseases %>% filter(Disease.Area!='')
top_diseases <- top_diseases %>% arrange(-counts)
top_diseases <- top_diseases[1:10,]

### Load predictions and true values--------------------------------
files_preds <- list.files('../results/A375_ensembles/preds/')
files_preds <- files_preds[grep('.csv',files_preds)]
files_preds <- files_preds[grep('mean_val',files_preds)]
# Load validation correlation data
df_preds_val <- data.frame()
for (file in files_preds){
  cell <- str_split_fixed(file,'_',3)[1,2]
  file <- paste0('../results/A375_ensembles/preds/',file)
  tmp <- data.table::fread(file,header = T)
  colnames(tmp)[1] <- 'sample'
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  
  df_preds_val <- rbind(df_preds_val,tmp)
}
df_preds_val <- df_preds_val %>% gather('TF','prediction',-cell,-sample)
TFoutput <- read.delim('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv')%>%
  column_to_rownames('X') %>% rownames_to_column('sample')
gc()
TFoutput <- TFoutput %>% gather('TF','activity',-sample) 
TFoutput <- TFoutput %>% filter(sample %in% df_preds_val$sample)
df_preds_val <- left_join(df_preds_val,TFoutput)

### Combine performance,activity and MOA info----------------------------
sig_info <- read.delim("../data/GSE92742_Broad_LINCS_sig_info.txt")
sig_info <- sig_info %>% filter(sig_id %in% df_preds_val$sample)
sig_info <- sig_info %>% select(sig_id,pert_id,pert_iname,pert_dose,pert_time) %>% unique()
pert_info <- read.delim(file = "../data/GSE92742_Broad_LINCS_pert_info.txt")
pert_info <- pert_info %>% filter(pert_type == "trt_cp")
pert_info <- pert_info %>% filter(canonical_smiles %in% cmap_drugs$canonical_smiles)
pert_info <- pert_info %>% select(canonical_smiles,pert_iname,pert_id) %>% unique()
sig_info <- left_join(pert_info,sig_info)
cmap_drugs <- left_join(cmap_drugs,sig_info)
cmap_drugs <- cmap_drugs %>% filter(sig_id %in% df_preds_val$sample)
data_all <- left_join(cmap_drugs,df_preds_val,by=c('sig_id'='sample'))
data_all <- data_all %>% filter(!is.na(prediction))
print(length(unique(data_all$canonical_smiles)))

### Start analyzing results------------------------------
data_all_plot_moas <- data_all %>% group_by(MOA,TF) %>% mutate(r = cor(prediction,activity)) %>% ungroup()
data_all_plot_moas <- data_all_plot_moas %>% filter(!is.na(r))
tmp <- data_all_plot_moas %>% group_by(MOA) %>% mutate(mean_r = mean(r)) %>% ungroup()
tmp <- tmp %>% select(MOA,mean_r) %>% unique()
data_all_plot_moas$MOA <- factor(data_all_plot_moas$MOA,levels = tmp$MOA[order(tmp$mean_r)])
ggboxplot(data_all_plot_moas %>% filter(MOA %in% c(top_moas$MOA,lestaurtinib_moas)) %>% select(MOA,TF,r) %>% unique(),
          x='MOA',y='r',color='MOA',add='jitter')+
  theme(text=element_text(family = 'Arial',size=24),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = 'none')+
  coord_flip()

## Check per disease area
data_all_plot_diseases <- data_all %>% group_by(Disease.Area,TF) %>% mutate(r = cor(prediction,activity)) %>% ungroup()
data_all_plot_diseases <- data_all_plot_diseases %>% filter(!is.na(r))
tmp <- data_all_plot_diseases %>% group_by(Disease.Area) %>% mutate(mean_r = mean(r)) %>% ungroup()
tmp <- tmp %>% select(Disease.Area,mean_r) %>% unique()
data_all_plot_diseases$Disease.Area <- factor(data_all_plot_diseases$Disease.Area,levels = tmp$Disease.Area[order(tmp$mean_r)])
ggboxplot(data_all_plot_diseases %>% filter(Disease.Area %in% top_diseases$Disease.Area) %>% select(Disease.Area,TF,r) %>% unique(),
          x='Disease.Area',y='r',color='Disease.Area',add='jitter')+
  theme(text=element_text(family = 'Arial',size=24),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = 'none')+
  coord_flip()
