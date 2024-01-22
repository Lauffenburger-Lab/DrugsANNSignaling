library(tidyverse)
library(readr)
library(GGally) 
library(ggplot2)
library(ggpubr)
library(PharmacoGx)
library(caret)
library(pheatmap)
library(ggrepel)

# Load cmap data---------------
drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]
data_cmap <- readRDS("../preprocessing/preprocessed_data/all_cmap_sigs_with_pert_info.rds")
TF_activities_2 <- read.delim('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv',row.names = 1)
data_cmap <- data_cmap %>% filter(sig_id %in% rownames(TF_activities_2))
data_cmap <- data_cmap %>% filter(canonical_smiles %in% rownames(drug_targets))
gc()
TF_activities_2 <- NULL
merged_interactions_all <- NULL
gc()

# Load lethality/viability datasets----------------------
columns <- c('drugid','cid','inchikey','smiles')
available_datasets <- availablePSets()
available_datasets <- available_datasets %>% filter(!(`PSet Name` %in% c('GDSC_2020(v1-8.2)','GBM_scr2','GBM_scr3','PDTX_2019','BeatAML_2018','Tavor_2020')))
datasets <- available_datasets$`PSet Name`
names <- available_datasets$`Dataset Name`
print(datasets)
data_all <- data.frame()
# data_all <- readRDS('../results/viability_data_all.rds')
for (j in 1:length(datasets)){
  # message(paste0('Begun dataset: ',datasets[j]))
  dataset <- PharmacoGx::downloadPSet(datasets[j])
  cell_line_sensitivity <- dataset@sensitivity[["profiles"]]
  pert_info <-dataset@sensitivity[["info"]]
  drug_info <- dataset@drug
  if (any(!(columns %in% colnames(drug_info)))){
    cols_add <- columns[which(!(columns %in% colnames(drug_info)))]
    for (column in cols_add){
      drug_info[[column]] <- NA
    }
  }
  drug_info <- drug_info %>% select(all_of(columns)) %>% unique()
  drug_info <- drug_info %>% mutate(cid = as.character(cid))
  if ('ec50' %in% colnames(cell_line_sensitivity)){
    ind <- which(colnames(cell_line_sensitivity)=='ec50')
    colnames(cell_line_sensitivity)[ind] <- 'EC50'
  }
  ## Keep perturbations only for those drugs
  pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
  pert_info <- pert_info %>% mutate(cellid=gsub('-','',cellid))
  pert_info <- pert_info %>% mutate(cellid=toupper(cellid))
  pert_info <- pert_info %>% filter(cellid %in% data_cmap$cell_id)
  drug_info <- drug_info %>% filter(drugid %in% pert_info$drugid)
  pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
  # filter sensitivity data
  cell_line_sensitivity <- cell_line_sensitivity[rownames(pert_info),]
  gc()
  if (!('exp_id' %in% colnames(pert_info))){
    pert_info <- pert_info %>% rownames_to_column('exp_id')
    rownames(pert_info) <- pert_info$exp_id
  }
  # filter to keep drugs that are in our data
  cell_line_sensitivity$EC50 <- log10(cell_line_sensitivity$EC50)
  pert_info <- pert_info %>% filter(drugid %in% drug_info$drugid)
  cell_line_sensitivity <- cell_line_sensitivity[rownames(pert_info),]
  cell_line_sensitivity <- cell_line_sensitivity %>% rownames_to_column('exp_id') %>%
    select(exp_id,EC50) %>% unique()
  # create first cell-line one-hot encoding matrix
  cell_info <- pert_info %>% select(exp_id,cellid) %>% unique()
  pert_info <- pert_info %>% select(exp_id,drugid,cellid)
  pert_info <- left_join(pert_info,drug_info %>% unique())
  pert_info <- left_join(pert_info,cell_info)
  
  # get all data
  data <- left_join(pert_info,cell_line_sensitivity)
  data <- data %>% group_by(drugid,cellid) %>% mutate(EC50=mean(EC50,na.rm=T)) %>% ungroup()
  data <- data %>% filter(!is.na(EC50))
  data <- distinct(data)
  
  # merge with other datasets
  data_all <- rbind(data_all,data %>% mutate(dataset = names[j]))
  # saveRDS(data_all,'../results/viability_data_all.rds')
  message(paste0('Finished dataset: ',datasets[j]))
}
# saveRDS(data_all,'../results/viability_data_all.rds')

# data_all <- rbind(data1 %>% select(drugid,cellid,EC50) %>% unique() %>% mutate(id = paste0(cellid,'_',drugid)) %>% mutate(dataset = 'data1'),
#                   data2 %>% select(drugid,cellid,EC50) %>% unique() %>% mutate(id = paste0(cellid,'_',drugid)) %>% mutate(dataset = 'data2'))
data_all <- data_all %>% select(drugid,cellid,EC50,dataset) %>% unique() %>% mutate(id = paste0(cellid,'_',drugid))
data_all <- data_all %>% select(id,dataset,EC50) %>% spread('dataset','EC50')
data_all <- data_all %>% column_to_rownames('id')

ggpairs(data_all,title  = "Scatterplot Matrix for multiple screening datasets",axisLabels = "show") +
  xlab('log10(EC50)') + ylab('log10(EC50)') +
  theme(text = element_text(family = 'Arial',size=20),
        plot.title = element_text(hjust = 0.5))
ggsave('../article_supplementary_info/drug_screening_correlation.png',
       width = 12,
       height = 12,
       units = 'in',
       dpi=600)

### See only common experiments and keep them as features
X <- data_all %>% rownames_to_column('id') %>% gather('dataset','EC50',-id) %>% spread('id','EC50')
count_nas <- function(x){
  return(length(which(is.na(x))))
}
X <- X %>% column_to_rownames('dataset')
inds <- apply(X,2,count_nas)
print(min(inds))
hist(inds)
inds_keep <- which(inds<=1)
X <- X[,inds_keep]
row_na <- apply(X,1,count_nas)
min(row_na)
hist(row_na)
p <- 0.2 * ncol(X)
rows <- which(row_na<=p)
X <- X[rows,]
medians <- apply(X,2,median,na.rm=T)
for (i in 1:nrow(X)){
  j <- which(is.na(X[i,]))
  X[i,j] <- medians[j]
}

#pca plot
pca <- prcomp(X,center = T,scale = T)
df <- pca$x[,1:2]
df <- as.data.frame(df) %>% rownames_to_column('dataset')
ggplot(df,aes(x=PC1,y=PC2,color=dataset)) + geom_point(size=3) + 
  geom_label_repel(aes(label = dataset),
                   size=6,
                   box.padding   = 0.75, 
                   point.padding = 0.1,
                   max.overlaps = 40) + 
  theme_minimal(base_family = 'Arial',base_size = 22) +
  theme(text = element_text(family = 'Arial',size=22),
        legend.position = 'top')
ggsave('../article_supplementary_info/viability_datasets_pca.eps',
       device = cairo_ps,
       width = 9,
       height = 9,
       units = 'in',
       dpi = 600)
# clustermap
pheatmap(X)
setEPS()
postscript('../article_supplementary_info/viability_datasets_clustermap.eps',width = 9,height = 9)
pheatmap(X)
dev.off()

### Save datasets to use for predicting viability-----------------------------------------
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

available_datasets <- availablePSets()
available_datasets <- available_datasets %>% filter(!(`PSet Name` %in% c('GDSC_2020(v1-8.2)','GBM_scr2','GBM_scr3',
                                                                         'PDTX_2019','BeatAML_2018','Tavor_2020',
                                                                         'UHNBreast_2019')))
datasets <- available_datasets$`PSet Name`
names <- available_datasets$`Dataset Name`
print(datasets)
j <- 1
for (pset in datasets){
  # PRISM has sensitivity data
  available_datasets <- availablePSets()
  # Load data through PharmacoDB
  dataset <- PharmacoGx::downloadPSet(pset)
  ## Plot Drug Dose response curves, using the same names for compounds and cell lines as PharmacoDB
  cell_line_sensitivity <- dataset@sensitivity[["profiles"]]
  pert_info <-dataset@sensitivity[["info"]]
  drug_info <- dataset@drug
  drug_info$cid <- as.character(drug_info$cid)
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
  if (!('exp_id' %in% colnames(pert_info))){
    cell_info <- pert_info %>% rownames_to_column('exp_id') %>% select(exp_id,cellid) %>% unique() %>% mutate(value=1) 
  }else{
    cell_info <- pert_info %>% select(exp_id,cellid) %>% unique() %>% mutate(value=1)
  }
  cell_info <- cell_info %>% spread('cellid','value')
  cell_info[is.na(cell_info)] <- 0
  
  if (!('exp_id' %in% colnames(pert_info))){
    pert_info <- pert_info %>% rownames_to_column('exp_id') %>% select(exp_id,drugid,cellid)
  }else{
    pert_info <- pert_info %>% select(exp_id,drugid,cellid)
  }
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
  saveRDS(data,paste0(pset,'.rds'))
  print(paste0('Finished ',names[j]))
  j <- j+1
}
