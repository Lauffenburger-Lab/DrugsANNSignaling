library(tidyverse)
cmap_drugs_unnested <- readRDS('preprocessed_data/l1000_drugs_with_targets_unnested_exemplar.rds')
pkn <- data.table::fread('preprocessed_data/PKN/pkn.tsv')
annot <- read.csv('../data/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t')
annot <- annot %>% select(Gene.names...primary..,Entry)
colnames(annot)[1] <- 'Target'
cmap_drugs_unnested <- left_join(cmap_drugs_unnested,annot)
cmap_drugs_unnested <- cmap_drugs_unnested %>% mutate(Entry=ifelse(Target=='None','None',Entry))
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(!is.na(Entry))
dorothea <- read.delim('preprocessed_data/TF_activities/l1000_allgenes_lvl3_tfs.tsv',row.names = 1)
print(length(which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source))))) # 5 TFs not in PKN for lands, 9 for all
dorothea <- dorothea[,which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source)))]
write.table(dorothea, file = 'preprocessed_data/TF_activities/Trimmed_l1000_allgenes_lvl3_tfs.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
tfs_targeted <- unique(cmap_drugs_unnested$Entry[which(cmap_drugs_unnested$Entry %in% colnames(dorothea))])
write.table(data.frame('Entry'=tfs_targeted), file = 'preprocessed_data/TF_activities/tfs_targetd_alls_genes_lvl3.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
pkn_tfs_targeted <- pkn %>% filter((source %in% tfs_targeted) | (target %in% tfs_targeted))
write.table(pkn_tfs_targeted, file = 'preprocessed_data/PKN/L1000_latest_Add_lvl3.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
write.table(cmap_drugs_unnested %>% select('source'='canonical_smiles','target'='Entry'), 
            file = 'preprocessed_data/PKN/L1000_lvl3_DT.tsv', quote=FALSE, sep = "\t", row.names = F)


#### After trimming net ###
pkn <- data.table::fread(('preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv'))
nodes <- unique(c(pkn$source,pkn$target))
cmap_drugs <- cmap_drugs_unnested %>% filter((Entry %in% nodes) | (Entry=='None') | (canonical_smiles=='CS(C)=O.CS(C)=O'))
print(length(unique(cmap_drugs$canonical_smiles))) # 893
print(length(unique(cmap_drugs$Target))) # 461
saveRDS(cmap_drugs,'preprocessed_data/TrimmedFinal_lvl3_allgenes_all_conditions.rds')
dorothea <- read.delim('preprocessed_data/TF_activities/Trimmed_l1000_allgenes_lvl3_tfs.tsv',row.names = 1)
print(length(which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source)))))
dorothea <- dorothea[,which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source)))]
print(length(which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source))))) # 101 out of 106 for alls
write.table(dorothea, file = 'preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)


### Split data according to cells ###
pkn <- data.table::fread(('preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv'))
data <- readRDS("preprocessed_data/all_cmap_sigs_with_pert_info.rds")
data <- data %>% filter(is_exemplar==1)  %>%
  filter(pert_itime %in% c("6 h"))
data <- data %>% dplyr::select(sig_id,pert_type,pert_idose,pert_itime,cell_id,pert_id,canonical_smiles) %>% unique()
gc()
data_merged <- data

cmap_drugs <- readRDS('preprocessed_data/TrimmedFinal_lvl3_allgenes_all_conditions.rds')
dorothea <- read.delim('preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv',row.names = 1)
data_merged <-  data_merged %>% filter(sig_id %in% rownames(dorothea))
cmap_ctrls <- data_merged %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt') %>% unique()
cmap_ctrls <- cmap_ctrls %>% mutate(canonical_smiles = ifelse(pert_id=='DMSO','CS(C)=O',
                                                              ifelse(pert_id=='H2O','O',
                                                                     ifelse(pert_id=='PBS','Ctrl_PBS','Ctrl_UnTrt'))))
cmap_ctrls <- cmap_ctrls %>% mutate(Entry=ifelse(pert_id=='DMSO',c('Q01344', 'Q8WXI7', 'P01106'),'None')) %>% unnest(Entry)
cmap_ctrls <- cmap_ctrls  %>%  mutate(Target=ifelse(Entry=='Q01344','IL5RA',ifelse(Entry=='P01106','MYC',ifelse(Entry=='Q8WXI7','MUC16',Entry)))) 
cmap_ctrls <- cmap_ctrls %>% mutate(MOA=NA,Phase=NA,Disease.Area='control',pubchem_cid=NA)


cmap_drugs <- left_join(cmap_drugs,data_merged)
cmap_drugs <- cmap_drugs %>% filter(sig_id %in% rownames(dorothea)) %>% filter(!is.na(sig_id))
cmap_ctrls <- cmap_ctrls %>% select(colnames(cmap_drugs))

cmap_drugs <- rbind(cmap_drugs,cmap_ctrls)
cmap_drugs <- cmap_drugs %>% filter(!is.na(Entry))
cmap_drugs <- cmap_drugs %>% filter(Entry!='')
cmap_drugs <- cmap_drugs %>% filter(Entry!=' ')
cmap_drugs <- cmap_drugs %>% filter((Entry %in% c(pkn$source,pkn$target)) | (Entry=='None')) %>% unique()
cmap_drugs <- cmap_drugs %>% filter(sig_id %in% rownames(dorothea))
saveRDS(cmap_drugs,'preprocessed_data/TrimmedFinal_allgenes_lvl3_all_conditions_with_targets.rds')

### Load conditions ###----------------------------
library(tidyverse)
cmap_drugs_unnested <- readRDS('preprocessed_data/TrimmedFinal_allgenes_lvl3_all_conditions_with_targets.rds')

### Find 5 cell-lines with most common drugs
cells <- unique(cmap_drugs_unnested$cell_id)
summary <- cmap_drugs_unnested %>% group_by(cell_id) %>% summarise(n_distinct(canonical_smiles))
summary <- summary %>% filter(`n_distinct(canonical_smiles)`>=400)
bestSample <-  unique(summary$cell_id)
print(bestSample)

cell_list <-  NULL
for (j in 1:length(bestSample)){
  cell <- cmap_drugs_unnested %>% filter(cell_id==bestSample[j]) %>% unique()
  cell_list[[j]] <-  unique(cell$canonical_smiles)
}

common <- Reduce(intersect,cell_list)

cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(cell_id %in% bestSample) %>% unique()
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter((Entry %in% c(pkn$source,pkn$target)) | (Entry=='None')) %>% unique()
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(canonical_smiles %in% common)
saveRDS(cmap_drugs_unnested,'preprocessed_data/TrimmedFinal_lvl3_conditions_for_pairedcells.rds')

#### Split sets------------------------
cmap_drugs_unnested <- readRDS('preprocessed_data/TrimmedFinal_lvl3_conditions_for_pairedcells.rds')
uniq_cells <- unique(cmap_drugs_unnested$cell_id)

### Make drug target matrix and cell matrix ###
uniq_targets <- unique(cmap_drugs_unnested$Entry)

binary_classes <- function(sig,df,targets,type=c('smiles','sig_id'),tar=c('Entry','cell_id')){
  if (length(tar)>1){
    tar <- tar[1]
  }
  if (type=='smiles'){
    df <- df %>% select(canonical_smiles,tar) %>% unique()
    df <- df %>% filter(canonical_smiles==sig)
  }else{
    df <- df %>% select(sig_id,tar)  %>% unique()
    df <- df %>% filter(sig_id==sig)
  }
  m <- as.matrix(df[,2])
  ind <- which(targets %in% m)
  n <- targets
  targets <- rep(0,length(targets))
  targets[ind] <- 1
  names(targets) <- n
  return(targets)
}

x <- t(sapply(unique(cmap_drugs_unnested$canonical_smiles),binary_classes,cmap_drugs_unnested,uniq_targets))
#x <- x[,-which(colnames(x)=='None')]
write.table(x,'preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

xc <- t(sapply(unique(cmap_drugs_unnested$sig_id),binary_classes,cmap_drugs_unnested,uniq_cells,type='sig_id',tar='cell_id'))
write.table(xc,'preprocessed_data/TrainingValidationData/L1000_lvl3-conditions_cells.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

print(unique(cmap_drugs_unnested$pert_idose))
cmap_drugs_unnested$pert_idose[which(cmap_drugs_unnested$pert_idose==-666)] <- '10 ÂµM'
cmap_drugs_unnested$pert_idose[which(cmap_drugs_unnested$pert_idose=='0.1 %')] <- '10 ÂµM'
print(unique(cmap_drugs_unnested$pert_idose))
cmap_drugs_unnested <- cmap_drugs_unnested %>% separate(pert_idose,c('dose','unit'),sep = " ") %>% mutate(dose=as.numeric(dose))

# Visualize how many drug-target interaction distribution
y <- cmap_drugs_unnested %>% group_by(canonical_smiles) %>% summarise(no_targets=n_distinct(Target)) 
png('targets_distribution.png',width = 9,height = 6,units = 'in',res = 600)
ggplot(y,aes(x=no_targets)) + geom_bar() + xlab('Number of targets') + ylab('Number of drugs')+
  ggtitle('Drugs with a specific number of known targets')+
  theme_minimal(base_family = "serif",base_size = 17)+
  theme(plot.title = element_text(hjust = 0.5,size=20))
dev.off()


min_val <- 3 # uM
cmap_drugs_unnested <- cmap_drugs_unnested %>% mutate(log10Dose = log10(dose+1))
hist(cmap_drugs_unnested$log10Dose,breaks=20)

conditions <- cmap_drugs_unnested %>% select(sig_id,canonical_smiles,log10Dose) %>% unique()
conditions <- conditions %>% spread(canonical_smiles,log10Dose)
conditions <- as.matrix(conditions %>% column_to_rownames('sig_id'))
conditions[which(is.na(conditions))] <- 0.0
write.table(conditions,'preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

smiles = colnames(conditions)
smiles = data.frame(smiles=smiles)
data.table::fwrite(smiles,'lvl3_smiles.csv')

#### smiles to left out for testing------------------
# Load ecfp4 similarities
conditions <- data.table::fread('preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv',sep = "\t") %>% column_to_rownames('V1')
ecfp4 <- data.table::fread('preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv',header = T) %>% column_to_rownames('V1')
png('../article_supplementary_info/chemical_similarity_all_distribution.png',width = 9,height = 6,units = 'in',res = 600)
hist(as.matrix(ecfp4[lower.tri(ecfp4,diag=F)]),
     main = 'Pairwise chemical similarity between all drugs',
     xlab = 'Chemical similarity',
     ylab = 'counts',
     freq = T,
     xlim = c(0,1),
     breaks = 50)
abline(v=0.6,col="red",lwd=2,lty='dashed')
text(0.75, 400, 'Similar drugs space')
text(0.45, 400, 'Dissimilar drugs space')
dev.off()
hist(ecfp4[lower.tri(ecfp4,diag=F)])
diag(ecfp4) <- 0
max_ecfp4 <- apply(ecfp4,2,max)
test_candidates <- names(max_ecfp4)[which(max_ecfp4<0.6)]
test_smiles <- sample(test_candidates,ceiling(length(max_ecfp4)*0.2))
test_smiles <- test_smiles[which(test_smiles!='CS(C)=O')]
#write.csv(test_smiles,'preprocessed_data/ChemicalSims/allcells_test_smiles.csv')
test_smiles <- read.csv('preprocessed_data/ChemicalSims/allcells_test_smiles.csv') %>% select(-X)
conditions <- as.data.frame(conditions) %>% rownames_to_column('sig_id') %>% gather('smiles','value',-sig_id) %>% filter(value!=0) %>% 
  select(-value) %>% unique()
test_conds <- left_join(test_smiles,conditions,by=c('x'='smiles'))
colnames(test_conds) <- c('smile','sig_id')
test_candidates <- test_conds %>% unique()
#write.csv(test_conds,'preprocessed_data/TrainingValidationData/smiles_and_sigs_to_drop.csv')

# Visualize chemical similarity of test and train smiles
ecfp4_test_train <- ecfp4[which(rownames(ecfp4) %in% test_conds$smile),
                          which(!(rownames(ecfp4) %in% test_conds$smile))]
png('../article_supplementary_info/chemical_similarity_test_train_distribution.png',width = 9,height = 6,units = 'in',res = 600)
hist(as.matrix(ecfp4_test_train),
     main = 'Chemical similarity between test drugs and train drugs',
     xlab = 'Chemical similarity',
     ylab = 'counts',
     freq = T,
     xlim = c(0,1),
     breaks = 50)
abline(v=0.6,col="red",lwd=2,lty='dashed')
text(0.75, 400, 'Similar drugs space')
text(0.45, 400, 'Dissimilar drugs space')
dev.off()

max_ecfp4_test_train <- apply(ecfp4_test_train,1,max)
png('../article_supplementary_info/max_chemical_similarity_test_train_distribution.png',width = 9,height = 6,units = 'in',res = 600)
hist(as.matrix(max_ecfp4_test_train),
     main = 'Maximum Chemical similarity of a test drug with some train drug',
     xlab = 'Max Chemical similarity',
     ylab = 'counts',
     freq = T,
     xlim = c(0,0.6),
     breaks = 10)
dev.off()
