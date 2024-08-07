# Load requeired packages
library("dorothea")
library(tidyverse)
library(doFuture)
# parallel: set number of workers
cores <- 16
registerDoFuture()
plan(multisession,workers = cores)

# Access command-line arguments
args <- commandArgs(trailingOnly = TRUE)

inputFile <- args[1]
outputFile = args[2]
# inputFile <- 'preprocessed_data/unmerged_l1000_allgenes_lvl3_tfs.rds'
# outputFile = "preprocessed_data/TF_activities/l1000_allgenes_lvl3_tfs.tsv"
# Read inferred TF activities
TF_activities <- readRDS(inputFile)

### First filter out TFS with consistently high variance across replicates in the samples 
### Then filter samples with low correlation between replicates.---------------
data <- readRDS("preprocessed_data/all_cmap_sigs_with_pert_info.rds")
data <- data %>% filter(is_exemplar==1) %>% #filter(tas>0.1) %>%
  filter(pert_itime %in% c("2 h",  "3 h", "4 h","6 h",'9 h','12 h'))
data <- data %>% dplyr::select(sig_id,distil_id,pert_type,pert_idose,pert_itime,cell_id,
                               pert_id,canonical_smiles,distil_nsample) %>% unique()
data <- data %>% mutate(inst_id=strsplit(distil_id,"\\|")) %>% unnest(inst_id)
gc()
data_lvl3 <- read.delim("../data/GSE92742_Broad_LINCS_inst_info.txt")
data_lvl3 <- data_lvl3 %>% dplyr::select(inst_id,pert_type) %>%
  filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt' | pert_type=='trt_cp') %>% dplyr::select(-pert_type) %>%
  unique()
data_merged <- left_join(data,data_lvl3,by="inst_id")
data_merged <- data_merged %>% dplyr::select(sig_id,inst_id,distil_nsample) %>% unique()

data_replicates <- data_merged %>% filter(distil_nsample>=2) %>% select(-distil_nsample)
data_replicates <- left_join(data_replicates,as.data.frame(TF_activities) %>% rownames_to_column('inst_id'))

### Filter samples----------------
replicatesCor <- function(data,s){
  x <- data %>% filter(sig_id==s)
  x <- t(x[,3:ncol(data)])
  c <- cor(x)
  c <- c[lower.tri(c, diag = FALSE)]
  return(mean(c))
}
sigs <- unique(data_replicates$sig_id)

TFs_cor <- foreach(sig = sigs) %dopar% {
  replicatesCor(data_replicates ,
                s = sig)
}

TFs_cor <-do.call(rbind,TFs_cor)
TFs_cor <- as.data.frame(TFs_cor)
TFs_cor <- TFs_cor %>% mutate(sig_id=sigs)
TFs_cor <- TFs_cor %>% select(sig_id,c('corr'='V1'))
TFs_cor <- left_join(TFs_cor,data_merged %>% select(sig_id,distil_nsample) %>% unique())

rep_size <- data_merged %>% filter(distil_nsample>=2)
rep_size <- unique(rep_size$distil_nsample)
rep_size <- rep_size[order(rep_size)]
random_corrs <- matrix(0,nrow=length(rep_size),ncol=1000)
for (i in 1:length(rep_size)){
  n <- rep_size[i]
  for (j in 1:1000){
    tfs <- data_replicates[sample(1:nrow(data_replicates),n),3:ncol(data_replicates)]
    c <- cor(t(tfs))
    c <- c[lower.tri(c, diag = FALSE)]
    random_corrs[i,j] <- mean(c)
  }
  print(paste0('Finished n: ',n))
}
rownames(random_corrs) <- rep_size
# saveRDS(random_corrs,'random_correlations.rds')


TFs_cor$p.value <- 1
for (i in 1:nrow(TFs_cor)){
  TFs_cor$p.value[i] <- sum(random_corrs[which(rownames(random_corrs)==TFs_cor$distil_nsample[i]),]>=TFs_cor$corr[i])/1000
  if (i %% 1000 ==0){
    print(paste0('Finished ',i))
  }
}
# saveRDS(TFs_cor,'preprocessed_data/TF_activities/TFs_cor_repls.rds')

TFs_cor <- TFs_cor %>% filter(p.value<0.05)
### Merge replicates------------------
data_merged <- data_merged %>% filter(sig_id %in% TFs_cor$sig_id | distil_nsample==1) %>% select(-distil_nsample) 

TF_activities <- as.data.frame(TF_activities) %>% rownames_to_column('inst_id')
TF_activities <- TF_activities %>% filter(inst_id %in% data_merged$inst_id)
TF_activities <- left_join(TF_activities,data_merged)
TF_activities <- TF_activities %>% select(-inst_id)
ind <- which(colnames(TF_activities)=='sig_id')
TF_activities <- aggregate(TF_activities[,1:115],by=list(TF_activities$sig_id),FUN=median)
TF_activities <- TF_activities %>% column_to_rownames('Group.1')
TF_activities <- as.matrix(TF_activities)
hist(TF_activities)
write.table(TF_activities, file = outputFile, quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)