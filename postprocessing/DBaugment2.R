library(tidyverse)
library("dorothea")
library(cmapR)

## Load all cmap_drugs and cmap+DrugBank drugs and targets---------------------
## And for cmap augment drug-targets with tarhets from drugbank
# #cmap_drugs <- readRDS('all_cmapdrugs_moa+target_v1.rds')
drug_targets_all <- readRDS('drug_targets_space.rds')
drug_targets_all <- distinct(as.data.frame(drug_targets_all) %>% rownames_to_column('canonical_smiles') %>% 
                               gather('Target','value',-canonical_smiles) %>%
                               filter(value!=0) %>% select(-value))
#print(length(which(cmap_drugs_unt$canonical_smiles %in% rownames(drug_targets_all))))
data <- readRDS("all_cmap_sigs_with_pert_info.rds")
dorothea <- read.delim('l1000_allgenes_lvl3_tfs.tsv',row.names = 1)

# Get specidically TPCA-1 experiment in HEPG2 to add it in the data
minNrOfGenes = 5
dorotheaData = read.table('TF activities/annotation/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]
# Get cmap data only for that experiment
tpca1_data <- data %>% filter(pert_iname=='TPCA-1') %>% 
  filter(cell_id %in% c('A375','A549','HA1E','HCC515','HT29','PC3','MCF7','HEPG2','VCAP'))
tpca1_data <- tpca1_data %>% mutate(inst_id=strsplit(distil_id,"\\|")) %>% unnest(inst_id)
# Path to raw data
ds_path <- '../GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
gene_info <- read.delim(file = "../deepSIBA/preprocessing/data_preprocessing/cmap/GSE92742_Broad_LINCS_gene_info.txt")
landmark_df <- read.delim(file = "../deepSIBA/preprocessing/data_preprocessing/cmap/cmap_landmark_genes.txt")
gene_info <- gene_info %>% filter(pr_is_bing==1)
gctx_file <- parse_gctx(ds_path ,rid = unique(as.character(gene_info$pr_gene_id)),cid = tpca1_data$inst_id)
E <- gctx_file@mat
E.standardized = E-rowMeans(E)
colnames(E.standardized) = colnames(E)
gene_info <- read.delim(file = "../deepSIBA/preprocessing/data_preprocessing/cmap/GSE92742_Broad_LINCS_gene_info.txt")
landmark_df <- read.delim(file = "../deepSIBA/preprocessing/data_preprocessing/cmap/cmap_landmark_genes.txt")
gene_info <- gene_info %>% filter(pr_is_bing==1)
print(all(rownames(E.standardized)==gene_info$gene_id))
rownames(E.standardized) <- gene_info$pr_gene_symbol
# Estimate TF activities
settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities = run_viper(E.standardized, dorotheaData, options =  settings)
annot <- read.csv('TF activities/annotation/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t')
annot <- annot %>% dplyr::select(Gene.names...primary..,Entry)
tfs <- data.frame("Gene.names...primary.."=rownames(TF_activities))
tfs <- left_join(tfs,annot)
print(all(rownames(TF_activities)==tfs$Gene.names...primary..))
rownames(TF_activities) <- tfs$Entry
TF_activities <- t(TF_activities)
hist(1/(1+exp(-TF_activities)),100,main='Inferred TF activities of TPCA-1 experiments',xlab='activity')
TF_activities <- 1/(1+exp(-TF_activities))
data_replicates <- tpca1_data %>% filter(distil_nsample>=2) %>% select(sig_id,inst_id)
data_replicates <- left_join(data_replicates,as.data.frame(TF_activities) %>% rownames_to_column('inst_id'))
TF_activities <- data_replicates %>% select(-inst_id)
ind <- which(colnames(TF_activities)=='sig_id')
TF_activities <- aggregate(TF_activities[,2:116],by=list(TF_activities$sig_id),FUN=median)
TF_activities <- TF_activities %>% column_to_rownames('Group.1')
TF_activities <- as.matrix(TF_activities)
hist(TF_activities)
print(all(colnames(dorothea) %in% colnames(TF_activities)))
dorothea <- rbind(dorothea,TF_activities)
####

data <- data %>% filter(sig_id %in% rownames(dorothea))
data <- data %>% mutate(canonical_smiles = ifelse(pert_id=='DMSO','CS(C)=O',
                                     ifelse(pert_id=='H2O','O',
                                            ifelse(pert_id=='PBS','Ctrl_PBS',
                                                   ifelse(pert_type=='ctl_untrt','Ctrl_UnTrt',canonical_smiles)))))
data_untreated <- data %>% filter(canonical_smiles=='O' | canonical_smiles=='PBS' | canonical_smiles=='Ctrl_UnTrt')
data_untreated <- data_untreated %>% mutate(Target = 'None')
cmap_drugs_unnested <- left_join(data,drug_targets_all)
cmap_drugs_unnested <- rbind(cmap_drugs_unnested,data_untreated)
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(!is.na(Target)) %>% filter(Target!='') %>% filter(Target!=' ')
dorothea <- dorothea[which(rownames(dorothea) %in% cmap_drugs_unnested$sig_id),]
gc()

## Start trimming of pkn with these new data-------
# Before trimming designate interactions to mannually add (again) after trimming #
pkn <- data.table::fread('Network Construction/model/pkn.tsv')
annot <- read.csv('TF activities/annotation/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t')
annot <- annot %>% select(Gene.names...primary..,Entry)
colnames(annot)[1] <- 'Target'
print(length(which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source))))) # 5 TFs not in PKN for lands, 9 for all
write.table(dorothea, file = 'Model/data/Trimmed_l1000_allgenes_lvl3_tfs_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
tfs_targeted <- unique(cmap_drugs_unnested$Target[which(cmap_drugs_unnested$Target %in% colnames(dorothea))])
write.table(data.frame('Entry'=tfs_targeted), file = 'Model/data/tfs_targetd_alls_genes_lvl3_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
pkn_tfs_targeted <- pkn %>% filter((source %in% tfs_targeted) | (target %in% tfs_targeted))
write.table(pkn_tfs_targeted, file = 'Model/data/L1000_latest_Add_lvl3_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
write.table(distinct(cmap_drugs_unnested %>% select('source'='canonical_smiles','target'='Target') %>% filter(target!='None')), 
            file = 'Network Construction/model/L1000_lvl3_DT_drugbank.tsv', quote=FALSE, sep = "\t", row.names = F)

# After trimming net #
pkn <- data.table::fread(('Network Construction/l1000_lvl3_withsignor-Model_drugbank.tsv'))
nodes <- unique(c(pkn$source,pkn$target))
cmap_drugs <- cmap_drugs_unnested %>% 
  filter((Target %in% nodes) | (Target=='None') | (canonical_smiles=='CS(C)=O.CS(C)=O') |
           (canonical_smiles=='CS(C)=O') | (canonical_smiles=='O'))
print(length(unique(cmap_drugs$canonical_smiles))) # 969
print(length(unique(cmap_drugs$Target))) # 446
saveRDS(cmap_drugs,'TrimmedFinal_lvl3_allgenes_all_conditions_drugbank.rds')
dorothea <- read.delim('Model/data/Trimmed_l1000_allgenes_lvl3_tfs_drugbank.tsv',row.names = 1)
print(length(which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source))))) # 101 out of 106 for alls
dorothea <- dorothea[,which(colnames(dorothea) %in% unique(c(pkn$target,pkn$source)))]
write.table(dorothea, file = 'Model/data/TrimmedFinal_l1000_allgenes_lvl3_tfs_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### Split data according to cells ----------------------
pkn <- data.table::fread(('Network Construction/l1000_lvl3_withsignor-Model_drugbank.tsv'))
dorothea <- read.delim('Model/data/TrimmedFinal_l1000_allgenes_lvl3_tfs_drugbank.tsv',row.names = 1)
data <- data %>% filter(sig_id %in% rownames(dorothea))
cmap_drugs_unnested <- cmap_drugs %>% filter(sig_id %in% rownames(dorothea))

### Make conditions ###----------------------------
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
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter((Target %in% c(pkn$source,pkn$target)) | (Target=='None')) %>% unique()
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(canonical_smiles %in% common)
saveRDS(cmap_drugs_unnested,'TrimmedFinal_lvl3_conditions_for_pairedcells_drugbank.rds')

#### Split sets------------------------
cmap_drugs_unnested <- readRDS('TrimmedFinal_lvl3_conditions_for_pairedcells_drugbank.rds')
uniq_cells <- unique(cmap_drugs_unnested$cell_id)

### Make drug target matrix and cell matrix ###
uniq_targets <- unique(cmap_drugs_unnested$Target)

binary_classes <- function(sig,df,targets,type=c('smiles','sig_id'),tar=c('Target','cell_id')){
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
write.table(x,'Model/data/L1000_lvl3_allcells-drugs_targets_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

xc <- t(sapply(unique(cmap_drugs_unnested$sig_id),binary_classes,cmap_drugs_unnested,uniq_cells,type='sig_id',tar='cell_id'))
write.table(xc,'Model/data/L1000_lvl3-conditions_cells_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

print(unique(cmap_drugs_unnested$pert_idose))
cmap_drugs_unnested$pert_idose[which(cmap_drugs_unnested$pert_idose==-666)] <- '10 ÂµM'
cmap_drugs_unnested$pert_idose[which(cmap_drugs_unnested$pert_idose=='0.1 %')] <- '10 ÂµM'
print(unique(cmap_drugs_unnested$pert_idose))
cmap_drugs_unnested <- cmap_drugs_unnested %>% separate(pert_idose,c('dose','unit'),sep = " ") %>% mutate(dose=as.numeric(dose))

min_val <- 3 # uM
#cmap_drugs_unnested <- cmap_drugs_unnested %>% mutate(log10Dose = ifelse(unit=="nM",log10(1000*dose/min_val+1),log10(dose/min_val+1)))
cmap_drugs_unnested <- cmap_drugs_unnested %>% mutate(log10Dose = log10(dose+1))
hist(cmap_drugs_unnested$log10Dose,breaks=20)

conditions <- cmap_drugs_unnested %>% select(sig_id,canonical_smiles,log10Dose) %>% unique()
conditions <- conditions %>% spread(canonical_smiles,log10Dose)
conditions <- as.matrix(conditions %>% column_to_rownames('sig_id'))
conditions[which(is.na(conditions))] <- 0.0
write.table(conditions,'Model/data/L1000_lvl3_allcells-conditions_drugs_drugbank.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

smiles = colnames(conditions)
smiles = data.frame(smiles=smiles)
data.table::fwrite(smiles,'lvl3_smiles_drugbank.csv')

#### smiles to left out for testing------------------
# Load ecfp4 similarities
conditions <- data.table::fread('Model/data/L1000_lvl3_allcells-conditions_drugs_drugbank.tsv',sep = "\t") %>% column_to_rownames('V1')
ecfp4 <- data.table::fread('out_lvl3_similaritiess_drugbank.csv',header = T) %>% column_to_rownames('V1')
png('chemical_similarity_all_distribution_drugbank.png',width = 9,height = 6,units = 'in',res = 600)
hist(as.matrix(ecfp4[lower.tri(ecfp4,diag=F)]),
     main = 'Pairwise chemical similarity between all drugs',
     xlab = 'Chemical similarity',
     ylab = 'counts',
     freq = T,
     xlim = c(0,1),
     breaks = 50)
# hist(as.matrix(ecfp4[lower.tri(ecfp4,diag=F)]),
#      breaks = 500,
#      freq = F,
#      col=rgb(0.8,0.8,0.8,0.5),
#      border = rgb(0.8,0.8,0.8,0.1),
#      add=TRUE)
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
test_smiles <- test_smiles[which(test_smiles!='O')]
test_smiles <- test_smiles[which(test_smiles!='PBS')]
test_smiles <- test_smiles[which(test_smiles!='Ctrl_UnTrt')]
#write.csv(test_smiles,'allcells_test_smiles_drugbank.csv')
test_smiles <- read.csv('allcells_test_smiles_drugbank.csv') %>% select(-X)
conditions <- as.data.frame(conditions) %>% rownames_to_column('sig_id') %>% gather('smiles','value',-sig_id) %>% filter(value!=0) %>% 
  select(-value) %>% unique()
test_conds <- left_join(test_smiles,conditions,by=c('x'='smiles'))
colnames(test_conds) <- c('smile','sig_id')
#write.csv(test_conds,'smiles_and_sigs_to_drop_drugbank.csv')

# Visualize chemical similarity of test and train smiles
ecfp4_test_train <- ecfp4[which(rownames(ecfp4) %in% test_conds$smile),
                          which(!(rownames(ecfp4) %in% test_conds$smile))]
png('chemical_similarity_test_train_distribution_drugbank.png',width = 9,height = 6,units = 'in',res = 600)
hist(as.matrix(ecfp4_test_train),
     main = 'Chemical similarity between test drugs and train drugs',
     xlab = 'Chemical similarity',
     ylab = 'counts',
     freq = T,
     xlim = c(0,1),
     breaks = 50)
# hist(as.matrix(ecfp4[lower.tri(ecfp4,diag=F)]),
#      breaks = 500,
#      freq = F,
#      col=rgb(0.8,0.8,0.8,0.5),
#      border = rgb(0.8,0.8,0.8,0.1),
#      add=TRUE)
abline(v=0.6,col="red",lwd=2,lty='dashed')
text(0.75, 400, 'Similar drugs space')
text(0.45, 400, 'Dissimilar drugs space')
dev.off()

max_ecfp4_test_train <- apply(ecfp4_test_train,1,max)
png('max_chemical_similarity_test_train_distribution_drugbank.png',width = 9,height = 6,units = 'in',res = 600)
hist(as.matrix(max_ecfp4_test_train),
     main = 'Maximum Chemical similarity of a test drug with some train drug',
     xlab = 'Max Chemical similarity',
     ylab = 'counts',
     freq = T,
     xlim = c(0,0.6),
     breaks = 10)
dev.off()
