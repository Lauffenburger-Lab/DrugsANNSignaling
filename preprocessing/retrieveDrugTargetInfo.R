library(tidyverse)
library(cmapR)
library(org.Hs.eg.db)
library(rhdf5)
library(doFuture)

####  Process level 3 data ---------------
sig_metrics <- read.delim('../data/GSE92742_Broad_LINCS_sig_metrics.txt')
data_lvl5 <- read.delim("../data/GSE92742_Broad_LINCS_sig_info.txt")
data_lvl5 <- data_lvl5 %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt' | pert_type=='trt_cp')
gc()

cells <- unique(data_lvl5$cell_id)
data <- data.frame()
for (cell in cells){
  drug_sigs_per_line(cell,data_lvl5,sig_metrics)
  data <- rbind(data,drug_sigs_per_line(cell,data_lvl5,sig_metrics))
}
data <- data %>% select(-quality)
ctrl_lvl5 <- data_lvl5 %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt')
ctrl_lvl5 <- left_join(ctrl_lvl5,sig_metrics)
ctrl_lvl5 <- ctrl_lvl5 %>% unique()
data <- rbind(data,ctrl_lvl5)
saveRDS(data,"preprocessed_data/sig_info_5_with_exempler.rds")


### read the initial signatures of cmap version 2018----------------------
data <- readRDS("preprocessed_data/sig_info_5_with_exempler.rds")

### read the GSE of the perturbations
pert_info <- read.delim(file = "../data/GSE92742_Broad_LINCS_pert_info.txt")
pert_info <- pert_info %>% filter(pert_type == "trt_cp")
pert_info <- pert_info %>% select(pert_id,inchi_key_prefix,inchi_key,canonical_smiles,pubchem_cid,pert_iname)
tt <- left_join(data,pert_info)
saveRDS(tt,'preprocessed_data/all_cmap_sigs_with_pert_info.rds')
cmap_drugs <- tt %>%  filter(pert_type == "trt_cp") 

### drugs in cmap
cmap_drugs <- cmap_drugs %>% select(pert_id,pert_iname) %>% unique()
### add drug info to drugs in cmap
cmap_drugs <- left_join(cmap_drugs,pert_info)
### add moa and target from the broad repo
broad_repo <- read.delim(file = "../data/Repurposing_Hub_export (1).txt" ,skip = 0)

broad_repo <- broad_repo %>%
  mutate(InChIKey = str_split(InChIKey,pattern = ",")) %>% unnest(InChIKey) %>%
  mutate(InChIKey = str_trim(InChIKey)) %>%
  filter(InChIKey != "") %>% unique()

broad_repo <- broad_repo  %>%
  mutate(broad_id = substr(x = as.character(Id),start = 1,stop = 13)) %>%
  mutate(broad_id_old = substr(x = as.character(Deprecated.ID),start = 1,stop = 13)) %>%
  mutate(Name = toupper(Name))


### lists of indices for mapping between cmap_drugs + broad_repo

new <- which(unique(cmap_drugs$pert_id) %in% broad_repo$broad_id)
old <- which(unique(cmap_drugs$pert_id) %in% broad_repo$broad_id_old)
inchi <- which(cmap_drugs$inchi_key %in% broad_repo$InChIKey)
smiles <- which(cmap_drugs$canonical_smiles %in% broad_repo$SMILES)
names <- which(cmap_drugs$pert_iname %in% broad_repo$Name)

a <- (unique(union(new,inchi)))
b <- unique(union(a,old))
c <- unique(union(b,smiles))
d <- unique(union(c,names))

### 2245 were mapped to broad repo

df1 <- left_join(cmap_drugs,broad_repo,by = c("pert_id"="broad_id"))
df2 <- left_join(cmap_drugs,broad_repo,by = c("pert_id"="broad_id_old"))
df3 <- left_join(cmap_drugs,broad_repo,by = c("inchi_key"="InChIKey"))
df4 <- left_join(cmap_drugs,broad_repo,by = c("canonical_smiles"="SMILES"))
df5 <- left_join(cmap_drugs,broad_repo,by = c("pert_iname"="Name"))

df <- bind_rows(df1,df2,df3,df4,df5) %>% unique()

df_new <- df %>%
  select(pert_id,Target,MOA,Phase,Disease.Area) %>%
  filter(!(is.na(Target)&is.na(MOA))) %>%
  mutate(nchar_target = nchar(as.character(Target))) %>%
  group_by(pert_id) %>% 
  filter(nchar_target == max(nchar_target)) %>% 
  ungroup() %>% 
  unique()


### add the moa and target to cmap_drugs
cmap_drugs <- left_join(cmap_drugs,df_new)
## add dmso and ctrls
cmap_ctrls <- tt%>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt') %>% 
  select(pert_id,pert_iname,inchi_key_prefix,inchi_key,canonical_smiles,pubchem_cid)%>%unique()
cmap_ctrls <- cmap_ctrls %>% mutate(canonical_smiles = ifelse(pert_id=='DMSO','CS(C)=O',
                                                              ifelse(pert_id=='H2O','O',
                                                                     ifelse(pert_id=='PBS','Ctrl_PBS','Ctrl_UnTrt'))))
cmap_ctrls <- cmap_ctrls %>% mutate(Target=ifelse(pert_id=='DMSO','IL5RA, MUC16, MYC','None'))
cmap_ctrls <- cmap_ctrls %>% mutate(MOA=NA,Phase=NA,Disease.Area='control') %>% mutate(nchar_target = nchar(as.character(Target)))
cmap_ctrls <- cmap_ctrls %>% select(colnames(cmap_drugs))
cmap_drugs <- rbind(cmap_drugs,cmap_ctrls)
### save the results
saveRDS(cmap_drugs,"preprocessed_data/all_cmapdrugs_moa+target_v1.rds")

cmap_drugs <- cmap_drugs %>% filter(!is.na(canonical_smiles))
cmap_drugs <- cmap_drugs %>% filter(canonical_smiles!='')
cmap_drugs <- cmap_drugs %>% filter(canonical_smiles!=' ')
cmap_drugs <- cmap_drugs %>% filter(canonical_smiles!='bad')
cmap_drugs <- cmap_drugs %>% filter(!is.na(Target))
cmap_drugs <- cmap_drugs %>% filter(Target!='')
cmap_drugs <- cmap_drugs %>% filter(Target!=' ')
cmap_drugs <- cmap_drugs %>%filter(pubchem_cid!=-666)
cmap_drugs <- cmap_drugs %>% select(canonical_smiles,pubchem_cid,Target,MOA,Phase,Disease.Area)
cmap_drugs <- cmap_drugs %>% group_by(canonical_smiles) %>% mutate(dupl = n()) %>% ungroup()
cmap_drugs1 <- cmap_drugs %>% filter(dupl==1)
cmap_drugs2 <- cmap_drugs %>% filter(dupl>1)
smis <- unique(cmap_drugs2$canonical_smiles)
for (i in 1:length(smis)){
  smi <- smis[i]
  inds <- which(cmap_drugs2$canonical_smiles==smi)
  pubchem <- paste(cmap_drugs2$pubchem_cid[inds],collapse = ',')
  cmap_drugs2$pubchem_cid[inds] <- pubchem
}
cmap_drugs2 <-cmap_drugs2 %>% unique()
cmap_drugs <- rbind(cmap_drugs1,cmap_drugs2) %>% select(-dupl) %>% unique()
saveRDS(cmap_drugs,"preprocessed_data/l1000_drugs_with_targets_all.rds") # 1223 unique drugs in different doses , times, cells, and qualities
cmap_drugs <- left_join(cmap_drugs,tt,by="canonical_smiles") %>% filter(is_exemplar==1) %>%
  select(canonical_smiles,'pubchem_cid'='pubchem_cid.x',Target,MOA,Phase,"Disease.Area")
cmap_drugs <- cmap_drugs %>% unique()
saveRDS(cmap_drugs,"preprocessed_data/l1000_drugs_with_targets_exemplar.rds") # 1214 unique drugs in different doses , times, cells

### Remove drugs with no target in the network at all ###
net <- read.delim('preprocessed_data/PKN/pkn.tsv')
df_net <- data.frame(Entry=unique(c(net$target,net$source)))
annot <- read.delim('../data/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab')
annot <- annot %>% select(Entry,c('Target'=Gene.names...primary..)) %>% unique()
df_net <- left_join(df_net,annot)

cmap_drugs_unnested <-  cmap_drugs %>% mutate(Target=strsplit(Target,", ")) %>% unnest(Target)
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(Target %in% df_net$Target)
print(length(unique(cmap_drugs_unnested$canonical_smiles))) #1015 drugs in the end
saveRDS(cmap_drugs_unnested,'preprocessed_data/l1000_drugs_with_targets_unnested_exemplar.rds')

### Get gene expression ###
cmap_drugs <- left_join(cmap_drugs %>% unique(),tt) %>% unique()
cmap_drugs <- cmap_drugs %>% filter(is_exemplar==1) %>% filter(!is.na(sig_id)) %>% unique()

print(nrow(cmap_drugs %>% filter(pert_itime=='6 h') %>% dplyr::select(canonical_smiles) %>% unique()))
### We have 768 unique drugs in 6 hours with different dose and different cell-line ###

## Read transcriptomics data------------------

gene_info <- read.delim(file = "../data/GSE92742_Broad_LINCS_gene_info.txt")
gene_info <- gene_info %>% filter(pr_is_bing==1)
data <- data %>% filter(is_exemplar==1) %>% 
  filter(pert_itime %in% c("2 h",  "3 h", "4 h","6 h",'9 h','12 h'))
data <- data %>% dplyr::select(sig_id,distil_id,pert_type,pert_idose,pert_itime,cell_id,pert_id) %>% unique()
data <- data %>% mutate(inst_id=strsplit(distil_id,"\\|")) %>% unnest(inst_id)
gc()
data_lvl3 <- read.delim("../data/GSE92742_Broad_LINCS_inst_info.txt")
data_lvl3 <- data_lvl3 %>% dplyr::select(inst_id,pert_type) %>%
  filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt' | pert_type=='trt_cp') %>% dplyr::select(-pert_type) %>%
  unique()

data_merged <- left_join(data,data_lvl3,by="inst_id")
# parallel: set number of workers
cores <- 16
registerDoFuture()
plan(multisession,workers = cores)

# Split sig_ids to run in parallel
sigIds <- unique(as.character(unique(data_merged$inst_id)))
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))

# Parallelize parse_gctx function

# Path to raw data
ds_path <- '../data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'

# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}

# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dopar% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(gene_info$pr_gene_id)),
                      cid = sigs)
}
gex_all <-do.call(cbind,cmap_gctx)
saveRDS(gex_all,'preprocessed_data/l1000_all_genes_lvl3_drugs_with_targets_exemplar.rds')
