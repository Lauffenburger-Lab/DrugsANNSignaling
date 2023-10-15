library(tidyverse)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggforce)
library(ggsignif)
library(ggstatsplot)
library(ggrepel)
library(cmapR)
library(org.Hs.eg.db)
library(rhdf5)
library(doFuture)
cores <- 16
registerDoFuture()
plan(multisession,workers = cores)
library(doRNG)
library(dorothea)

shRNA_kd_significance <- function(sig,sigInfo,cmap){
  gene <- sigInfo$pert_iname[which(sigInfo$sig_id==sig)]
  cell <- sigInfo$cell_id[which(sigInfo$sig_id==sig)]
  if (gene %in% rownames(cmap)){
    value <- cmap[gene,sig]
    data <- sigInfo %>% filter(cell_id==cell)
    null_cell_distribution <- cmap[gene,data$sig_id]
    pval <- sum(null_cell_distribution<=value)/length(null_cell_distribution)
  }else{
    pval <- NA
  }
  return(pval)
}

####  Load Level 5 KOs ---------------
sig_metrics <- read.delim('../data/GSE92742_Broad_LINCS_sig_metrics.txt')
sig_metrics <- sig_metrics %>% filter(is_exemplar==1)
data_lvl5_shs <- read.delim("../data/GSE92742_Broad_LINCS_sig_info.txt")
data_lvl5_shs <- data_lvl5_shs %>% filter(pert_type=='trt_sh')
data_lvl5_shs <- data_lvl5_shs %>% filter(sig_id %in% sig_metrics$sig_id)
data_lvl5_shs <- data_lvl5_shs %>% filter(cell_id %in% c('A375','A549','MCF7',
                                                         'HEPG2','HT29','HA1E',
                                                         'HCC515','VCAP','PC3'))
sig_metrics <- sig_metrics %>% filter(distil_nsample>=3)
data_lvl5_shs <- data_lvl5_shs %>% filter(sig_id %in% sig_metrics$sig_id)
gene_info <- read.delim(file = "../data/GSE92742_Broad_LINCS_gene_info.txt")
gene_info <- gene_info %>% filter(pr_is_bing==1)
### Load cmap
sigIds <- unique(as.character(unique(data_lvl5_shs$sig_id)))
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))
ds_path <- '../data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dorng% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(gene_info$pr_gene_id)),
                      cid = sigs)
}
cmap_lvl5_shrna <-do.call(cbind,cmap_gctx)

df_annotation = data.frame(gene_id=rownames(cmap_lvl5_shrna))
gene_info$gene_id <- as.character(gene_info$pr_gene_id)
df_annotation <- left_join(df_annotation,gene_info)
rownames(cmap_lvl5_shrna) <- df_annotation$pr_gene_symbol

p.values <- sapply(data_lvl5_shs$sig_id,shRNA_kd_significance,data_lvl5_shs,cmap_lvl5_shrna)
data_lvl5_shs <- data_lvl5_shs[which(p.values<0.05),]
cmap_lvl5_shrna <- cmap_lvl5_shrna[,which(colnames(cmap_lvl5_shrna) %in% data_lvl5_shs$sig_id)]
#saveRDS(data_lvl5_shs,'../results/ExperimentalValidation/sigInfo_shrna_filtered.rds')
# data_lvl5_shs <- readRDS('../results/ExperimentalValidation/sigInfo_shrna_filtered.rds')

####  Load Level 5 Ligands ---------------
sig_metrics <- read.delim('../data/GSE92742_Broad_LINCS_sig_metrics.txt')
sig_metrics <- sig_metrics %>% filter(is_exemplar==1)
sig_metrics <- sig_metrics%>% filter(pert_type=='trt_lig')
sig_metrics <- sig_metrics %>% filter(distil_nsample>=3)
data_lvl5_ligands <- read.delim("../data/GSE92742_Broad_LINCS_sig_info.txt")
data_lvl5_ligands <- data_lvl5_ligands %>% filter(pert_type=='trt_lig')
data_lvl5_ligands <- data_lvl5_ligands %>% filter(sig_id %in% sig_metrics$sig_id)
data_lvl5_ligands <- data_lvl5_ligands %>% filter(pert_time<=12)
data_lvl5_ligands <- data_lvl5_ligands %>% filter(cell_id %in% c('A375','A549','MCF7',
                                                                 'HEPG2','HT29','HA1E',
                                                                 'HCC515','VCAP','PC3'))
### Load cmap
sigIds <- unique(as.character(unique(data_lvl5_ligands$sig_id)))
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))
ds_path <- '../data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dorng% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(gene_info$pr_gene_id)),
                      cid = sigs)
}
cmap_lvl5_ligands <-do.call(cbind,cmap_gctx)
df_annotation = data.frame(gene_id=rownames(cmap_lvl5_ligands))
#gene_info$gene_id <- as.character(gene_info$pr_gene_id)
df_annotation <- left_join(df_annotation,gene_info)
print(all(rownames(cmap_lvl5_ligands)==df_annotation$gene_id))
rownames(cmap_lvl5_ligands) <- df_annotation$pr_gene_symbol
####  Load Level 5 drugs ---------------
sig_metrics <- read.delim('../data/GSE92742_Broad_LINCS_sig_metrics.txt')
sig_metrics <- sig_metrics %>% filter(is_exemplar==1)
sig_metrics <- sig_metrics%>% filter(pert_type=='trt_cp')
data_lvl5 <- read.delim("../data/GSE92742_Broad_LINCS_sig_info.txt")
data_lvl5 <- data_lvl5 %>% filter(pert_type=='trt_cp')
data_lvl5 <- data_lvl5 %>% filter(sig_id %in% sig_metrics$sig_id)
data_lvl5 <- data_lvl5 %>% filter(pert_time<=12)
data_lvl5 <- data_lvl5 %>% filter(cell_id %in% c('A375','A549','MCF7',
                                                 'HEPG2','HT29','HA1E',
                                                 'HCC515','VCAP','PC3'))
### Load cmap
sigIds <- unique(as.character(unique(data_lvl5$sig_id)))
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))
ds_path <- '../data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dorng% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(gene_info$pr_gene_id)),
                      cid = sigs)
}
cmap_lvl5 <-do.call(cbind,cmap_gctx)
df_annotation = data.frame(gene_id=rownames(cmap_lvl5))
#gene_info$gene_id <- as.character(gene_info$pr_gene_id)
df_annotation <- left_join(df_annotation,gene_info)
print(all(rownames(cmap_lvl5)==df_annotation$gene_id))
rownames(cmap_lvl5) <- df_annotation$pr_gene_symbol

####  Load Level 5 controls ---------------
sig_metrics <- read.delim('../data/GSE92742_Broad_LINCS_sig_metrics.txt')
sig_metrics <- sig_metrics %>% filter(is_exemplar==1)
sig_metrics <- sig_metrics %>% filter(distil_nsample>3)
sig_metrics <- sig_metrics %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt')
data_lvl5_ctrl <- read.delim("../data/GSE92742_Broad_LINCS_sig_info.txt")
data_lvl5_ctrl <- data_lvl5_ctrl %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt')
data_lvl5_ctrl <- data_lvl5_ctrl %>% filter(sig_id %in% sig_metrics$sig_id)
data_lvl5_ctrl <- data_lvl5_ctrl %>% filter(pert_time<=12)
data_lvl5_ctrl <- data_lvl5_ctrl %>% filter(cell_id %in% c('A375','A549','MCF7',
                                                           'HEPG2','HT29','HA1E',
                                                           'HCC515','VCAP','PC3'))
### Load cmap
sigIds <- unique(as.character(unique(data_lvl5_ctrl$sig_id)))
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))
ds_path <- '../data/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dorng% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(gene_info$pr_gene_id)),
                      cid = sigs)
}
cmap_lvl5_ctrl <-do.call(cbind,cmap_gctx)
df_annotation = data.frame(gene_id=rownames(cmap_lvl5_ctrl))
#gene_info$gene_id <- as.character(gene_info$pr_gene_id)
df_annotation <- left_join(df_annotation,gene_info)
print(all(rownames(cmap_lvl5_ctrl)==df_annotation$gene_id))
rownames(cmap_lvl5_ctrl) <- df_annotation$pr_gene_symbol

# Combine all data---------
cmap_lvl5 <- cbind(cmap_lvl5,cmap_lvl5_ligands,cmap_lvl5_shrna,cmap_lvl5_ctrl)
cmap_lvl5 <- t(cmap_lvl5)
data_lvl5 <- rbind(data_lvl5,data_lvl5_ligands,data_lvl5_shs,data_lvl5_ctrl)
gc()
# saveRDS(data_lvl5,'../results/ExperimentalValidation/sigInfo_filtered.rds')
# saveRDS(cmap_lvl5,'../results/ExperimentalValidation/cmap_filtered.rds')
#data_lvl5 <- readRDS('../results/ExperimentalValidation/sigInfo_filtered.rds')
#cmap_lvl5 <- readRDS('../results/ExperimentalValidation/cmap_filtered.rds')

### Infer transcription factors with Dorothea----
cmap_lvl5 <- t(cmap_lvl5)
cmap_standardized = cmap_lvl5-rowMeans(cmap_lvl5)
colnames(cmap_standardized) = colnames(cmap_lvl5)
minNrOfGenes = 5

dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]

# Estimate TF activities
settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities = run_viper(cmap_standardized, dorotheaData, options =  settings)

# Transform results
TF_activities <- t(TF_activities)
#hist(1/(1+exp(-TF_activities)),main='Inferred TF activities',xlab='activity')
#TF_activities <- 1/(1+exp(-TF_activities))

### Experiment --------------
cell = 'A375'
tf = 'FOXM1'
node = 'CDK2' #CDK2
drug_sig <- "CPC014_A375_6H:BRD-K23192422-001-01-1:10"
drug_name <- "lestaurtinib"
drug <- "C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13"

# shrnas info
shRNA_sig <- data_lvl5_shs %>% filter(cell_id==cell) %>% 
  filter(pert_iname==node | pert_iname %in% c('ADRB1','EGFR','NTRK1','FLT3','CDK1','CDK3','CDK4','CDK5','CDK6'))
### Select dose that makes more sense
shRNA_sig$pert_dose <- as.numeric(shRNA_sig$pert_dose)
shRNA_sig <- shRNA_sig %>% filter(pert_dose==1.0)

# see ligand info
lig_sig <-data_lvl5_ligands %>% filter(cell_id==cell) %>% 
  filter(pert_iname==node | pert_iname %in% c('ADRB1','EGFR','NTRK1','FLT3','CDK1','CDK3','CDK4','CDK5','CDK7','CDK6'))
### Select dose that makes more sense
lig_sig$pert_dose <- as.numeric(lig_sig$pert_dose)
#lig_sig <- lig_sig %>% filter(pert_dose==1.0)

# see dmso info
dmso_sig <-data_lvl5_ctrl %>% filter(pert_type=='ctl_vehicle') %>% filter(cell_id==cell)
### Select DMSO or PBS that makes more sense
dmso_sig <- dmso_sig %>% filter(pert_iname=='DMSO')

# see untreated info
untreated_sig <-data_lvl5_ctrl %>% filter(pert_type!='ctl_vehicle') %>% filter(cell_id==cell)

# know inhibitor of node if it exists
sig_inhibitor   <- 'CPC017_A375_6H:BRD-K71860425-001-01-1:10'
inhi_name <- 'CDK2-5-inhibitor'

# in-silico validation
insilico_val <- TF_activities[c(drug_sig,lig_sig$sig_id,shRNA_sig$sig_id,dmso_sig$sig_id,untreated_sig$sig_id,sig_inhibitor),tf]
insilico_val <- as.data.frame(insilico_val)
colnames(insilico_val) <- tf
insilico_val <- insilico_val %>% rownames_to_column('sig_id')
insilico_val <- left_join(insilico_val,data_lvl5 %>% dplyr::select(sig_id,pert_iname)) 
insilico_val <- insilico_val %>% mutate(perturbation= ifelse(sig_id %in% lig_sig$sig_id,paste0(pert_iname,' ligand'),
                                                             ifelse(sig_id %in% shRNA_sig$sig_id,paste0(pert_iname,' shRNA'),
                                                                    ifelse(sig_id %in% dmso_sig$sig_id,unique(dmso_sig$pert_iname),
                                                                           ifelse(sig_id %in% untreated_sig$sig_id,'Untreated',
                                                                                  ifelse(sig_id==sig_inhibitor,inhi_name,drug_name))))))
insilico_val_2plot <-insilico_val %>% filter(!grepl('shRNA',perturbation))
insilico_val_2plot$perturbation <- factor(insilico_val_2plot$perturbation,
                                          levels = c('DMSO','Untreated','CDK2-5-inhibitor','EGFR ligand','FLT3 ligand','lestaurtinib'))
ggdotplot(insilico_val_2plot,x='perturbation',y=tf,fill = 'perturbation')+
  #ylim(c(0,1)) +
  ylab(paste0(tf,' TF activity')) +
  geom_hline(yintercept = 0,linetype='dashed',color='black',linewidth=1) 
ggsave('../results/ExperimentalValidation/l1000_ligands_leustartinib_foxm1.png',
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)

## See distribution of FOXM1 activity of CK2 inhibitors----------------------
drug_targets <- readRDS('../preprocessing/preprocessed_data/drug_targets_space.rds')
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]
cdk2_inhibitors <- rownames(drug_targets)[which(drug_targets[,'P24941']==1)][which(rownames(drug_targets)[which(drug_targets[,'P24941']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
cdk2_inhibitors <- names(which(rowSums(drug_targets[cdk2_inhibitors,])<10))
cdk1_inhibitors <- rownames(drug_targets)[which(drug_targets[,'P06493']==1)][which(rownames(drug_targets)[which(drug_targets[,'P06493']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
cdk1_inhibitors <- names(which(rowSums(drug_targets[cdk1_inhibitors,])<10))
cdk6_inhibitors <- rownames(drug_targets)[which(drug_targets[,'Q00534']==1)][which(rownames(drug_targets)[which(drug_targets[,'Q00534']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
cdk6_inhibitors <- names(which(rowSums(drug_targets[cdk6_inhibitors,])<10))
egfr_inhibitors <- rownames(drug_targets)[which(drug_targets[,'P00533']==1)][which(rownames(drug_targets)[which(drug_targets[,'P00533']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
egfr_inhibitors <- names(which(rowSums(drug_targets[egfr_inhibitors,])<10))
flt3_inhibitors <- rownames(drug_targets)[which(drug_targets[,'P36888']==1)][which(rownames(drug_targets)[which(drug_targets[,'P36888']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
flt3_inhibitors <- names(which(rowSums(drug_targets[flt3_inhibitors,])<10))
jak2_inhibitors <- rownames(drug_targets)[which(drug_targets[,'O60674']==1)][which(rownames(drug_targets)[which(drug_targets[,'O60674']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
jak2_inhibitors <- names(which(rowSums(drug_targets[jak2_inhibitors,])<10))
ntrk1_inhibitors <- rownames(drug_targets)[which(drug_targets[,'P04629']==1)][which(rownames(drug_targets)[which(drug_targets[,'P04629']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
ADRB1_inhibitors <- rownames(drug_targets)[which(drug_targets[,'P08588']==1)][which(rownames(drug_targets)[which(drug_targets[,'P08588']==1)]!='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')]
ADRB1_inhibitors <- names(which(rowSums(drug_targets[ADRB1_inhibitors,])<10))

#put all in a data frame
inhibitors <- rbind(data.frame(drug=cdk1_inhibitors,type=rep('CDK1 inhibitors',length(cdk1_inhibitors))),
                    data.frame(drug=cdk2_inhibitors,type=rep('CDK2 inhibitors',length(cdk2_inhibitors))),
                    data.frame(drug=cdk6_inhibitors,type=rep('CDK6 inhibitors',length(cdk6_inhibitors))),
                    data.frame(drug=egfr_inhibitors,type=rep('EGFR inhibitors',length(egfr_inhibitors))),
                    data.frame(drug=flt3_inhibitors,type=rep('FLT3 inhibitors',length(flt3_inhibitors))),
                    data.frame(drug=ADRB1_inhibitors,type=rep('ADRB1 inhibitors',length(ADRB1_inhibitors))),
                    data.frame(drug=ntrk1_inhibitors,type=rep('NTRK1 inhibitors',length(ntrk1_inhibitors))))

conditions <- data.table::fread('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv',sep = "\t") %>% column_to_rownames('V1')
conditions <- conditions %>% rownames_to_column('sample') %>% gather('drug','value',-sample) %>% filter(value>0) %>%
  dplyr::select(-value) %>% unique()
inhibitors <- left_join(inhibitors,conditions)
inhibitors <- inhibitors %>% filter(sample %in% rownames(TF_activities))
inhibitors <- inhibitors %>% filter(grepl('A375',sample))

## Merge with other 
insilico_val_inhibitors  <- as.data.frame(TF_activities) %>% rownames_to_column('sample') %>% filter(sample %in% inhibitors$sample)
insilico_val_inhibitors <- insilico_val_inhibitors %>% dplyr::select(sample,all_of(c(tf)))
insilico_val_inhibitors <- left_join(inhibitors,insilico_val_inhibitors) %>% 
  dplyr::select(c('sig_id'='sample'),FOXM1,c('perturbation'='type'))
insilico_val_2plot <- insilico_val_2plot %>% dplyr::select(sig_id,FOXM1,perturbation)
insilico_val_2plot <- insilico_val_2plot %>% filter(perturbation!='CDK2-5-inhibitor')

# saveRDS(insilico_val_2plot,'../results/ExperimentalValidation/insilico_val_2plot_figure4.rds')
# saveRDS(insilico_val_inhibitors,'../results/ExperimentalValidation/insilico_val_inhibitors_figure4.rds')
# insilico_val_2plot <-  readRDS('../results/ExperimentalValidation/insilico_val_2plot_figure4.rds')
# insilico_val_inhibitors <-  readRDS('../results/ExperimentalValidation/insilico_val_inhibitors_figure4.rds')

ggdotplot(rbind(insilico_val_2plot,insilico_val_inhibitors),x='perturbation',y=tf,fill = 'perturbation',shape='group')+
  ylim(c(-5.5,4))+
  geom_segment(aes(x = 'lestaurtinib', xend = 'lestaurtinib',
                   y = -5.5,yend = -2.6),
               linewidth=1,color='red',arrow = arrow())+
  #ylab(paste0(tf,' activity')) +
  ylab('') + 
  geom_hline(yintercept = 0,linetype='dashed',color='black',linewidth=1) +
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=14,angle=90),
        legend.position = 'none')
ggsave('../figures/figure4B.eps',
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)

# Interesting resutls to do some in-silico network simulations----------------------------------------------------
# Load the simulated in-silico perturbations
insilico_perturbations <- data.table::fread('../results/ExperimentalValidation/inSilicoKOs_minus10.csv')
colnames(insilico_perturbations)[1] <- 'model'
insilico_perturbations <- insilico_perturbations %>% gather('perturbation','FOXM1 activity',-model)
insilico_perturbations <- insilico_perturbations %>% 
  filter(perturbation %in% c('EGFR','NTRK1','FLT3','ADRB1','CDK1','CDK2','lestaurtinib','DMSO'))

# # run the line bellow only when using the reduced net
# insilico_perturbations <- insilico_perturbations %>%
#   mutate(`FOXM1 activity`=ifelse(perturbation %in% c('EGFR','NTRK1','FLT3','ADRB1',
#                                                      'CDK1+EGFR','CDK1+NTRK1','CDK1+FLT3','CDK1+ADRB1'),
#                                  NA,`FOXM1 activity`))

stats.tests = insilico_perturbations %>% filter(!is.na(`FOXM1 activity`)) %>%
  rstatix::wilcox_test(`FOXM1 activity` ~ perturbation) %>% 
  rstatix::adjust_pvalue(method = 'bonferroni') %>% ungroup()

insilico_perturbations$perturbation <- factor(insilico_perturbations$perturbation,
                                              levels = c("lestaurtinib",
                                                         "DMSO",
                                                         #"on-target",
                                                         #"off-target",
                                                         "ADRB1",
                                                         "EGFR",
                                                         "NTRK1",
                                                         "FLT3",
                                                         "CDK1",
                                                         "CDK2"))
mu_dmso <- insilico_perturbations %>% filter(perturbation=='DMSO')
mu_dmso <- median(mu_dmso$`FOXM1 activity`)
mu_lestaurtinib <- insilico_perturbations %>% filter(perturbation=='lestaurtinib')
mu_lestaurtinib <- median(mu_lestaurtinib$`FOXM1 activity`)

# #THIS IS ALSO FOR THE REDUCED SUBNETWORK
# ggboxplot(insilico_perturbations,
#           x='perturbation',y='FOXM1 activity',
#           color='perturbation',
#           add='jitter') +
#   annotate('text',x=c(3,4,5,6),y=0.2,label='N/A')+
#   scale_color_manual(values = c("black","#E76BF3FF",
#                                 "#FF62BCFF","#DE8C00",
#                                 'black', 'black',
#                                 "#F8766D", 'black',
#                                 "#7CAE00","#B79F00"))+
#   #geom_hline(yintercept = 0.5,linetype='dashed',color='black',linewidth=1) +
#   geom_hline(yintercept = mu_dmso,linetype='dashed',color="#DE8C00",linewidth=1) +
#   annotate('text',x=3.5,y=0.6,label='DMSO-induced median activity',color="#DE8C00",size=6)+
#   geom_hline(yintercept = mu_lestaurtinib,linetype='dashed',color="#F8766D",linewidth=1) +
#   annotate('text',x=3.5,y=0.37,label='Lestaurtinib-induced median activity',color="#F8766D",size=6)+
#   theme(text=element_text(size=24),
#         axis.text.x = element_text(angle = 0),
#         legend.position = 'none')+
#   stat_pvalue_manual(stats.tests %>% mutate(x_position = ifelse(group1=='lestaurtinib',group2,group1)) %>%
#                        filter(group1=='lestaurtinib' | group2=='lestaurtinib'),
#                      label = "{p.adj.signif}",
#                      x = "x_position",
#                      y.position = 0.85,
#                      size = 10)

## Comment this if you are analyzing the reuced subnetwork
ggboxplot(insilico_perturbations %>% filter(!grepl('-target',perturbation)),
          x='perturbation',y='FOXM1 activity',
          color='perturbation',
          add='jitter') +
  # geom_hline(yintercept = 0.5,linetype='dashed',color='black',linewidth=1) +
  geom_hline(yintercept = mu_dmso,linetype='dashed',color="#DE8C00",linewidth=1) +
  annotate('text',x=1.1,y=0.83,label='DMSO-induced median activity',color="#DE8C00",size=6)+
  geom_hline(yintercept = mu_lestaurtinib,linetype='dashed',color="#F8766D",linewidth=1) +
  annotate('text',x=2.25,y=0.38,label='Lestaurtinib-induced median activity',color="#F8766D",size=6)+
  theme(text=element_text(size=24),
        axis.text.x = element_text(angle = 0),
        legend.position = 'none') +
  stat_pvalue_manual(stats.tests %>% mutate(x_position = ifelse(group1=='lestaurtinib',group2,group1)) %>%
                       filter(group1=='lestaurtinib' | group2=='lestaurtinib'),
                     label = "{p.adj.signif}",
                     x = "x_position",
                     y.position = 0.85,
                     size = 10)

ggsave('../figures/figure4C.eps',
       device= cairo_ps,
       scale = 1,
       width = 20,
       height = 6,
       units = "in",
       dpi = 600)

### KO level versus activity for knockdown in-silico experiments---------------
insilico_perturbations_all <- data.table::fread('../results/ExperimentalValidation/inSilicoKOs_dose_response.csv')
colnames(insilico_perturbations_all)[1] <- 'model'
insilico_perturbations_all <- insilico_perturbations_all %>% dplyr::select(model,level,ADRB1,EGFR,NTRK1,FLT3,CDK1,CDK2)
insilico_perturbations_all <- insilico_perturbations_all %>% gather('knockdown','FOXM1 activity',-model,-level)
insilico_perturbations_all <- insilico_perturbations_all %>% group_by(knockdown,level) %>% 
  mutate(act_sd = sd(`FOXM1 activity`))%>% 
  mutate(`FOXM1 activity` = mean(`FOXM1 activity`))  %>% 
  ungroup()
insilico_perturbations_all$knockdown <- factor(insilico_perturbations_all$knockdown,
                                               levels = c('ADRB1','EGFR','NTRK1','FLT3','CDK1','CDK2')) 
p1 <- ggplot(insilico_perturbations_all,
             aes(x=level,
                 y=`FOXM1 activity`,
                 color=knockdown)) +
  geom_line(linewidth=1)+
  geom_point() +
  geom_errorbar(aes(ymin = `FOXM1 activity` - act_sd/sqrt(50), 
                    ymax = `FOXM1 activity` + act_sd/sqrt(50)), 
                size = 1 ,width = 3)+
  xlab('KO strength level')+
  scale_color_manual(values = c("#00BA38", "#00C08B",
                                "#00BFC4", "#00B4F0",
                                "#619CFF", "#C77CFF"))+
  geom_hline(yintercept = 0.5,linetype='dashed',color='black',linewidth=1) +
  geom_hline(yintercept = mu_dmso,linetype='dashed',color="#DE8C00",linewidth=1) +
  annotate('text',x=100,y=0.76,label='DMSO-induced median activity',color="#DE8C00",size=5)+
  geom_hline(yintercept = mu_lestaurtinib,linetype='dashed',color="#F8766D",linewidth=1) +
  annotate('text',x=100,y=0.44,label='Lestaurtinib-induced median activity',color="#F8766D",size=5)+
  ggtitle('FOXM1 activity as a function of the knockdown strength')+
  theme_minimal() +
  theme(text=element_text(size=19),
        plot.title = element_text(hjust=0.5))

# Do the same to check that the node has been knocked successfully
insilico_nodes_ko_all <- data.table::fread('../results/ExperimentalValidation/inSilicoKOs_node_response.csv')
colnames(insilico_nodes_ko_all)[1] <- 'model'
insilico_nodes_ko_all <- insilico_nodes_ko_all %>% dplyr::select(model,level,ADRB1,EGFR,NTRK1,FLT3,CDK1,CDK2)
insilico_nodes_ko_all <- insilico_nodes_ko_all %>% gather('knockdown','activity',-model,-level)
insilico_nodes_ko_all <- insilico_nodes_ko_all %>% group_by(knockdown,level) %>% 
  mutate(act_sd = sd(activity))%>% 
  mutate(activity = mean(activity))  %>% 
  ungroup()
insilico_nodes_ko_all$knockdown <- factor(insilico_nodes_ko_all$knockdown,
                                          levels = c('ADRB1','EGFR','NTRK1','FLT3','CDK1','CDK2')) 
p2 <- ggplot(insilico_nodes_ko_all %>% filter(level<=20),
             aes(x=level,
                 y=activity,
                 color=knockdown)) +
  geom_line(linewidth=1)+
  geom_point() +
  geom_errorbar(aes(ymin = activity - act_sd/sqrt(50), 
                    ymax = activity + act_sd/sqrt(50)), 
                size = 1 ,width = 0.5)+
  xlab('KO strength level')+
  scale_color_manual(values = c("#00BA38", "#00C08B",
                                "#00BFC4", "#00B4F0",
                                "#619CFF", "#C77CFF"))+
  geom_hline(yintercept = 0.5,linetype='dashed',color='black',linewidth=1) +
  ggtitle('Signaling node activity after knockdown experiments of given strength')+
  theme_minimal() +
  theme(text=element_text(size=19),
        plot.title = element_text(hjust=0.5))

# Combine into one plot
p1+p2
ggsave('../article_supplementary_info/supple_figure11CD.eps',
       device = cairo_ps,
       scale = 1,
       width = 18,
       height = 9,
       units = "in",
       dpi = 600)
