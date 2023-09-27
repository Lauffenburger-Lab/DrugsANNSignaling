library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(viridis)
library(patchwork)
library(caret)
library(scales)
library(doFuture)
#options(future.globals.maxSize= 891289600)
#options(future.globals.maxSize= 4831838208)
# parallel set number of workers
cores <- 16
registerDoFuture()
plan(multisession,workers = cores)
library(doRNG)

library(dbparser)
library(XML)

### For speed you can skip the drugbank pre-processing section bellow and load the following files for the rest of the analysis ###

# drug_targets <- readRDS('../preprocessing/preprocessed_data/drug_targets_space.rds')
# df_drugbank <- readRDS('../data/df_drugbank.rds')


### Load drug target data and augment drug target interaction space with drugbank---------
cmap_drugs <- readRDS('../preprocessing/preprocessed_data/l1000_drugs_with_targets_all.rds')
annot <- read.delim('../data/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab')
annot <- annot %>% select(Entry,c('Target'=Gene.names...primary..)) %>% unique()
cmap_drugs <-  cmap_drugs %>% mutate(Target=strsplit(Target,", ")) %>% unnest(Target)
cmap_drugs <- left_join(cmap_drugs,annot)
cmap_drugs <- cmap_drugs %>% filter(!is.na(Entry))
data <- readRDS("../preprocessing/preprocessed_data/all_cmap_sigs_with_pert_info.rds")
data <- data %>% filter(is_exemplar==1)  %>%
  filter(pert_itime %in% c("6 h"))
data <- data %>% dplyr::select(sig_id,pert_iname,pert_type,pert_idose,pert_itime,cell_id,pert_id,canonical_smiles) %>% unique()
cmap_ctrls <- data %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt') %>% unique()
cmap_ctrls <- cmap_ctrls %>% mutate(canonical_smiles = ifelse(pert_id=='DMSO','CS(C)=O',
                                                              ifelse(pert_id=='H2O','O',
                                                                     ifelse(pert_id=='PBS','Ctrl_PBS','Ctrl_UnTrt'))))
cmap_ctrls <- cmap_ctrls %>% mutate(Entry=ifelse(pert_id=='DMSO',c('Q01344', 'Q8WXI7', 'P01106'),'None')) %>% unnest(Entry)
cmap_ctrls <- cmap_ctrls  %>%  mutate(Target=ifelse(Entry=='Q01344','IL5RA',ifelse(Entry=='P01106','MYC',ifelse(Entry=='Q8WXI7','MUC16',Entry)))) 
cmap_ctrls <- cmap_ctrls %>% mutate(MOA=NA,Phase=NA,Disease.Area='control',pubchem_cid=NA)
cmap_drugs <- left_join(cmap_drugs,data)
cmap_ctrls <- cmap_ctrls %>% select(colnames(cmap_drugs))
cmap_drugs <- rbind(cmap_drugs,cmap_ctrls)
cmap_drugs <- cmap_drugs %>% filter(!is.na(Entry))
cmap_drugs <- cmap_drugs %>% filter(Entry!='')
cmap_drugs <- cmap_drugs %>% filter(Entry!=' ')

cmap_drugs_unnested <-  cmap_drugs %>% mutate(Target=strsplit(Target,", ")) %>% unnest(Target)
pkn <- data.table::fread('../preprocessing/preprocessed_data/PKN/pkn.tsv')
annot <- read.delim('../data/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab')
annot <- annot %>% select(Gene.names...primary..,Entry)
colnames(annot)[1] <- 'Target'
cmap_drugs_unnested <- left_join(cmap_drugs_unnested,annot)
cmap_drugs_unnested <- cmap_drugs_unnested %>% mutate(Entry=ifelse(Target=='None','None',Entry))
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(!is.na(Entry))
cmap_drugs_unnested <- cmap_drugs_unnested %>% select(canonical_smiles,Entry) %>% unique()
# Read training data (it has pbs and dmso)
long_drug_target <-  read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv')
colnames(long_drug_target)[1] <- 'canonical_smiles'
long_drug_target <- long_drug_target %>% gather('Entry','value',-canonical_smiles)
long_drug_target <- long_drug_target %>% filter(value==1) %>% select(-value) %>% unique()
long_drug_target <- rbind(long_drug_target,data.frame(canonical_smiles='Ctrl_PBS',Entry='None'))
cmap_drugs_unnested <- rbind(cmap_drugs_unnested,long_drug_target)
cmap_drugs_unnested <- cmap_drugs_unnested %>% unique()

## Manually add dmso from drugbank and also add more targets for other drugs too

### SOS : The following commented code was used to parse the xml file of DrugBank.
### Unfortunately this read_drugbank_xml_db function is not supported anymore in the dbparser library
### Instead we provide the extracted file in .rds format in the data folder. 
### We provided the code for transparency and completeness ,
### But please skip the commented lines and just read the .rds file.
df_drugbank <- readRDS('../data/df_drugbank.rds')
#read_drugbank_xml_db("../../full database.xml")
#drugs <- drugs()
#snp_effects <- drugs$snp_effects
#drug_interactions <- drugs$interactions
#snp_effects <- snp_effects %>% select(`protein-name`,`gene-symbol`,`uniprot-id`,`pubmed-id`,parent_key)
#drug_interactions <- drug_interactions %>% select(-description)
#drug_interactions <- left_join(drug_interactions,snp_effects)
#df_drugbank <- drugs$general_information
#df_drugbank <- left_join(df_drugbank,drug_interactions,by=c('primary_key'='drugbank-id','name'='name'))
#df_drugbank$name <- tolower(df_drugbank$name)
#df_drugbank <- df_drugbank %>%
#  select(c('drugbank_id'='primary_key'),'name','cas_number',
#         c('target_symbol'='gene-symbol'),c('target_uniprot'='uniprot-id')) %>% 
#  filter(!is.na(target_symbol) | !is.na(target_uniprot)) %>% unique()

df_drugbank_cmap <- left_join(cmap_drugs %>% select(-Target) %>% mutate(cmap_name=tolower(pert_iname)) %>% unique(),
                              df_drugbank,
                              by=c('pert_iname'='name')) %>%
  select(-target_symbol,-cas_number) %>% filter(!is.na(target_uniprot)) %>% 
  filter(target_uniprot!='') %>% filter(target_uniprot!=' ')%>%
  unique()
df_drugbank_cmap <-  df_drugbank_cmap %>% select("pert_iname","canonical_smiles",
                                                 "Entry"='target_uniprot') %>% unique()
df_drugbank_cmap <-  df_drugbank_cmap %>% filter(!is.na(Entry)) %>% filter(Entry!='') %>% filter(Entry!=' ')
df_drugbank_cmap$Entry <- toupper(df_drugbank_cmap$Entry) 
df_drugbank_cmap <- df_drugbank_cmap %>% select(-pert_iname) %>% unique()
df_drugbank_not_in_cmap <- df_drugbank %>% filter(!(name %in% tolower(cmap_drugs$pert_iname))) %>%
  select(c('canonical_smiles'='name'),c('Entry'='target_uniprot')) %>% unique()

drug_targets_space <- rbind(cmap_drugs_unnested,df_drugbank_cmap,df_drugbank_not_in_cmap) %>% unique()


uniq_targets <- unique(drug_targets_space$Entry)
binary_classes <- function(sig,df,targets){

  df <- df %>% select(canonical_smiles,Entry) %>% unique()
  df <- df %>% filter(canonical_smiles==sig)
  m <- as.matrix(df[,2])
  ind <- which(targets %in% m)
  n <- targets
  targets <- rep(0,length(targets))
  targets[ind] <- 1
  names(targets) <- n
  return(targets)
}
drug_targets <- t(sapply(unique(drug_targets_space$canonical_smiles),binary_classes,drug_targets_space,uniq_targets))
ind <- which(colnames(drug_targets)=='None')
drug_targets <- drug_targets[,-ind]
### Save the file and for spead load it for the rest of the analysis
# saveRDS(drug_targets,'../preprocessing/preprocessed_data/drug_targets_space.rds')

#### Compare inferred with drugbank------------------------------
# ensemble_table <- readRDS('../results/A375_ensembles/ensemble_table.rds')
# new_predictions <- readRDS(../results/A375_ensembles/new_predictions.rds')
# pvals <- readRDS('../results/A375_ensembles/pvals.rds')
# drug_bank_tp <- readRDS('../results/A375_ensembles/drug_bank_tp.rds')
# ensemble_table_null <- readRDS(../results/A375_ensembles/ensemble_table_null.rds')
no_models <- 50
drug_targets <- readRDS('../preprocessing/preprocessed_data/drug_targets_space.rds')
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]

### See if I consider as ground-truth cmap plus drugbank
### how precision and recall change with respect to including predictions
### based on how many model they appear
no_models <- 50
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
cmap_prior <- merged_interactions_all %>% select(drug,variable,`Prior knowledge`) %>% unique()
cmap_prior <- cmap_prior %>% filter(`Prior knowledge`=='No interaction') %>% select(-`Prior knowledge`) %>% unique()
drug_targets_prior <- as.data.frame(drug_targets) %>% rownames_to_column('drug') %>% 
  gather('variable','Prior knowledge',-drug)
merged_interactions_all <- merged_interactions_all %>% select(-`Prior knowledge`)
merged_interactions_all <- left_join(merged_interactions_all,drug_targets_prior)
merged_interactions_all <- merged_interactions_all %>% filter(!is.na(`Prior knowledge`))
merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = 1*(Inferred=='Interaction'))
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
  mutate(model_counts = sum(Inferred==1)) %>% ungroup()
merged_interactions_all <- merged_interactions_all %>% mutate(model_counts=model_counts/no_models)
print(all(merged_interactions_all %>% filter(model_counts==0) %>% mutate(logic = (Inferred==`Prior knowledge`)) %>% select(logic)))

thresholds <- seq(1,no_models)/no_models
prec <- NULL
rec <- NULL
acc <- NULL
F1 <- NULL
pvals <- NULL
for (i in 1:length(thresholds)){
  # For those inferred keep everything above the threshold
  # Those never inferred in the model are indeed never in prior knowledge
  th <- thresholds[i]
  #thresholded_predictions <- distinct(merged_interactions_all %>% filter(model_counts>=th | model_counts==0) %>% 
  #  select("drug","variable","Prior knowledge","Inferred"))
  thresholded_predictions <- distinct(merged_interactions_all %>% mutate(Inferred = ifelse(model_counts>=th,1,0)) %>% 
                                        select("drug","variable","Prior knowledge","Inferred"))
  confusion <- confusionMatrix(data=factor(thresholded_predictions$Inferred,levels = c(1,0)), 
                               reference = factor(thresholded_predictions$`Prior knowledge`,levels = c(1,0)))
  #print(confusion$table)
  #print(confusion$byClass)
  acc[i] <- confusion$overall['Accuracy']
  prec[i] <- confusion$byClass['Precision']
  rec[i] <- confusion$byClass['Recall']
  F1[i] <- confusion$byClass['F1']
  pvals[i] <- confusion$overall['AccuracyPValue']
}
df_res <- data.frame(thresholds=thresholds,accuracy = acc,precision=prec,recall=rec,F1=F1,pvalue = pvals)
df_res <- df_res %>% gather('metric','value',-thresholds)

## Calculate p.values and F1 score of individual models
p_vals_one <- NULL
f1_one_model <- NULL
acc_one_model <- NULL
for (i in 1:no_models){
  model_predictions <- distinct(merged_interactions_all %>% filter(model_no==i-1) %>% 
                                  select("drug","variable","Prior knowledge","Inferred"))
  confusion <- confusionMatrix(data=factor(model_predictions$Inferred,levels = c(1,0)), 
                               reference = factor(model_predictions$`Prior knowledge`,levels = c(1,0)))
  p_vals_one[i] <- confusion$overall['AccuracyPValue']
  f1_one_model[i] <- confusion$byClass['F1']
  acc_one_model[i] <- confusion$overall['Accuracy']
}
hist(p_vals_one)
# hist(f1_one_model)
hist(acc_one_model)

# visualize
p1 <- ggplot(df_res %>% filter(metric!='pvalue'),
             aes(x=thresholds,y=value,color=metric)) + 
  geom_point() +geom_line(linewidth=1) + 
  xlab('minimum frequency to consider an interaction') + 
  ggtitle('Discovery of new interactions as the frequency of appearing in multiple models increases')+
  theme_minimal() +
  geom_ribbon(inherit.aes = FALSE,
              xmin=0,xmax=1,
              aes(x=seq(0,1,length.out=length(thresholds)),
                  ymin = min(acc_one_model), 
                  ymax = max(acc_one_model)), 
              fill = "red", alpha = 0.1) +  # Shaded area
  geom_hline(yintercept = min(acc_one_model), 
             color = "red", linwidth = 1) +
  geom_hline(yintercept = max(acc_one_model), 
             color = "red", linwidth = 1) + 
  annotate('text',x=0.01,y=mean(c(min(acc_one_model),max(acc_one_model))),
           label=paste0("Range:","\n","accuracy of","\n","individual models"), 
           hjust = 0 , size=7)+
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24),
        legend.text = element_text(family = 'Arial',size = 24),
        plot.title = element_text(family = 'Arial',size = 20, hjust = 0.5))
print(p1)

ggsave('../figures/figure2E.eps',
       plot = p1,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 8,
       units = "in",
       dpi = 600)

p2 <- ggplot(df_res %>% filter(metric=='pvalue'),
             aes(x=thresholds,y=value)) + 
  #ylim(c(0,1))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #ylim(c(0,1))+
  geom_point() +geom_line(linewidth=1,linetype='dashed') + 
  xlab('minimum frequency to consider prediction') + 
  ylab('p-value') +
  ggtitle('P-value for accuracy to be greater than no-information rate')+
  theme_minimal() +
  geom_hline(yintercept = 0.01,linetype='dashed',linewidth=1.5,color='red')+
  annotate('text',x=0.25,y=1e-07,
           label='p-value = 0.01',size=10)+
  #geom_hline(yintercept = mean(p_vals_one),linetype='dashed',linewidth=1.5,color='orange')+
  # annotate('text',x=0.25,y=10000,
  #          label='p-value of individual models',size=10)+
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24),
        legend.text = element_text(family = 'Arial',size = 24),
        plot.title = element_text(family = 'Arial',size = 22, hjust = 0.5))
print(p2)

ggsave('../figures/figure2D.eps',
       plot = p2,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 8,
       units = "in",
       dpi = 600)

# plot confusion matrix for all ensembles A375-----------------
no_models <- 50
#th <- 40/no_models
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
ensemble_table <- NULL
for (i in 1:no_models){
  merged_interactions <- merged_interactions_all %>% filter(model_no == i-1) %>% select(-model_no)
  confusion <- confusionMatrix(data=as.factor(merged_interactions$Inferred), 
                               reference = as.factor(merged_interactions$`Prior knowledge`))
  ensemble_table[[i]] <-  confusion$table
  message(paste0('Finished model ',i))
}
ensemble_mat <- do.call(cbind,ensemble_table)
ensemble_mat <- array(ensemble_mat,c(dim=dim(ensemble_table[[1]]),length(ensemble_table)))
ensemble_mat_mean <- apply(ensemble_mat, c(1,2), mean, na.rm = TRUE)
colnames(ensemble_mat_mean) <- colnames(ensemble_table[[1]])
rownames(ensemble_mat_mean) <- rownames(ensemble_table[[1]])
ensemble_mat_sd <- apply(ensemble_mat, c(1,2), sd, na.rm = TRUE)
colnames(ensemble_mat_sd) <- colnames(ensemble_table[[1]])
rownames(ensemble_mat_sd) <- rownames(ensemble_table[[1]])
ensemble_mat_mean <- melt(ensemble_mat_mean, variable.name = c("Reference", "Prediction"), value.name = "Count")
ensemble_mat_sd <- melt(ensemble_mat_sd, variable.name = c("Reference", "Prediction"), value.name = "Count_sd")
ensemble_final <- left_join(ensemble_mat_mean,ensemble_mat_sd)
plot_confusion <- ggplot(ensemble_final %>% mutate(Count = ifelse(Count>=1000,round(Count,-3),ifelse(Count>=10,round(Count,-1),round(Count,0)))) %>%
                           mutate(Count_se = Count_sd/sqrt(50)) %>%
                           mutate(Count_se = ifelse(Count_se>=1000,round(Count_se,-3),ifelse(Count_se>=10,round(Count_se,-1),round(Count_se,0)))) , 
                         aes(x = Var1, y = Var2, fill = Count)) +
  geom_tile() +
  scale_fill_viridis(trans='log10', limits = c(1,round(max(ensemble_final$Count),-3)),breaks=c(1,10,100,1000,10000)) +
  labs(fill='Count')+
  geom_text(aes(label = paste0('~',paste(Count, Count_se, sep = "\u00B1"))), size = 14) +
  labs(title = "Confusion matrix of inferred and prior knowledge drug-target interactions",
       x = "Inferred",
       y = "Prior knowledge") +
  theme_minimal() +
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24, face = "bold"),
        legend.text = element_text(family = 'Arial',size = 20),
        plot.title = element_text(family = 'Arial',size = 24, hjust = 0.5))
print(plot_confusion)

ggsave('../figures/figure2B.eps',
       plot = plot_confusion,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)
## Plot drug-target matrix filled by frequency
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
  mutate(model_counts = sum(Inferred=='Interaction')) %>% ungroup()

ann_text <- ensemble_final %>% filter(Var1=='Interaction') %>% select(c('Prior knowledge'='Var2'),Count,Count_sd)
ann_text <- ann_text %>% mutate(Count_se = Count_sd/sqrt(no_models)) %>% select(`Prior knowledge`,Count,Count_se)
total_prior <- distinct(merged_interactions_all %>% select(drug,variable,`Prior knowledge`))
total_prior_negatives <-  sum(total_prior$`Prior knowledge`!='Interaction')
total_prior <- sum(total_prior$`Prior knowledge`=='Interaction')
ann_text <- ann_text %>% mutate(total_prior = ifelse(`Prior knowledge`=='Interaction',total_prior,total_prior_negatives))
ann_text <- ann_text %>% mutate(Count=Count/total_prior) %>% mutate(Count_se=Count_se/total_prior)
ann_text <- ann_text %>% mutate(lab=paste0('~',paste(round(100*Count,2),round(100*Count_se,2), sep = "\u00B1"),'%'))
ann_text <- ann_text %>% mutate(`Prior knowledge`=paste0('Prior : ',`Prior knowledge`))
ann_text <- ann_text %>% select(`Prior knowledge`,lab)
ann_text <- ann_text %>% mutate(lab = ifelse(`Prior knowledge`=='Prior : Interaction',
                                             paste0('retrieved',"\n",lab,"\n",' of prior knowledge'),
                                             paste0(lab,"\n",' of prior negatives')))

p <- ggplot(distinct(merged_interactions_all %>% filter(Inferred=='Interaction') %>%
                       select(drug,variable,`Prior knowledge`,model_counts) %>% mutate(`Prior knowledge`=paste0('Prior : ',`Prior knowledge`))),
            aes(x=model_counts)) +
  geom_histogram(aes(y = ..count..),fill='#0F52BA') +
  scale_x_continuous(limits = c(0,NA))+
  facet_wrap(~`Prior knowledge`,scales = "free") + xlab('number of models')+
  ggtitle('Inferred interactions by the model') +
  #scale_y_log10()+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(panel.grid.major.y = element_line(linetype = 'dashed',colour = 'lightgrey'),
        panel.grid.minor.y = element_line(linetype = 'dashed',colour = 'lightgrey'),
        panel.grid.major.x = element_line(linetype = 'dashed',colour = 'lightgrey'),
        panel.grid.minor.x = element_line(linetype = 'dashed',colour = 'lightgrey'),
        plot.title = element_text(size=20,hjust = 0.5))
p <- p + geom_text(data = ann_text,mapping = aes(x = c(20,30), y = c(280,5200), label = lab),
                   size=9
)
print(p)

ggsave('../figures/figure2C.eps',
       plot = p,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)
