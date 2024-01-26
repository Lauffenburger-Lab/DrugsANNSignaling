library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(viridis)
library(patchwork)
library(RColorBrewer)

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
cmap_drugs <- cmap_drugs %>% mutate(Disease.Area = str_split(Disease.Area,pattern = ", ")) %>% unnest(Disease.Area) %>% unique()

top_moas <- cmap_drugs %>% group_by(MOA) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup() 
top_moas <- top_moas %>% select(MOA,counts) %>% unique()
top_moas <- top_moas %>% arrange(-counts)
all_moas <- top_moas %>% select(MOA,counts) %>% unique()
top_moas <- top_moas[1:10,]

top_diseases <- cmap_drugs %>% group_by(Disease.Area) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_diseases <- cmap_drugs %>% group_by(Disease.Area) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_diseases <- top_diseases %>% select(Disease.Area,counts) %>% unique()
top_diseases <- top_diseases %>% filter(Disease.Area!='')
top_diseases <- top_diseases %>% arrange(-counts)
all_diseases <- top_diseases %>% select(Disease.Area,counts) %>% unique()
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
tmp <- data_all_plot_moas %>% group_by(MOA) %>% mutate(mean_r = mean(r,na.rm = T)) %>% ungroup()
tmp <- tmp %>% select(MOA,mean_r) %>% unique()
data_all_plot_moas$MOA <- factor(data_all_plot_moas$MOA,levels = tmp$MOA[order(tmp$mean_r)])
moa_colors <- get_palette(palette = "default", length(unique(data_all_plot_moas$MOA)))
# top_moas2 <- data_all_plot_moas %>% group_by(MOA) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup() 
# top_moas2 <- top_moas2 %>% select(MOA,counts) %>% unique()
# top_moas2 <- top_moas2 %>% arrange(-counts)
# top_moas2 <- top_moas2[1:10,]
p_moa <- ggboxplot(data_all_plot_moas %>% select(MOA,TF,r) %>% unique(),
          x='MOA',y='r',color='MOA',add='jitter')+
  scale_color_manual(values = moa_colors) +
  ylab('pearson`s r') + xlab('mechanism of action') +
  theme(text=element_text(family = 'Arial',size=24),
        axis.text.y = element_text(family = 'Arial',size=14),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = 'none')+
  stat_compare_means(method = 'kruskal.test',label.y = -0.4,label.x = 48.5,size=5)+
  coord_flip()
moa_colors_zoomed <- moa_colors[which(levels(data_all_plot_moas$MOA) %in% top_moas$MOA)]
p_moa_zoomed <- ggboxplot(data_all_plot_moas %>% select(MOA,TF,r) %>% filter(MOA %in% top_moas$MOA) %>% unique(),
                           x='MOA',y='r',color='MOA',add='jitter')+
  scale_color_manual(values = moa_colors_zoomed) +
  ylab('pearson`s r')+ xlab('mechanism of action') +
  theme(text=element_text(family = 'Arial',size=24),
        axis.text.y = element_text(family = 'Arial',size=24),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = 'none')+
  stat_compare_means(method = 'kruskal.test',label.y = -0.2,label.x = 7.3,size=8)+
  coord_flip()
print(p_moa)
ggsave('../article_supplementary_info/performance_per_moa.eps',
       device = cairo_ps,
       plot = p_moa,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)
print(p_moa_zoomed)
ggsave('../article_supplementary_info/performance_per_moa_zoomed.eps',
       device = cairo_ps,
       plot = p_moa_zoomed,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)
## Check per disease area
data_all_plot_diseases <- data_all %>% group_by(Disease.Area,TF) %>% mutate(r = cor(prediction,activity)) %>% ungroup()
data_all_plot_diseases <- data_all_plot_diseases %>% filter(!is.na(r))
tmp <- data_all_plot_diseases %>% group_by(Disease.Area) %>% mutate(mean_r = mean(r,na.rm = T)) %>% ungroup()
tmp <- tmp %>% select(Disease.Area,mean_r) %>% unique()
data_all_plot_diseases$Disease.Area <- factor(data_all_plot_diseases$Disease.Area,levels = tmp$Disease.Area[order(tmp$mean_r)])
p_disease <- ggboxplot(data_all_plot_diseases %>% filter(Disease.Area %in% top_diseases$Disease.Area) %>% select(Disease.Area,TF,r) %>% unique(),
          x='Disease.Area',y='r',color='Disease.Area',add='jitter')+
  ylab('pearson`s r') + xlab('disease area') +
  theme(text=element_text(family = 'Arial',size=24),
        axis.text.y = element_text(family = 'Arial',size=24),
        legend.position = 'none')+
  stat_compare_means(method = 'kruskal.test',label.y = -0.1,label.x = 10,size=8)+
  coord_flip()
print(p_disease)
ggsave('../article_supplementary_info/performance_per_disease.eps',
       device = cairo_ps,
       plot = p_disease,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)
### Check off-target distribution across MoA and Disease.Area-------------------------------------
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
cmap_drugs <- cmap_drugs %>% mutate(Disease.Area = str_split(Disease.Area,pattern = ", ")) %>% unnest(Disease.Area) %>% unique()

top_moas <- cmap_drugs %>% group_by(MOA) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup() 
top_moas <- top_moas %>% select(MOA,counts) %>% unique()
top_moas <- top_moas %>% arrange(-counts)
all_moas <- top_moas %>% select(MOA,counts) %>% unique()
top_moas <- top_moas[1:10,]

top_diseases <- cmap_drugs %>% group_by(Disease.Area) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_diseases <- cmap_drugs %>% group_by(Disease.Area) %>% mutate(counts = n_distinct(canonical_smiles)) %>% ungroup()
top_diseases <- top_diseases %>% select(Disease.Area,counts) %>% unique()
top_diseases <- top_diseases %>% filter(Disease.Area!='')
top_diseases <- top_diseases %>% arrange(-counts)
all_diseases <- top_diseases %>% select(Disease.Area,counts) %>% unique()
top_diseases <- top_diseases[1:10,]
performance <- data.table::fread('../results/A375_ensembles/meanCorrPerTFEnsembleVal_lamda6.csv',
                                 header=T)
performance <- performance %>% dplyr::select(-model) %>% unique()
performance <- performance %>% group_by(TF) %>% mutate(mean_r=mean(r)) %>%
  ungroup() %>% dplyr::select(-cell,-r) %>% unique()
base_cell_performance <- data.table::fread('../results/A375_ensembles/A375TrainEnsemblePerformance.csv')
colnames(base_cell_performance)[1] <- 'TF'
Delta <- data.table::fread('../results/A375_ensembles/DeltaTF1.csv',header=T)
# Load TF activities
TFoutput <- read.delim('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv') %>% filter(X %in% Delta$V1) %>%
  column_to_rownames('X') %>% rownames_to_column('sample')
gc()
TFoutput <- TFoutput %>% gather('TF','activity',-sample) 
Delta <- Delta  %>% column_to_rownames('V1') %>%  rownames_to_column('sample') %>% gather('TF','delta',-sample)
# merge everything
df <- left_join(TFoutput,Delta,by=c('sample','TF'))
df <- left_join(df,performance)
df <- left_join(df,base_cell_performance)
df <- df %>% mutate(score=0.5*(mean_r+r))
# Load conditions to get rid of DMSO
conditions <- data.table::fread('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv',sep = "\t") %>% column_to_rownames('V1')
conditions <- conditions %>% rownames_to_column('sample') %>% gather('drug','value',-sample) %>% filter(value>0) %>%
  select(-value) %>% unique()
conditions <- conditions %>% filter(sample %in% df$sample) %>% filter(drug!='CS(C)=O')
annotation <- read.delim('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv') %>% dplyr::select(c('TF'='code'),name)
annotation <- annotation %>% filter(TF %in% df$TF)
df <- left_join(df,annotation)
df <- left_join(df,conditions)
df <- df %>% filter(!is.na(drug))

### Merge the moa data with effect data
effect_data <- left_join(cmap_drugs,
                         df %>% filter(!is.na(drug))  %>% filter(mean_r>0.4 & score>=0.5) %>% 
                           select(sample,c('canonical_smiles'='drug'),TF,name,delta,activity) %>% unique())
effect_data <- distinct(effect_data %>% select(sample,canonical_smiles,TF,name,MOA,Disease.Area,Phase,delta,activity))

effect_data_moa_plot <- effect_data
tmp <- effect_data_moa_plot %>% group_by(MOA) %>% mutate(mean_effect = mean(abs(delta),na.rm = T)) %>% ungroup()
tmp <- tmp %>% select(MOA,mean_effect) %>% unique()
effect_data_moa_plot$MOA <- factor(effect_data_moa_plot$MOA,levels = tmp$MOA[order(tmp$mean_effect)])
moa_colors <- get_palette(palette = "default", length(unique(effect_data_moa_plot$MOA)))
p_effect_moa <- ggboxplot(effect_data_moa_plot %>% select(MOA,TF,delta) %>% 
                            mutate(delta=abs(delta)) %>% unique(),
                   x='MOA',y='delta',color='MOA',add='jitter')+
  scale_color_manual(values = moa_colors) +
  ylab(expression(abs(Delta*"TF"))) + xlab('mechanism of action') +
  theme(text=element_text(family = 'Arial',size=24),
        axis.text.y = element_blank(),
        # axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = 'none')+
  stat_compare_means(method = 'kruskal.test',label.y = 0.25,label.x = 5,size=5)+
  # stat_compare_means(ref.group = 'lysophospholipid receptor antagonist',label.y = 0.3,size=6,label = 'p.signif')+
  coord_flip()
moa_colors_zoomed <- moa_colors[which(levels(effect_data_moa_plot$MOA) %in% top_moas$MOA)]
p_effect_moa_zoomed <- ggboxplot(effect_data_moa_plot %>% select(MOA,TF,delta) %>% 
                                           mutate(delta=abs(delta)) %>% filter(MOA %in% top_moas$MOA) %>% unique(),
                           x='MOA',y='delta',color='MOA',add='jitter')+
  scale_color_manual(values = moa_colors_zoomed) +
  ylab(expression(abs(Delta*"TF"))) + xlab('mechanism of action') +
  theme(text=element_text(family = 'Arial',size=24),
        # axis.text.y = element_text(family = 'Arial',size=24),
        # axis.title.y = element_blank(),
        # axis.text.x = element_blank(),
        legend.position = 'none')+
  stat_compare_means(method = 'kruskal.test',label.y = 0.27,label.x = 2,size=8)+
  coord_flip()
print(p_effect_moa)
ggsave('../article_supplementary_info/effect_per_moa.eps',
       device = cairo_ps,
       plot = p_effect_moa,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)

print(p_effect_moa_zoomed)
ggsave('../article_supplementary_info/effect_per_moa_zoomed.eps',
       device = cairo_ps,
       plot = p_effect_moa_zoomed,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)

## Check per disease area
effect_data_disease_plot <- effect_data
effect_data_disease_plot <- effect_data_disease_plot %>% filter(!is.na(Disease.Area))
effect_data_disease_plot <- effect_data_disease_plot %>% filter(Disease.Area!='')
tmp <- effect_data_disease_plot %>% group_by(Disease.Area) %>% mutate(mean_effect = mean(abs(delta),na.rm = T)) %>% ungroup()
tmp <- tmp %>% select(Disease.Area,mean_effect) %>% unique()
effect_data_disease_plot$Disease.Area <- factor(effect_data_disease_plot$Disease.Area,levels = tmp$Disease.Area[order(tmp$mean_effect)])
p_disease <- ggboxplot(effect_data_disease_plot %>% select(Disease.Area,TF,delta) %>% 
                         mutate(delta=abs(delta)) %>% unique(),
                       x='Disease.Area',y='delta',color='Disease.Area',add='jitter')+
  ylab(expression(abs(Delta*"TF"))) + xlab('disease area') +
  theme(text=element_text(family = 'Arial',size=24),
        axis.text.y = element_text(family = 'Arial',size=24),
        legend.position = 'none')+
  # stat_compare_means(ref.group = 'urology',label.y = 0.3,size=6,label = 'p.signif')+
  stat_compare_means(method = 'kruskal.test',label.y = 0.25,label.x = 2,size=8)+
  coord_flip()
print(p_disease)
ggsave('../article_supplementary_info/effect_per_disease.eps',
       device = cairo_ps,
       plot = p_disease,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)
### Test if there are large off-target effects in disease areas
one_area_drugs <- effect_data_disease_plot %>% group_by(canonical_smiles) %>% mutate(disease_counts = n_distinct(Disease.Area)) %>% ungroup() %>%
  filter(disease_counts==1) %>% select(-disease_counts) %>% unique()
interest <- unique(one_area_drugs$Disease.Area) #c('genetics','hematologic malignancy','oncology')
pvalue <- NULL
i <- 1
for (area in interest){
  tmp <- one_area_drugs %>% filter(Disease.Area==area) %>% unique()
  test <- wilcox.test(abs(tmp$delta))
  pvalue[i] <- test$p.value
  i <- i+1
}
p.adj <- p.adjust(pvalue,'bonferroni')
names(p.adj) <- interest
print(p.adj)
