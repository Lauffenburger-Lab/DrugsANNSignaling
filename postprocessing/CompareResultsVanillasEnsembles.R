library(tidyverse)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(RColorBrewer)
library(patchwork)
library(ggforce)
library(ggsignif)
library(ggstatsplot)
library(factoextra)
library(gridExtra)
# library(cowplot)

### Load data per TF analysis-----------
files <- list.files('../results/FinalEnsemble/test/')
files <- files[grep('.csv',files)]

file_vanilla <- list.files('../results/vanilla_ensembles/',recursive = T)
file_vanilla <- file_vanilla[grep('.csv',file_vanilla)]
all_files <- list(files,file_vanilla)

# Load training performance
df_corr_train_ensemble <- data.frame()
for (i in 1:length(all_files)){
  if (i==1){
    files_corr <- all_files[[i]][grep('trainEnsemblePerformance',all_files[[i]])]
  }else{
    files_corr <- all_files[[i]][grep('train_Performance',all_files[[i]])]
  }
  for (file in files_corr){
    if (i==1){
      cell <- str_split_fixed(file,'_',2)[1,2]
      cell <- str_split_fixed(cell,'[.]',2)[1,1]
      file <- paste0('../results/FinalEnsemble/test/',file)
      tmp <- data.table::fread(file,header = T)
      colnames(tmp)[1] <- 'TF'
      colnames(tmp)[2] <- 'r'
      tmp <- tmp %>% mutate(model='DT-LEMBAS')
      tmp <- tmp %>% mutate(cell=cell)
      names <- colnames(tmp)
    }else{
      model <- str_split_fixed(file,'/',4)[1,1]
      file <- paste0('../results/vanilla_ensembles/',file)
      tmp <- data.table::fread(file,) %>% mutate(model=model)
      tmp <- tmp %>% mutate(cell=V1) %>% select(-V1)
      tmp <- tmp %>% gather('TF','r',-model,-cell)
      tmp <- tmp %>% select(all_of(names))
    }
    df_corr_train_ensemble <- rbind(df_corr_train_ensemble,tmp)
  }
}
df_corr_train <-  df_corr_train_ensemble %>% select(model,cell,TF,r)

# Load validation correlation data
df_corr_val <- data.frame()
for (i in 1:length(all_files)){
  if (i==1){
    files_corr <- all_files[[i]][grep('valEnsemblePerformance',all_files[[i]])]
  }else{
    files_corr <- all_files[[i]][grep('val_Performance',all_files[[i]])]
  }
  for (file in files_corr){
    if (i==1){
      cell <- str_split_fixed(file,'_',2)[1,2]
      cell <- str_split_fixed(cell,'[.]',2)[1,1]
      file <- paste0('../results/FinalEnsemble/test/',file)
      tmp <- data.table::fread(file,header = T)
      colnames(tmp)[1] <- 'TF'
      colnames(tmp)[2] <- 'r'
      tmp <- tmp %>% mutate(model='DT-LEMBAS')
      tmp <- tmp %>% mutate(cell=cell)
      names <- colnames(tmp)
    }else{
      model <- str_split_fixed(file,'/',4)[1,1]
      file <- paste0('../results/vanilla_ensembles/',file)
      tmp <- data.table::fread(file,) %>% mutate(model=model)
      tmp <- tmp %>% mutate(cell=V1) %>% select(-V1)
      tmp <- tmp %>% gather('TF','r',-model,-cell)
      tmp <- tmp %>% select(all_of(names))
    }
    df_corr_val <- rbind(df_corr_val,tmp)
  }
}
df_pvalue_val <- data.frame()
for (i in 1:length(all_files)){
  if (i==1){
    files_pvalue <- all_files[[i]][grep('valEnsemblePvalues',all_files[[i]])]
    for (file in files_pvalue){
      cell <- str_split_fixed(file,'_',2)[1,2]
      cell <- str_split_fixed(cell,'[.]',2)[1,1]
      file <- paste0('../results/FinalEnsemble/test/',file)
      tmp <- data.table::fread(file,header = T)
      colnames(tmp)[1] <- 'TF'
      colnames(tmp)[2] <- 'pvalue'
      tmp <- tmp %>% mutate(p.adj = p.adjust(pvalue,'BH'))
      tmp <- tmp %>% mutate(cell=cell)
      df_pvalue_val <- rbind(df_pvalue_val,tmp)
    }
  }
}

# Add random shuffle
random_files <- list.files('../results/RandomSuffleY/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesY <- random_files[grep('Performance',random_files)]
random_filesY <- random_filesY[which(!grepl('class',random_filesY))]
random_filesY <- random_filesY[which(!grepl('quantile',random_filesY))]
random_filesY <- random_filesY[which(!grepl('R2',random_filesY))]
random_files <- list.files('../results/RandomSuffleX/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesX <- random_files[grep('Performance',random_files)]
random_filesX <- random_filesX[which(!grepl('class',random_filesX))]
random_filesX <- random_filesX[which(!grepl('quantile',random_filesX))]
random_filesX <- random_filesX[which(!grepl('R2',random_filesX))]
all_random_files <- list(random_filesX,random_filesY)

df_corr_val_random <- data.frame()
for (i in 1:length(all_random_files)){
  files_corr <- all_random_files[[i]][grep('valPerformance',all_random_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==1) {
      file <- paste0('../results/RandomSuffleX/test/',file)
    }else{
      file <- paste0('../results/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'model'
    # tmp <- tmp %>% select(- lamda) %>% unique()
    # if (i==1){
    #   tmp <- tmp %>% mutate(model = 'shuffle X')
    # }else{
    #   tmp <- tmp %>% mutate(model = 'shuffle Y')
    # }
    if (i==1){
      tmp <- tmp %>% mutate(model = gsub("Y", "X", model))
    }
    tmp <- tmp %>% mutate(cell=cell)
    tmp <- tmp %>% gather('TF','r',-model,-cell)
    tmp <- tmp %>% select(all_of(names))
    df_corr_val_random <- rbind(df_corr_val_random,tmp)
  }
}

# Add random shuffle TRAINING
df_corr_train_random <- data.frame()
for (i in 1:length(all_random_files)){
  files_corr <- all_random_files[[i]][grep('trainPerformance',all_random_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==1) {
      file <- paste0('../results/RandomSuffleX/test/',file)
    }else{
      file <- paste0('../results/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'model'
    # tmp <- tmp %>% select(- lamda) %>% unique()
    # if (i==1){
    #   tmp <- tmp %>% mutate(model = 'shuffle X')
    # }else{
    #   tmp <- tmp %>% mutate(model = 'shuffle Y')
    # }
    if (i==1){
      tmp <- tmp %>% mutate(model = gsub("Y", "X", model))
    }
    tmp <- tmp %>% mutate(cell=cell)
    tmp <- tmp %>% gather('TF','r',-model,-cell)
    tmp <- tmp %>% select(all_of(names))
    df_corr_train_random <- rbind(df_corr_train_random,tmp)
  }
}

# Visualize
df_corr_val_random <- df_corr_val_random %>% group_by(model,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val <- df_corr_val %>% group_by(model,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val_random <- df_corr_val_random %>% group_by(model,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val <- rbind(df_corr_val,df_corr_val_random)
# df_corr_val$model <- factor(df_corr_val$model,
#                             levels = c('DT-LEMBAS',
#                                        'ANN','GCNN',
#                                        'SVM','KNN',
#                                        'shuffle X','shuffle Y'))
p3_2_1 <- ggboxplot(df_corr_val %>% mutate(model=ifelse(grepl('shuffleY',model),'shuffle Y',ifelse(grepl('shuffleX',model),'shuffle X',model))) %>%
                      mutate(model = factor(model,levels= c('DT-LEMBAS','ANN','GCNN','SVM','KNN','shuffle X','shuffle Y'))), 
                    x = "model", y = "r",outlier.shape = NA,add='jitter')+
  #ggtitle('Ensembles` performance comparison with standard machine learning approaches') +
  #geom_jitter(aes(color=cell),fill='black',alpha=0.2)+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5),legend.position = '') +
  xlab('model') + ylab('per TFs pearson`s r')+ 
  scale_y_continuous(breaks = c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,seq(0.1,1,0.1)),limits = c(-0.635,1)) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_1 <- p3_2_1 + stat_compare_means(comparisons = list(c('DT-LEMBAS','ANN'),
                                                         c('DT-LEMBAS','GCNN'),
                                                         c('DT-LEMBAS','SVM'),
                                                         c('DT-LEMBAS','KNN'),
                                                         c('DT-LEMBAS','shuffle X'),
                                                         c('DT-LEMBAS','shuffle Y')),
                                      label = 'p.format',
                                      label.y = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95),
                                      method = 'wilcox.test',
                                      tip.length=0.05)
print(p3_2_1)
png('../MIT/LauffenburgerLab/SBHD 2023/model_vanillas_ensembles_comparison_validation_perTF_1.png',units = 'in',width = 12,height = 12,res=600)
print(p3_2_1)
dev.off()
p3_2_2 <- ggboxplot(df_corr_val %>% mutate(model=ifelse(grepl('shuffleY',model),'shuffle Y',ifelse(grepl('shuffleX',model),'shuffle X',model))) %>%
                      group_by(model,cell) %>%mutate(mean_pear = ifelse(grepl('shuffle',model),mean(mean_pear),mean_pear)) %>% ungroup() %>%
                      select(model,mean_pear,cell) %>% unique()%>% 
                      mutate(model = factor(model,
                                            levels= c('DT-LEMBAS','ANN','GCNN','SVM','KNN','shuffle X','shuffle Y'))), 
                    x = "model", y = "mean_pear",
                    add='jitter',outlier.shape = NA)+
  ggtitle('') + #geom_jitter(aes(color=cell))+
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5),legend.position = 'bottom') +
  xlab('model') + ylab('average per TFs pearson`s r') + scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_2 <- p3_2_2 + stat_compare_means(comparisons = list(c('DT-LEMBAS','ANN'),
                                                         c('DT-LEMBAS','GCNN'),
                                                         c('DT-LEMBAS','SVM'),
                                                         c('DT-LEMBAS','KNN'),
                                                         c('DT-LEMBAS','shuffle X'),
                                                         c('DT-LEMBAS','shuffle Y')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05)
print(p3_2_1/p3_2_2)
png('model_vanillas_ensembles_comparison_validation_perTF.png',units = 'in',width = 12,height = 12,res=600)
print(p3_2_1/p3_2_2)
dev.off()

lembas <- distinct(df_corr_val %>% filter(model=='DT-LEMBAS') %>% select(-model))
p3_2_3 <- ggboxplot(lembas, x = "cell", y = "r",
                    add='jitter',outlier.shape = NA)+
  ggtitle('DT-LEMBAS model per cell line') + #geom_jitter(aes(color=cell))+
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5),legend.position = 'bottom') +
  xlab('cell line') + ylab('per TFs pearson`s r') + scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
print(p3_2_1/p3_2_3)
png('model_vanillas_ensembles_comparison_validation_perTF_figure1.png',units = 'in',width = 12,height = 12,res=600)
print(p3_2_1/p3_2_3)
dev.off()


### Get list of TFs well predicted in validation and train-----------------
df_corr_train <- rbind(df_corr_train,df_corr_train_random %>% select(all_of(colnames(df_corr_train))))
df_corr_train <- df_corr_train %>% group_by(cell,model) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
df_corr_train <- df_corr_train %>% mutate(tf_rank=100 * tf_rank/length(unique(df_corr_train$TF)))
df_val <- df_corr_val

#df_corr_train_ensemble_filtered <- df_corr_train %>% filter(tf_rank<=15)
#tfs_to_keep <- df_corr_train_ensemble_filtered %>% select(cell,TF) %>% unique()
df_corr_val_ensemble_filtered <- left_join(df_val,df_corr_train %>% select(model,cell,TF,tf_rank))
df_corr_val_ensemble_filtered <- df_corr_val_ensemble_filtered %>% 
  mutate(model=ifelse(grepl('shuffleY',model),'shuffle Y',ifelse(grepl('shuffleX',model),'shuffle X',model))) %>%
  mutate(model = factor(model,levels= c('DT-LEMBAS','ANN','GCNN','SVM','KNN','shuffle X','shuffle Y')))
df_corr_val_ensemble_filtered <- distinct(df_corr_val_ensemble_filtered)
df_corr_val_ensemble_filtered <- df_corr_val_ensemble_filtered %>% mutate(noisy = ifelse(tf_rank<=25,'well-fitted',
                                                                                         ifelse(tf_rank<75,'NA',
                                                                                                'poorly-fitted')))
df_corr_val_ensemble_filtered <- df_corr_val_ensemble_filtered %>% mutate(noisy=ifelse(model=='DT-LEMBAS',noisy,'NA'))
df_corr_val_ensemble_filtered$noisy <- factor(df_corr_val_ensemble_filtered$noisy,
                                              levels=c('well-fitted','poorly-fitted','NA')) #'medium-fitted'

p3_2_1 <- ggboxplot(df_corr_val_ensemble_filtered,
                    x = "model", y = "r" ,outlier.shape = NA)+
  #ggtitle('Ensembles` performance comparison with standard machine learning approaches') +
  #geom_jitter(data = df_corr_val_ensemble_filtered %>% filter(tf_rank<=25 | tf_rank>=50),
  #            aes(color=noisy),alpha=0.2)+
  geom_point(aes(color=noisy,alpha=noisy), position = position_jitter(width = 0.2))+
  scale_color_manual(labels = c('well-fitted','poorly-fitted',''),
                     values =  c("#00AFBB","#FC4E07", "black"))+
  scale_alpha_manual(values=c(1,1,0.5))+
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA))),
         alpha = "none") +
  theme(text = element_text(family='Arial',size=22),plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=22)) +
  xlab('model') + ylab('per TF pearson`s r')+ 
  scale_y_continuous(breaks = c(-0.6,-0.4,-0.2,0,seq(0.2,1,0.2)),limits = c(-0.635,1.1)) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_1 <- p3_2_1 + stat_compare_means(comparisons = list(c('DT-LEMBAS','ANN'),
                                                         c('DT-LEMBAS','GCNN'),
                                                         c('DT-LEMBAS','SVM'),
                                                         c('DT-LEMBAS','KNN'),
                                                         c('DT-LEMBAS','shuffle X'),
                                                         c('DT-LEMBAS','shuffle Y')),
                                      label = 'p.signif',
                                      label.y = c(0.5,0.6,0.7,0.8,0.9,1),
                                      method = 'wilcox.test',
                                      tip.length=0.05)
print(p3_2_1)
ggsave('../figures/figure1C.eps',
       plot = p3_2_1,
       device = cairo_ps,
       width = 12,
       height = 4,
       units = 'in',
       dpi = 600)
### Boxplot the p-values of individual TFs
plot_df <- left_join(df_pvalue_val,
                     df_corr_val %>% filter(model=='DT-LEMBAS') %>% select(TF,cell,r,mean_pear) %>% unique())
plot_df <- plot_df %>% mutate(logPadj = -log10(p.adj))
# Function to create a custom ordered plot for a single cell
# Create a custom function to create and plot a facet for each cell
create_custom_ordered_plot <- function(data, cell_line) {
  data <- data %>%
    filter(cell == cell_line) %>%
    arrange(desc(logPadj))  # Sort by logPadj in descending order
  data$TF <- factor(data$TF, levels = data[order(data$cell, -data$logPadj), ]$TF)
  
  
  p <- ggplot(data,
              aes(x = TF, y = logPadj)) +
    geom_point(aes(color = r), size = 1) +
    scale_color_gradient2(low = "blue",
                          mid = 'white',
                          high = "red",
                          midpoint = 0,
                          limits = c(-0.65, 0.8)) +
    ylab('-log10(p.adjusted)') +
    xlab('transcription factors') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'black', size = 1) +
    annotate('text',x=65,y=1.8,label='p.adjusted = 0.05',size=5)+
    theme_pubr(base_family = 'Arial', base_size = 20) +
    theme(text = element_text(size = 20, family = 'Arial'),
          axis.text.x = element_blank(),  # Rotate x-axis labels
          legend.position = 'right') +
    facet_wrap(~cell, scales = "free_x")
  
  #print(p)
}
# Create a list of custom-ordered plots for each cell
unique_cells <- unique(plot_df$cell)
plots_list <- lapply(unique_cells, function(cell) {
  create_custom_ordered_plot(plot_df, cell)
})
# Combine individual legends into one common legend
combined_plot <- wrap_plots(plots_list) + #plot_annotation(theme = theme(legend.position = 'top')) +
  plot_layout(guides = "collect") 

print(combined_plot)
ggsave('../article_supplementary_info/ajusted_pvals_perTF.png',
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)

lembas_noisy <- left_join(lembas,df_corr_train %>% filter(model=='DT-LEMBAS') %>%select(cell,TF,tf_rank) %>%unique())
# lembas_noisy <- lembas_noisy %>% filter(!(cell %in% c('MCF7','HEPG2')))
lembas_noisy <- lembas_noisy %>% mutate(noisy = ifelse(tf_rank<=25,'well-fitted',
                                                       ifelse(tf_rank<75,'NA',
                                                              'poorly-fitted')))
lembas_noisy$noisy <- factor(lembas_noisy$noisy,
                             levels=c('well-fitted','poorly-fitted','NA')) #'medium-fitted'
median_cell_values <- aggregate(r ~ cell, lembas_noisy, mean)
lembas_noisy$cell <- factor(lembas_noisy$cell, levels = median_cell_values$cell[order(median_cell_values$r)])
p3_2_3 <- ggboxplot(lembas_noisy, x = "cell", y = "r",outlier.shape = NA)+
  geom_point(aes(color=noisy,alpha=noisy), position = position_jitter(width = 0.2))+
  scale_color_manual(labels = c('well-fitted','poorly-fitted',''),
                     values =  c("#00AFBB","#FC4E07", "black"))+
  scale_alpha_manual(values=c(1,1,0.5))+
  guides(color = guide_legend(override.aes = list(shape = c(16, 16, NA))),
         alpha = "none") +
  #ggtitle('DT-LEMBAS model per cell line') + #geom_jitter(aes(color=cell))+
  theme(text = element_text(family='Arial',size=22),plot.title = element_text(size=20,hjust = 0.5,vjust=1),
        legend.position = 'top',legend.title = element_blank()) +
  xlab('cell line') + ylab('per TF pearson`s r') + scale_y_continuous(n.breaks = 10,expand = c(0, 0.2)) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_3 <- p3_2_3 + stat_compare_means(method='kruskal.test',size = 7,label.y = 0.78)
stat.test <- lembas_noisy %>% 
  rstatix::wilcox_test(r ~ cell, comparisons = list(c('A375','HA1E')))
p3_2_3 <- p3_2_3  + stat_pvalue_manual(stat.test, label = "Wilcox test p = {p}",y.position = 0.75,size = 7)
print(p3_2_3)
ggsave('../figures/figure1B.eps',
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 4,
       units = "in",
       dpi = 600)

#visualize only well-fitted
p3_2_1 <- ggboxplot(df_corr_val_ensemble_filtered  %>% filter(tf_rank<=10),
                    x = "model", y = "r" ,outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2))+
  theme(text = element_text(family = 'Arial',size=20),
        legend.position = 'right',
        legend.text = element_text(family = 'Arial',size = 20),
        legend.title = element_text(family = 'Arial',size = 20),
        legend.title.align=0.5,
        plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('per TFs pearson`s r')+ 
  scale_y_continuous(breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,seq(0.1,1,0.1)),limits = c(-0.45,1)) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_1 <- p3_2_1 + stat_compare_means(comparisons = list(c('DT-LEMBAS','ANN'),
                                                         c('DT-LEMBAS','GCNN'),
                                                         c('DT-LEMBAS','SVM'),
                                                         c('DT-LEMBAS','KNN'),
                                                         c('DT-LEMBAS','shuffle X'),
                                                         c('DT-LEMBAS','shuffle Y')),
                                      label = 'p.format',
                                      label.y = c(0.65,0.7,0.75,0.8,0.85,0.9,0.95),
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=6)
print(p3_2_1)
ggsave('../figures/figure1D.eps',
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 8,
       units = "in",
       dpi = 600)
