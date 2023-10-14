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

### Load data per TF analysis-----------
files <- list.files('../results/FinalEnsemble/test/')
files <- files[grep('.csv',files)]

# Load train correlation data
all_corr_train <- data.frame()
files_corr <- files[grep('trainPerformance',files)]
for (file in files_corr){
  cell <- str_split_fixed(file,'_',4)[1,2]
  file <- paste0('../results/FinalEnsemble/test/',file)
  tmp <- data.table::fread(file)
  colnames(tmp)[1] <- 'no'
  tmp <- tmp %>% select(- no) %>% unique()
  tmp <- tmp %>% mutate(model='Individual models')
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  all_corr_train <- rbind(all_corr_train,tmp)
}
# Load validation correlation data
all_corr_val <- data.frame()
files_corr <- files[grep('valPerformance',files)]
for (file in files_corr){
  cell <- str_split_fixed(file,'_',4)[1,2]
  file <- paste0('../results/FinalEnsemble/test/',file)
  tmp <- data.table::fread(file)
  colnames(tmp)[1] <- 'no'
  tmp <- tmp %>% select(- no) %>% unique()
  tmp <- tmp %>% mutate(model='Individual models')
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  all_corr_val <- rbind(all_corr_val,tmp)
}

# Add ensemble performance
# Load train correlation data
df_corr_train_ensemble <- data.frame()
files_corr <- files[grep('trainEnsemblePerformance',files)]
for (file in files_corr){
  cell <- str_split_fixed(file,'_',2)[1,2]
  cell <- str_split_fixed(cell,'[.]',2)[1,1]
  file <- paste0('../results/FinalEnsemble/test/',file)
  tmp <- data.table::fread(file,header = T)
  colnames(tmp)[1] <- 'TF'
  colnames(tmp)[2] <- 'r'
  tmp <- tmp %>% mutate(model='Ensembles')
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  df_corr_train_ensemble <- rbind(df_corr_train_ensemble,tmp)
}
# Load validation correlation data
df_corr_val_ensemble <- data.frame()
files_corr <- files[grep('valEnsemblePerformance',files)]
for (file in files_corr){
  cell <- str_split_fixed(file,'_',2)[1,2]
  cell <- str_split_fixed(cell,'[.]',2)[1,1]
  file <- paste0('../results/FinalEnsemble/test/',file)
  tmp <- data.table::fread(file,header = T)
  colnames(tmp)[1] <- 'TF'
  colnames(tmp)[2] <- 'r'
  tmp <- tmp %>% mutate(model='Ensembles')
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  df_corr_val_ensemble <- rbind(df_corr_val_ensemble,tmp)
}

# Visualize
df_corr_train <- all_corr_train %>% gather('TF','r',-model,-cell)
df_corr_train_ensemble <-  df_corr_train_ensemble %>% select(model,cell,TF,r)
df_corr_train <- rbind(df_corr_train,df_corr_train_ensemble)
df_corr_train <- df_corr_train %>% group_by(model,cell) %>%  mutate(mean_r=mean(r)) %>% ungroup()
df_corr_train$model <- factor(df_corr_train$model,
                              levels = c('Individual models','Ensembles'))
p3_1_1 <- ggboxplot(df_corr_train, x = "model", y = "r",add='jitter')+
  ggtitle('Train performance of individual runs and ensemble predictions') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('per TF pearson`s r')+ scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_1_1 <- p3_1_1 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
p3_1_2 <- ggboxplot(df_corr_train %>% select(model,mean_r) %>% unique(), x = "model", y = "mean_r",add='jitter')+
  ggtitle('') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('average per TF pearson`s r') + scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_1_2 <- p3_1_2 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
print(p3_1_1/p3_1_2)

png('../article_supplementary_info/suppl_figure_2A_training.png',units = 'in',width = 12,height = 12,res=600)
print(p3_1_1/p3_1_2)
dev.off()

#validation visualize
df_corr_val <- all_corr_val %>% gather('TF','r',-model,-cell)
df_corr_val_ensemble <-  df_corr_val_ensemble %>% select(model,cell,TF,r)
df_corr_val <- rbind(df_corr_val,df_corr_val_ensemble)
df_corr_val <- df_corr_val %>% group_by(model,cell) %>%  mutate(mean_r=mean(r)) %>% ungroup()
df_corr_val$model <- factor(df_corr_val$model,
                              levels = c('Individual models','Ensembles'))
p3_2_1 <- ggboxplot(df_corr_val, x = "model", y = "r",add='jitter')+
  ggtitle('Validation performance of individual runs and ensemble predictions') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('per TF pearson`s r')+ scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_1 <- p3_2_1 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
p3_2_2 <- ggboxplot(df_corr_val %>% select(model,mean_r) %>% unique(), x = "model", y = "mean_r",add='jitter')+
  ggtitle('') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('average per TF pearson`s r') + scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_2 <- p3_2_2 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
print(p3_2_1/p3_2_2)

png('../article_supplementary_info/suppl_figure2A.png',units = 'in',width = 12,height = 12,res=600)
print(p3_2_1/p3_2_2)
dev.off()

data.table::fwrite(df_corr_val_ensemble,'../results/FinalEnsemble/perTFValidationPerformance.csv')

### Check relationship between goodness of fit and validation performance---------------------------
### Find median rank cut-off to include versus train correlation
df_corr_train_ensemble <- df_corr_train_ensemble %>% group_by(cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
df_corr_train_ensemble <- df_corr_train_ensemble %>% mutate(tf_rank=100 * tf_rank/length(unique(df_corr_train_ensemble$TF)))

p6 <- ggplot(df_corr_train_ensemble,aes(x=tf_rank,y=r)) + geom_point() + geom_smooth(color='black') +
  xlab('% median rank based on performace in train') + ylab('per TF average correlation in train')+ 
  theme_pubr() +
  ggtitle('Ranking of best fitted TFs') + theme(text = element_text(size=13),
                                                plot.title = element_text(hjust = 0.5),
                                                panel.grid.major = element_line(colour="black",
                                                                                linetype = 'dashed', 
                                                                                linewidth=0.25)) +
  scale_y_continuous(breaks = seq(-0.75,1,0.25))
print(p6)

## visualize
df_corr <- left_join(df_corr_val_ensemble,df_corr_train_ensemble %>% select(-r),by=c('model','cell','TF'))
df_corr <- df_corr %>% group_by(cell,TF) %>%  
  mutate(mean_rank = mean(tf_rank)) %>% 
  mutate(mean_r = mean(r)) %>% ungroup()
p7 <- ggplot(df_corr,aes(x=tf_rank,y=r,color=cell)) + geom_point() + geom_smooth(color='black') +
  xlab('% median rank based on performace in training set') + ylab('per TF correlation in validation')+ 
  theme_pubr() +
  ggtitle('Ranking of best fitted TFs VS validation performance') + theme(text = element_text(size=22),
                                                                          plot.title = element_text(hjust = 0.5),
                                                                          panel.grid.major = element_line(colour="black",
                                                                                                          linetype = 'dashed', 
                                                                                                          linewidth=0.25)) +
  scale_y_continuous(limits = c(-0.85,0.85),breaks = seq(-0.85,0.85,0.2))
png('../article_supplementary_info/suppl_figure6B.png',
    units = 'in',width = 12,height = 12,res = 600)
print(p7)
dev.off()

cutPear <- data.frame()
thresholds <- c(1,2.5,5,10,15,20,25,30,40,50,70,90,100)
for (i in 1:length(thresholds)){
  tmp <- df_corr %>% filter(tf_rank<=thresholds[i])
  tmp <- tmp %>% group_by(model,cell) %>% mutate(r= mean(r)) %>% ungroup()
  tmp <- tmp %>% select(model,cell,r) %>% unique() %>% mutate(thresh=thresholds[i])
  cutPear <- rbind(cutPear,tmp)
}
p8 <- ggplot(cutPear,aes(x=thresh,y=r,color=cell)) + geom_point() + geom_smooth(se=F,color='black') +
  xlab('% median rank based on performace in training set') + ylab('average correlation in validation')+ 
  theme_pubr() +
  ggtitle('Ranking of best fitted TFs VS validation performance') + theme(text = element_text(size=22),
                                                                          plot.title = element_text(hjust = 0.5),
                                                                          panel.grid.major = element_line(colour="black",
                                                                                                          linetype = 'dashed', 
                                                                                                          linewidth=0.25)) +
  scale_y_continuous(limits = c(0,0.8),breaks = seq(0,0.8,0.1))
print(p8)
png('../article_supplementary_info/suppl_figure6A.png'
    ,units = 'in',width = 12,height = 12,res = 600)
print(p8)
dev.off()

### Visualise only well-fitted TFs-------------------------
# Load train correlation data
all_corr_train <- data.frame()
files_corr <- files[grep('trainPerformance',files)]
for (file in files_corr){
  cell <- str_split_fixed(file,'_',4)[1,2]
  file <- paste0('../results/FinalEnsemble/test/',file)
  tmp <- data.table::fread(file)
  colnames(tmp)[1] <- 'no'
  tmp <- tmp %>% unique()
  tmp <- tmp %>% mutate(model='Individual models')
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  all_corr_train <- rbind(all_corr_train,tmp)
}
# Load validation correlation data
all_corr_val <- data.frame()
files_corr <- files[grep('valPerformance',files)]
for (file in files_corr){
  cell <- str_split_fixed(file,'_',4)[1,2]
  file <- paste0('../results/FinalEnsemble/test/',file)
  tmp <- data.table::fread(file)
  colnames(tmp)[1] <- 'no'
  tmp <- tmp %>% unique()
  tmp <- tmp %>% mutate(model='Individual models')
  tmp <- tmp %>% mutate(cell=cell)
  names <- colnames(tmp)
  all_corr_val <- rbind(all_corr_val,tmp)
}
df_corr_val <- all_corr_val %>% gather('TF','r',-model,-cell,-no)
df_corr_train <- all_corr_train %>% gather('TF','r',-model,-cell,-no)
df_corr_train <- df_corr_train %>% group_by(cell,no) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
df_corr_train <- df_corr_train %>% mutate(tf_rank=100 * tf_rank/length(unique(df_corr_train$TF)))
df_corr_train <- df_corr_train %>% filter(tf_rank<=10)
# saveRDS(tfs_to_keep,'Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/tfs_to_keep_10perc.rds')
tfs_to_keep <- df_corr_train %>% select(no,cell,TF) %>% unique()
df_corr_val <- left_join(tfs_to_keep,df_corr_val)
df_corr_val <- df_corr_val %>% select(-no) %>% unique()
df_corr_train <- df_corr_train %>% select(-no) %>% unique()

df_corr_train_ensemble_filtered <- df_corr_train_ensemble %>% filter(tf_rank<=10)
tfs_to_keep <- df_corr_train_ensemble_filtered %>% select(cell,TF) %>% unique()
df_corr_val_ensemble_filtered <- left_join(tfs_to_keep,df_corr_val_ensemble)

df_corr_train <- rbind(df_corr_train_ensemble_filtered,df_corr_train)
df_corr_val <- rbind(df_corr_val_ensemble_filtered,df_corr_val)

df_corr_train <- df_corr_train %>% group_by(model,cell) %>%  mutate(mean_r=mean(r)) %>% ungroup()
df_corr_train$model <- factor(df_corr_train$model,
                              levels = c('Individual models','Ensembles'))
df_corr_val <- df_corr_val %>% group_by(model,cell) %>%  mutate(mean_r=mean(r)) %>% ungroup()
df_corr_val$model <- factor(df_corr_val$model,
                              levels = c('Individual models','Ensembles'))

p3_1_1 <- ggboxplot(df_corr_train, x = "model", y = "r",add='jitter')+
    ggtitle('Individual runs and ensemble predictions in top 10% best-fitted TFs') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('per TFs pearson`s r')+ scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_1_1 <- p3_1_1 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
p3_1_2 <- ggboxplot(df_corr_train %>% select(model,mean_r) %>% unique(), x = "model", y = "mean_r",add='jitter')+
  ggtitle('') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('average pearson`s r') + scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_1_2 <- p3_1_2 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
print(p3_1_1/p3_1_2)
png('../article_supplementary_info/suppl_figure_2B_training.png',units = 'in',width = 12,height = 12,res=600)
print(p3_1_1/p3_1_2)
dev.off()

p3_2_1 <- ggboxplot(df_corr_val, x = "model", y = "r",add='jitter')+
  ggtitle('Individual runs and ensemble predictions in top 10% best-fitted TFs') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('per TFs pearson`s r')+ scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_1 <- p3_2_1 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
p3_2_2 <- ggboxplot(df_corr_val %>% select(model,mean_r) %>% unique(), x = "model", y = "mean_r",add='jitter')+
  ggtitle('') +
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5)) +
  xlab('model') + ylab('average pearson`s r') + scale_y_continuous(n.breaks = 10) +
  geom_hline(yintercept = 0,linetype='dashed',linewidth=1,color='black')
p3_2_2 <- p3_2_2 + stat_compare_means(comparisons = list(c('Individual models','Ensembles')),
                                      label = 'p.format',
                                      method = 'wilcox.test',
                                      tip.length=0.05,
                                      size=5)
print(p3_2_1/p3_2_2)

png('../article_supplementary_info/suppl_figure_2B.png',units = 'in',width = 12,height = 12,res=600)
print(p3_2_1/p3_2_2)
dev.off()
