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

### Load data
files <- list.files('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/')
files <- files[grep('.csv',files)]

# # Load train correlation data
# all_corr_train <- data.frame()
# files_corr <- files[grep('train',files)]
# i <- 1
# for (file in files_corr){
#   cell <- str_split_fixed(str_split_fixed(file,'_',4)[1,4],'.csv',2)[1,1]
#   ensembles_no <- str_split_fixed(str_split_fixed(file,'_',4)[1,2],'ensembles',2)[1,2]
#   combo_no <- str_split_fixed(str_split_fixed(file,'_',4)[1,3],'combo',2)[1,2]
#   file <- paste0('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/',file)
#   tmp <- data.table::fread(file,header=T)
#   colnames(tmp) <- c('TF','r')
#   tmp <- tmp %>% mutate(ensembles=ensembles_no)
#   tmp <- tmp %>% mutate(cell=cell)
#   tmp <- tmp %>% mutate(combo_id = combo_no)
#   tmp <- tmp %>% select(ensembles,combo_id,cell,TF,r)
#   all_corr_train <- rbind(all_corr_train,tmp)
#   if (i %% 500 == 0){
#     print(ensembles_no)
#   }
#   i <- i+1
# }
# # Load validation correlation data
# all_corr_val <- data.frame()
# files_corr <- files[grep('valEnsemblePerformance',files)]
# i <- 1
# for (file in files_corr){
#   cell <- str_split_fixed(str_split_fixed(file,'_',4)[1,4],'.csv',2)[1,1]
#   ensembles_no <- str_split_fixed(str_split_fixed(file,'_',4)[1,2],'ensembles',2)[1,2]
#   combo_no <- str_split_fixed(str_split_fixed(file,'_',4)[1,3],'combo',2)[1,2]
#   file <- paste0('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/',file)
#   tmp <- data.table::fread(file,header=T)
#   colnames(tmp) <- c('TF','r')
#   tmp <- tmp %>% mutate(ensembles=ensembles_no)
#   tmp <- tmp %>% mutate(cell=cell)
#   tmp <- tmp %>% mutate(combo_id = combo_no)
#   tmp <- tmp %>% select(ensembles,combo_id,cell,TF,r)
#   all_corr_val <- rbind(all_corr_val,tmp)
#   if (i %% 500 == 0){
#     print(ensembles_no)
#   }
#   i <- i+1
# }
# saveRDS(all_corr_val,
#         'Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/figures/all_corr_val.rds')
# saveRDS(all_corr_train,
#         'Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/figures/all_corr_train.rds')
all_corr_val <- data.table::fread('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/valEnsemblePerformance_ensembles_all.csv')
colnames(all_corr_val)[1] <- 'TF' 
all_corr_train <- data.table::fread('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/trainEnsemblePerformance_ensembles_all.csv')
colnames(all_corr_train)[1] <- 'TF' 
### Get average across all TFs and create data frame for visualization
all_corr_val <- all_corr_val %>% group_by(ensembles,trial,cell) %>% mutate(mean_across_tfs = mean(r,na.rm = T)) %>% ungroup()
all_corr_val <- all_corr_val %>% group_by(ensembles,cell) %>% mutate(mean_r = mean(mean_across_tfs,na.rm = T)) %>% 
  mutate(std_r = sd(mean_across_tfs,na.rm = T)) %>% ungroup()
all_corr_train <- all_corr_train %>% group_by(ensembles,trial,cell) %>% mutate(mean_across_tfs = mean(r,na.rm = T)) %>% ungroup()
all_corr_train <- all_corr_train %>% group_by(ensembles,cell) %>% mutate(mean_r = mean(mean_across_tfs,na.rm = T)) %>% 
  mutate(std_r = sd(mean_across_tfs,na.rm = T)) %>% ungroup()
all_corr_val <- all_corr_val %>% mutate(ensembles=factor(ensembles,levels=seq(1,51)))
all_corr_val$ensembles <- as.numeric(all_corr_val$ensembles)
all_corr_train <- all_corr_train %>% mutate(ensembles=factor(ensembles,levels=seq(1,51)))
all_corr_train$ensembles <- as.numeric(all_corr_train$ensembles)
gc()

# remove some information
all_corr_train <- distinct(all_corr_train %>% select(ensembles,cell,mean_r,std_r))
all_corr_val <- distinct(all_corr_val %>% select(ensembles,cell,mean_r,std_r))
gc()

p <- ggplot(all_corr_val,aes(x=ensembles,y=mean_r,color = cell)) + geom_point() +
  geom_line() + geom_errorbar(aes(ymin= mean_r-2*std_r , ymax=mean_r+2*std_r),width = 0.1)
p <- p + ggtitle('Performance for different number of ensembles') + 
  xlab('number of ensembles') + ylab('average pearson`s r across TFs')
p <- p + xlim(c(0,max(all_corr_val$ensembles))) + ylim(c(0,max(all_corr_val$mean_r)+0.1))
p <- p + theme_minimal(base_family = "Arial",base_size = 16)+theme(plot.title = element_text(hjust = 0.5,size=16))
print(p)

png('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/figures/validationPerformance_vs_numberOfEnsembles.png',
    width=8,height=8,units = "in",res = 600)
print(p)
dev.off()

p2 <- ggplot(all_corr_train,aes(x=ensembles,y=mean_r,color = cell)) + geom_point() +
  geom_line() + geom_errorbar(aes(ymin= mean_r-2*std_r , ymax=mean_r+2*std_r),width = 0.1)
p2 <- p2 + ggtitle('Training performance for different number of ensembles') + 
  xlab('number of ensembles') + ylab('average pearson`s r across TFs')
p2 <- p2 + xlim(c(0,max(all_corr_train$ensembles))) + ylim(c(0,max(all_corr_train$mean_r)+0.1))
p2 <- p2 + theme_minimal(base_family = "Arial",base_size = 16)+theme(plot.title = element_text(hjust = 0.5,size=16))
print(p2)

png('Model/CVL1000_Paper/FinalEnsemble/EstimateEnsembles/figures/trainPerformance_vs_numberOfEnsembles.png',
    width=8,height=8,units = "in",res = 600)
print(p2)
dev.off()
