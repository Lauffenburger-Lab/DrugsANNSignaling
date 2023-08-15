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
# library(doFuture)
# # parallel set number of workers
# registerDoFuture()
# plan(multiprocess,workers = 9)

# Load performance results-----------------
performance <- read.xlsx('Model/CVL1000_Paper/HyperParamInteractionResults/hyperParamPerformanceResults.xlsx',
                         sheetIndex = 1)
performance <- performance %>% filter(!is.na(corr))
performance <- performance %>% group_by(lambda) %>% mutate(mean_corr=mean(corr)) %>%
  mutate(sd_corr=sd(corr)) %>% ungroup()
#performance <- performance[-c(33,34,35,36,37,38,39,40),]

# Load performance at 0 and inf regularization
performance_zero <- read.xlsx('Model/CVL1000_Paper/MultipleRuns/modeltype4.xlsx',sheetIndex=1)
performance_zero <- performance_zero %>% mutate(lambda=0) %>% group_by(lambda) %>% mutate(mean_corr=mean(corr)) %>%
  mutate(sd_corr=sd(corr)) %>% ungroup()
performance_inf <-  read.xlsx('Model/CVL1000_Paper/MultipleRuns/modeltype1.xlsx',sheetIndex=1)
performance_inf <- performance_inf %>% filter(!is.na(corr))
performance_inf <- performance_inf %>% mutate(lambda=Inf) %>% group_by(lambda) %>% mutate(mean_corr=mean(corr)) %>%
  mutate(sd_corr=sd(corr)) %>% ungroup()

performance <- rbind(performance_zero,performance,performance_inf)
# df <-
#   pairwise_comparisons(performance, lambda, corr,type ="robust" ,var.equal = F) %>%
#   dplyr::mutate(groups = purrr::pmap(.l = list(group1, group2), .f = c)) %>%
#   dplyr::arrange(group1)

### Read random performance results
shuffleX <- read.xlsx('Model/CVL1000_Paper/RandomSuffleX/test/randomShufflePerformance.xlsx',sheetIndex=1)
shuffleX <- shuffleX %>% select(-no) %>% mutate(lambda = 'shuffle X')
shuffleX <-shuffleX %>% group_by(lambda) %>% mutate(mean_corr=mean(corr)) %>%
mutate(sd_corr=sd(corr)) %>% ungroup()
shuffleY <- read.xlsx('Model/CVL1000_Paper/RandomSuffleY/test/randomShufflePerformance.xlsx',sheetIndex=1)
shuffleY <- shuffleY %>% select(-no) %>% mutate(lambda = 'shuffle Y')
shuffleY <-shuffleY %>% group_by(lambda) %>% mutate(mean_corr=mean(corr)) %>%
  mutate(sd_corr=sd(corr)) %>% ungroup()

performance <- rbind(performance,shuffleX,shuffleY)


### Plot
N <- length(unique(performance$cell))
p1 <- ggboxplot(performance, x = "lambda", y = "corr",add='jitter')+
 theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  ylim(c(0.0,1)) + xlab('λ')
p1 <- p1 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                 c('0','shuffle X'),
                                                 c('0','shuffle Y')),
                              label = 'p.format',
                              label.y =  c(0.8,0.85,0.9),
                              method = 'wilcox.test',
                              tip.length=0.05) 
p1 <- p1 + stat_compare_means(method = "kruskal.test", label.y = 0.97,label.x = 1, 
                              data = performance %>% filter(!(lambda %in% c('shuffle X','shuffle Y'))) %>%
                                mutate(lambda = factor(lambda,levels = c(0,1e-8,1e-7,1e-6,
                                                                         1e-5,1e-4,0.001,
                                                                         0.005,0.01,0.05,'Inf'))),
                              aes(x= lambda, y = corr, color = lambda))
print(p1)


# Load regularization results
regularization <- read.xlsx('Model/CVL1000_Paper/HyperParamInteractionResults/hitDiscoveryResults.xlsx',
                            sheetIndex = 1)[,1:13]
#regularization <- regularization %>% filter(threshold>=5 & threshold<=50)
#regularization <- regularization %>% filter(lamda>=1e-7 & lamda<=0.001 & lamda!=0.001)
regularization$threshold <- factor(regularization$threshold,
                                   levels = unique(regularization$threshold))
# p2 <- ggplot(regularization,aes(x=lamda,y=new_hits)) +
#   geom_point(aes(color=threshold),size=2)+
#   #geom_smooth(aes(color=threshold),se = F) +
#   geom_line(aes(color=threshold)) +
#   scale_x_log10(breaks=c(1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1),
#                                 labels=c(1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1))+
#   # scale_y_log10()+
#   ylab('New drug-target interactions')+xlab('λ')+
#   ggtitle('New discovery rate as a function of regularization')+
#   theme_minimal()+
#   theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
# print(p2)


#regularization$new_hits <- regularization$new_hits +1
#ggboxplot(regularization, x = "lamda", y = "G",add='jitter')
p2 <- ggplot(regularization %>% group_by(lamda) %>% mutate(meanG=mean(G)) %>%
               mutate(sdG=sd(G)/sqrt(n_distinct(G))) %>% ungroup(),
             aes(x=lamda,y=meanG))+
  geom_point() + geom_smooth(color='black')+ scale_x_log10()+
  geom_errorbar(aes(ymin=meanG-sdG, ymax=meanG+sdG), width=.3)+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))  + 
  xlab('λ') + ylab('G-mean') + 
  theme(legend.position = "none") + ylim(c(0,1))+  theme_pubr() + 
  (ggplot(regularization,aes(x=lamda,optimalG)) + geom_smooth(color='black',se=F) + geom_point()+
     ylim(c(0,1)) + scale_x_log10() + 
     xlab('λ') + ylab('optimal G-mean') +  theme_pubr()) +
  plot_annotation(title = "Trade-off of sensitivity/specificity for different levels of regularization") &
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2)

png('regularizationTuning.png',units = 'in',width = 12,height = 12,res = 600)
p1/p2
dev.off()


# Compare per TF performance for different levels of regularization---------------
files <- list.files('Model/CVL1000_Paper/MultipleRuns/test/')
files <- files[grep('.csv',files)]
files_inf <- files[grep('modeltype1',files)]
files_inf <- files_inf[which(!grepl('class',files_inf))]
files_inf <- files_inf[which(!grepl('quantile',files_inf))]
files_zero <- files[grep('modeltype4',files)]
files_zero <- files_zero[which(!grepl('class',files_zero))]
files_zero <- files_zero[which(!grepl('quantile',files_zero))]
files_zero <- files_zero[which(!grepl('R2',files_zero))]

files_reg <-  list.files('Model/CVL1000_Paper/hyperParamTune/test/')
files_reg <- files_reg[grep('.csv',files_reg)]
files_reg <- files_reg[which(!grepl('class',files_reg))]
files_reg <- files_reg[which(!grepl('quantile',files_reg))]

all_files  <- list(files_inf,files_zero,files_reg)

# Load variance data
all_vars <- data.frame()
for (i in 1:length(all_files)){
  files_var <- all_files[[i]][grep('Variance',all_files[[i]])]
  for (file in files_var){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_vars <- rbind(all_vars,tmp)
  }
}
df_var <- all_vars %>% gather('TF','std',-lamda,-cell)
# #df_var <- df_var %>% group_by(lamda,TF) %>% mutate(mean_std = mean(std))
# p3 <- ggboxplot(df_var, x = "lamda", y = "std",add='jitter')+
#   theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
#   xlab('λ')
# 
# p3 <- p3 + stat_compare_means(method = "kruskal.test")
# print(p3)

# Load train correlation data
all_corr_train <- data.frame()
for (i in 1:length(all_files)){
  files_corr <- all_files[[i]][grep('trainPerformance',all_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_train <- rbind(all_corr_train,tmp)
  }
}
df_corr_train <- all_corr_train %>% gather('TF','r',-lamda,-cell)
# #### to add ranking ###
# df_corr_train <- df_corr_train %>% group_by(cell,TF) %>% mutate(r=ifelse(lamda==0 | lamda==Inf,mean(r),r)) %>% ungroup() %>% unique()
# df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
# df_corr_train <- df_corr_train %>% filter(tf_rank<=15)
# tfs_to_keep <- df_corr_train %>% select(lamda,cell,TF) %>% unique()
# colnames(tfs_to_keep)[3] <- 'TF_train'
# ### Filtering rank end ###
df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_train <- df_corr_train %>% select(-TF,-r) %>% unique()
p3 <- ggboxplot(df_corr_train, x = "lamda", y = "mean_pear",add='jitter')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('per TF pearson`s r') + scale_y_continuous(limits = c(0,1))
#p3 <- p3 + stat_compare_means(comparisons = list(c('0','Inf')),
#                              label = 'p.format',
#                              method = 'wilcox.test',
#                              tip.length=0.05) 
##p3 <- p3 + stat_compare_means(method = "kruskal.test")
print(p3)

# Load validation correlation data
all_corr_val <- data.frame()
for (i in 1:length(all_files)){
  files_corr <- all_files[[i]][grep('valPerformance',all_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_val <- rbind(all_corr_val,tmp)
  }
}
df_corr_val <- all_corr_val %>% gather('TF','r',-lamda,-cell)
df_corr_val <- df_corr_val %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
# ### filter based on training ###
# df_corr_val <- left_join(df_corr_val,tfs_to_keep)
# df_corr_val <- df_corr_val %>% mutate(keep= (TF==TF_train)) %>% filter(keep==TRUE) %>% 
#   select(-TF_train,-keep) %>% unique()
# ####filtering ended ###
df_corr_val <- df_corr_val %>% select(-TF,-r) %>% unique()
p3 <- ggboxplot(df_corr_val, x = "lamda", y = "mean_pear",add='jitter')+
  ggtitle('Average per TF pearson`s correlation across all TFs in diffent cell-lines for increasing level of regularization')+
  xlab(expression(lambda)) + ylab('average per TF pearson`s r') + 
  theme(text = element_text(family = 'Arial',size=20),plot.title = element_text(hjust = 0.5,size=17)) +
  scale_y_continuous(limits = c(0,0.5))
p3 <- p3 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                 c('0','0.005'),
                                                 c('0.005','0.01'),
                                                 c('0.005','0.001')),
                              label = 'p.format',
                              method = 'wilcox.test',
                              tip.length=0.05,
                              size=7) 
##p3 <- p3 + stat_compare_means(method = "kruskal.test")
print(p3)

## redo that for figure1 with line plot and boxplot
p3_1 <- ggplot(df_corr_val %>% filter(!(lamda=='0' | lamda=='Inf')) %>% group_by(lamda) %>% 
                 mutate(mu = mean(mean_pear)) %>% mutate(se = sd(mean_pear)/sqrt(length(unique(df_corr_val$cell)))) %>%
                 ungroup() %>% mutate(lamda=as.numeric(lamda)),
               aes(x=lamda,y=mu)) +
  geom_point(color='black') + geom_line(linewidth=0.75) + geom_errorbar(aes(ymin=mu-se, ymax=mu+se), width=.3)+
  scale_x_log10()+
  ggtitle('Average per TF pearson`s correlation across all TFs across all cell-lines')+
  xlab(expression(lambda)) + ylab('average pearson`s r') + 
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),plot.title = element_blank(),
        axis.title.x =  element_text(family = 'Arial',size=36),
        plot.margin = margin(1,10,1,1, "mm")) +
  scale_y_continuous(limits = c(0,0.35))
print(p3_1)
p3_2 <- ggboxplot(df_corr_val %>% filter(lamda=='0' | lamda=='Inf'), x = "lamda", y = "mean_pear",add='jitter')+
  ggtitle('Average per TF pearson`s correlation across all TFs in diffent cell-lines')+
  xlab(expression(lambda)) + ylab('average pearson`s r') + 
  theme(text = element_text(family = 'Arial',size=24),plot.title = element_blank(),
        axis.title.x =  element_text(family = 'Arial',size=36),
        axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0,0.35))
stat.test <- df_corr_val %>% filter(lamda=='0' | lamda=='Inf') %>% 
  rstatix::wilcox_test(mean_pear ~ lamda, comparisons = list(c('0','Inf')))
p3_2 <- p3_2  + stat_pvalue_manual(stat.test, label = "Wilcoxon test p = {p}",y.position = 0.3,size = 9)
p3 <- p3_1 + p3_2 + plot_annotation(title = 'Average per TF pearson`s correlation across all TFs across all cell-lines',
                                    theme = theme(plot.title =element_text(family='Arial',hjust = 0.5,size=24)))
print(p3)
# df_corr_val <- df_corr_val %>% select(-TF,-r) %>% unique()
ggsave('averagePerTFPear_vs_regularization.eps',
       plot = p3,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 8,
       units = "in",
       dpi = 600)

### Find under-fitted TFs
df_corr_all <- left_join(all_corr_train %>% gather('TF','r',-lamda,-cell),
                         all_corr_val %>% gather('TF','r',-lamda,-cell),
                         by=c('lamda','cell','TF'))
p4 <- ggscatter(df_corr_all, x = 'r.x', y ='r.y',
          rug = TRUE,
          alpha = 0.5,size=1,
          cor.coef=T,cor.coef.size = 5,cor.coef.coord = c(-0.3, 0.7)) + 
  geom_abline(slope=1,intercept=0,color='red',lty=2,linewidth=1)+
  ggtitle('Correlation between train performance and validation performance per TF') +
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  xlab('per TF pearson`s r in train') + ylab('per TF pearson`s r in validation')
print(p4)
png('performance_train_vs_validation.png',units = 'in',width = 12,height = 12,res = 600)
print(p4)
dev.off()

df <- left_join(df_var,
                all_corr_train %>% gather('TF','r',-lamda,-cell),
                by=c('lamda','cell','TF'))
#df <- df %>% group_by(lamda,cell,TF) %>% mutate(mean_r=mean(r)) %>% ungroup()
p5 <- ggscatter(df, x = 'std', y ='r',
                rug = TRUE,
                alpha = 0.5,size=1,
                cor.coef=T,cor.coef.size = 5,cor.coef.coord = c(0.23, -0.33)) + 
  geom_abline(slope=1,intercept=0,color='red',lty=2,linewidth=1)+
  ggtitle('Correlation between standard deviation and train performance per TF') +
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  xlab('per TF std in train') + ylab('per TF pearson`s r in train')
print(p5)
png('performance_train_vs_std.png',units = 'in',width = 12,height = 12,res = 600)
print(p5)
dev.off()

### Find median rank cut-off to include versus train correlation
df_corr_train <- all_corr_train %>% gather('TF','r',-lamda,-cell)
df_corr_train <- df_corr_train %>% group_by(cell,TF) %>% mutate(r=ifelse(lamda==0 | lamda==Inf,mean(r),r)) %>% ungroup() %>% unique()
df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
#df_corr_train <- df_corr_train %>% group_by(TF) %>% mutate(median_rank=median(tf_rank)) %>% ungroup()
#df_corr_train <- df_corr_train %>% group_by(TF) %>% mutate(mean_r=mean(r)) %>% ungroup()
#df_corr_train <- df_corr_train %>% mutate(median_rank=100 * median_rank/length(unique(df_corr_train$TF)))
df_corr_train <- df_corr_train %>% mutate(tf_rank=100 * tf_rank/length(unique(df_corr_train$TF)))

p6 <- ggplot(df_corr_train,aes(x=tf_rank,y=r)) + geom_point() + geom_smooth(color='black') +
  xlab('% median rank based on performace in train') + ylab('per TF average correlation in train')+ 
  theme_pubr() +
  ggtitle('Ranking of best fitted TFs') + theme(text = element_text(size=13),
                                                plot.title = element_text(hjust = 0.5),
                                                panel.grid.major = element_line(colour="black",
                                                                                linetype = 'dashed', 
                                                                                linewidth=0.25)) +
  scale_y_continuous(breaks = seq(-0.75,1,0.25))
print(p6)
df_corr_val <- all_corr_val %>% gather('TF','r',-lamda,-cell)
# df_corr_val <- df_corr_val %>% group_by(lamda,cell) %>% mutate(max_r = max(r)) %>% ungroup()
df_corr_val <- df_corr_val %>% group_by(cell,TF) %>% mutate(r=ifelse(lamda==0 | lamda==Inf,mean(r),r)) %>% ungroup() %>% unique()
# df_corr_val <- df_corr_val %>% group_by(TF) %>% mutate(mean_r=mean(r)) %>% ungroup()
# df_corr_val <- df_corr_val %>% select(lamda,cell,TF,c('mean_r_val'='mean_r')) %>% unique()

df_corr <- left_join(df_corr_val,df_corr_train %>% select(-r),by=c('lamda','cell','TF'))
df_corr <- df_corr %>% group_by(cell,TF) %>%  
  mutate(mean_rank = mean(tf_rank)) %>% 
  mutate(mean_r = mean(r)) %>% ungroup()
p7 <- ggplot(df_corr,aes(x=tf_rank,y=r,color=cell)) + geom_point() + geom_smooth(color='black') +
  xlab('% median rank based on performace in train') + ylab('per TF correlation in validation')+ 
  theme_pubr() +
  ggtitle('Ranking of best fitted TFs VS validation performance') + theme(text = element_text(size=13),
                                                plot.title = element_text(hjust = 0.5),
                                                panel.grid.major = element_line(colour="black",
                                                                                linetype = 'dashed', 
                                                                                linewidth=0.25)) +
  scale_y_continuous(limits = c(-0.8,0.8),breaks = seq(-0.8,0.8,0.2))+
  facet_wrap(~lamda)
#png('rank_vs_performance.png',units = 'in',width = 12,height = 12,res = 600)
print(p7)
#dev.off()

cutPear <- data.frame()
thresholds <- c(1,2.5,5,10,15,20,25,30,40,50,70,90,100)
for (i in 1:length(thresholds)){
  tmp <- df_corr %>% filter(tf_rank<=thresholds[i])
  tmp <- tmp %>% group_by(lamda,cell) %>% mutate(r= mean(r)) %>% ungroup()
  tmp <- tmp %>% select(lamda,cell,r) %>% unique() %>% mutate(thresh=thresholds[i])
  cutPear <- rbind(cutPear,tmp)
}
p8 <- ggplot(cutPear,aes(x=thresh,y=r,color=cell)) + geom_point() + geom_smooth(se=F,color='black') +
  xlab('% threshold') + ylab('average per TF correlation in validation')+ 
  theme_pubr() +
  ggtitle('Ranking of best fitted TFs VS validation performance') + theme(text = element_text(size=13),
                                                                          plot.title = element_text(hjust = 0.5),
                                                                          panel.grid.major = element_line(colour="black",
                                                                                                          linetype = 'dashed', 
                                                                                                          linewidth=0.25)) +
  scale_y_continuous(limits = c(-0.2,0.8),breaks = seq(-0.2,0.8,0.2))+
  facet_wrap(~lamda)
print(p8)
#png('rankThreshold_vs_performance.png',units = 'in',width = 12,height = 12,res = 600)
#print(p8)
#dev.off()

### Visualize based on filtering
# Add random shuffle
random_files <- list.files('Model/CVL1000_Paper/RandomSuffleY/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesY <- random_files[grep('Performance',random_files)]
random_filesY <- random_filesY[which(!grepl('class',random_filesY))]
random_filesY <- random_filesY[which(!grepl('quantile',random_filesY))]
random_filesY <- random_filesY[which(!grepl('R2',random_filesY))]
random_files <- list.files('Model/CVL1000_Paper/RandomSuffleX/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesX <- random_files[grep('Performance',random_files)]
random_filesX <- random_filesX[which(!grepl('class',random_filesX))]
random_filesX <- random_filesX[which(!grepl('quantile',random_filesX))]
random_filesX <- random_filesX[which(!grepl('R2',random_filesX))]
all_random_files <- list(random_filesX,random_filesY)
all_corr_random_val <- data.frame()
for (i in 1:length(all_random_files)){
  files_corr <- all_random_files[[i]][grep('valPerformance',all_random_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==1) {
      file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==1){
      tmp <- tmp %>% mutate(lamda = 'shuffleX')
    }else{
      tmp <- tmp %>% mutate(lamda = 'shuffleY')
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_random_val <- rbind(all_corr_random_val,tmp)
  }
}
all_corr_random_train <- data.frame()
for (i in 1:length(all_random_files)){
  files_corr <- all_random_files[[i]][grep('trainPerformance',all_random_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==1) {
      file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==1){
      tmp <- tmp %>% mutate(lamda ='shuffleX')
    }else{
      tmp <- tmp %>% mutate(lamda = 'shuffleY')
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_random_train <- rbind(all_corr_random_train,tmp)
  }
}

df_corr_train <- all_corr_train %>% gather('TF','r',-lamda,-cell)
df_corr_train_random <-  all_corr_random_train %>% gather('TF','r',-lamda,-cell)
df_corr_train_random <- df_corr_train_random %>% group_by(cell,TF) %>% mutate(r=mean(r)) %>% ungroup() %>% unique()
df_corr_train_random <- df_corr_train_random %>% group_by(lamda,cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
#### to add ranking ###
df_corr_train <- df_corr_train %>% group_by(cell,TF) %>% 
  mutate(r=ifelse(lamda==0 | lamda==Inf,mean(r),
                  r)) %>% ungroup() %>% unique()
df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
df_corr_train <- rbind(df_corr_train,df_corr_train_random)
df_corr_train <- df_corr_train %>% filter(tf_rank<=15)
tfs_to_keep <- df_corr_train %>% select(lamda,cell,TF) %>% unique()
#colnames(tfs_to_keep)[3] <- 'TF_train'
### Filtering rank end ###
df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
#df_corr_train <- df_corr_train %>% select(-TF,-r) %>% unique()
df_corr_train$lamda <- factor(df_corr_train$lamda,
                              levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                         "0.001","0.005","0.01","0.05","Inf",
                                         "shuffleX","shuffleY"))
p3 <- ggboxplot(df_corr_train %>% select(lamda,mean_pear) %>% unique(), x = "lamda", y = "mean_pear",add='jitter')+
  ggtitle('Top 15% ranked TFs in each individual training set') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(limits = c(0,1))
p3_1 <- p3 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                   c('0','shuffleX'),
                                                   c('0','shuffleY')),
                              label = 'p.format',
                              label.y = c(0.88,0.92,0.96),
                              method = 'wilcox.test',
                              tip.length=0.05) 
#p3 <- p3 + stat_compare_means(method = "kruskal.test")
print(p3_1)

df_corr_val_random <-  all_corr_random_val %>% gather('TF','r',-lamda,-cell)
df_corr_val <- all_corr_val %>% gather('TF','r',-lamda,-cell)
df_corr_val_random <- df_corr_val_random %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val <- df_corr_val %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val <- rbind(df_corr_val,df_corr_val_random)
### filter based on training ###
df_corr_val <- left_join(tfs_to_keep,df_corr_val)
#df_corr_val <- df_corr_val %>% mutate(keep= (TF==TF_train)) %>% filter(keep==TRUE) %>% 
#  select(-TF_train,-keep) %>% unique()
####filtering ended ###
df_corr_val$lamda <- factor(df_corr_val$lamda,
                              levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                         "0.001","0.005","0.01","0.05","Inf",
                                         "shuffleX","shuffleY"))
#df_corr_val <- df_corr_val %>% select(-TF,-r) %>% unique()
p3 <- ggboxplot(df_corr_val %>% select(lamda,mean_pear) %>% unique(), x = "lamda", y = "mean_pear",add='jitter')+
  ggtitle('Filtered TFs to keep only 15% ranked TFs from training') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(breaks = seq(-0.4,0.8,0.2))
p3_2 <- p3 + stat_compare_means(comparisons = list(c('1e-04','0.005'),
                                                   c('0.001','0.005'),
                                                   c('0.005','0.01'),
                                                   c('0','Inf'),
                                                   c('0','shuffleX'),
                                                   c('0','shuffleY')),
                                label = 'p.format',
                                method = 'wilcox.test',
                                tip.length=0.05) 
#p3 <- p3 + stat_compare_means(method = "kruskal.test")
print(p3_2)


png('regularization_vs_filteredTFPerformance.png',units = 'in',width = 16,height = 12,res = 600)
p3_1 + p3_2
dev.off()

# Load sensitivity results------------------------
regs <- c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02)
path <- 'Model/CVL1000_Paper/hyperParamSensitivity/'
files <- list.files(path)
files <- data.frame(files)
files <- files %>% mutate(type=grepl('.csv',files)) %>% filter(type==T)
for (i in 1:nrow(files)){
  file <- files$files[i]
  if (i==1){
    sensitivity <- data.table::fread(paste0(path,file)) %>% select(-V1) %>% mutate(lambda=regs[i]) %>% unique()
  }else{
    sensitivity <- rbind(sensitivity,
                         data.table::fread(paste0(path,file)) %>% select(-V1) %>% mutate(lambda=regs[i]) %>% unique())
  }
}
# sensitivity <- sensitivity %>% filter(lambda>=0.00001)
sensitivity$lambda <- factor(sensitivity$lambda,
                             levels = regs)

minor_breaks <- rep(1:9, 5)*(10^rep(-4:1, each=9))
p3 <- ggplot(sensitivity,aes(x=`noise level`,y=`average pearson correlation`)) + 
  geom_point(aes(color=lambda))+ 
  geom_line(aes(color=lambda),linetype = 'dashed',lwd = 1) + 
  #geom_smooth(aes(color=lambda,fill=lambda),alpha=0.1,se=F) + 
  geom_errorbar(aes(ymin=`average pearson correlation`-2 * std_pear, 
                    ymax=`average pearson correlation`+2 * std_pear,
                    color=lambda),
                width=.05) +
  scale_x_log10(breaks=c(1e-4,1e-3,1e-2,1e-1,1,10),
                minor_breaks = minor_breaks,
                labels=c(1e-4,1e-3,1e-2,1e-1,1,10))+
  scale_y_continuous(breaks = seq(0,1,0.2),
                     labels = seq(0,1,0.2),limits = c(0,1.1))+
  ylab('average pearson correlation')+
  xlab('noise level') +
  ggtitle('Correlation between predicted and noise-perturbed data for increasing levels of noise')+
  guides(color=guide_legend(title="λ")) +
  annotate(geom="text", x=1e-03, y=1.03, label="Error bars demonstrate mean ± 2*SDs",size=4)+
  theme_minimal()+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(p3)

png('NoiseSensitivity4.png',units = 'in',width = 9,height = 9,res = 600)
print(p3)
dev.off()

#p3 <- p3 +geom_mark_circle(aes(filter = (`average pearson correlation`<=0.8 & `average pearson correlation`>=0.6)))

sensitivityFocused <- sensitivity# %>%  filter(`noise level`==0.5)
sensitivityFocused$`noise level` <- as.factor(sensitivityFocused$`noise level`)
minor_breaks <- rep(1:9, 8)*(10^rep(-8:-1, each=9))
p4 <- ggplot(sensitivityFocused,aes(x=as.numeric(as.character(lambda)),y=`average pearson correlation`)) + 
  geom_point(aes(color=`noise level`),size=2)+ 
  geom_line(aes(color=`noise level`),linetype = 'dashed',lwd = 1) + 
  #scale_colour_manual(values =c('black')) +
  scale_x_log10(breaks=c(1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1),
                minor_breaks = minor_breaks,
                labels=c(1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1))+
  scale_y_continuous(breaks = seq(0,1,0.2),
                     labels = seq(0,1,0.2),limits = c(0,1))+
  ylab('average pearson correlation')+
  xlab('λ') +
  ggtitle('Correlation between predicted and noise-perturbed data for increasing regularization')+
  theme_minimal()+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(p4)

png('NoiseSensitivity.png',units = 'in',width = 12,height = 9,res = 600)
p3/p4
dev.off()


# Correctly classifying inhibition/activation---------

## Calculate the probability of observing randomly an inhibited or activated TF in our data in each cell-line
TFoutput <- read.delim('Model/data/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv') %>% column_to_rownames('X')
cellInput <- read.delim('Model/data/L1000_lvl3-conditions_cells.tsv') %>% column_to_rownames('X')
TFoutput <- TFoutput[which(rownames(TFoutput) %in% rownames(cellInput)),]
gc()
inh <- sum(TFoutput<=0.25)/(nrow(TFoutput)*ncol(TFoutput))
act <- sum(TFoutput>=0.75)/(nrow(TFoutput)*ncol(TFoutput))
ggplot(TFoutput %>% gather('TF','activity'),aes(x=activity)) + geom_histogram(fill='#125b80',color='black',bins = 50)+
  geom_vline(xintercept = 0.25,color='red',lty='dashed',linewidth = 1)+
  geom_vline(xintercept = 0.75,color='red',lty='dashed',linewidth = 1)+
  annotate(geom = 'text',x=0.125,y=7600,label = paste0('Inhibited ',round(100*inh,2),'%'))+
  annotate(geom = 'text',x=0.875,y=7600,label = paste0('Activated ',round(100*act,2),'%'))+
  ggtitle('Data distribution') + ylab('counts')+
  theme(plot.title = element_text(hjust=0.5),text = element_text(size=12))
  
#TFoutput_binary <- -1*(1*TFoutput<=0.25) + 1*(1*TFoutput>=0.75)
max_iter <- 10000
for (j in 1:ncol(cellInput)){
  conds <- rownames(cellInput)[which(cellInput[,j]==1)]
  prec_inh_null <- NULL
  prec_mid_null <- NULL
  prec_act_null <- NULL
  f1_null <- NULL
  mcc_null <- NULL
  TFoutput_binary <- -1*(1*TFoutput[which(rownames(TFoutput) %in% conds),]<=0.25) + 1*(1*TFoutput[which(rownames(TFoutput) %in% conds),]>=0.75)
  for (i in 1:max_iter){
    df <- TFoutput[which(rownames(TFoutput) %in% conds),]
    df <- df [sample(1:nrow(df)),]
    df <- -1*(1*df<=0.25) + 1*(1*df>=0.75)
    confusion <- caret::confusionMatrix(factor(c(df),levels = c(-1,0,1)),
                                        factor(c(TFoutput_binary),levels = c(-1,0,1)),
                                        mode = "everything")
    mcc_null[i] <- mltools::mcc(factor(c(TFoutput_binary),levels = c(-1,0,1)),
                                factor(c(df),levels = c(-1,0,1)))
    
    prec_inh_null[i] <- confusion$byClass[1,5]
    prec_mid_null[i] <- confusion$byClass[2,5]
    prec_act_null[i] <- confusion$byClass[3,5]
    f1_null[i] <- MLmetrics::F1_Score_micro(factor(c(TFoutput_binary),levels = c(-1,0,1)),
                                            factor(c(df),levels = c(-1,0,1)))
    if (i %% 500 == 0 ){
      print(paste0('Finished iteration ',i,'/',max_iter))
    }
  }
  df_null <- data.frame(prec_inh_null,prec_mid_null,prec_act_null,f1_null,mcc_null)
  saveRDS(df_null,paste0('Model/CVL1000_Paper/random_classification_performance_',colnames(cellInput)[j],'.rds'))
  print(paste0('Finished cell ',colnames(cellInput)[j]))
}
#df_null <- data.frame(prec_inh_null,prec_mid_null,prec_act_null,f1_null,mcc_null)
#saveRDS(df_null,'Model/CVL1000_Paper/random_classification_performance.rds')
# Read saved simulations
df_null <- data.frame()
for (j in 1:ncol(cellInput)){
  df_null <- rbind(df_null,
                   readRDS(paste0('Model/CVL1000_Paper/random_classification_performance_',
                                  colnames(cellInput)[j],
                                  '.rds')) %>% mutate(cell=colnames(cellInput)[j]))
}

colnames(df_null) <- c('precision in inhibition','precision in middle bin',
                       'precision in activation',
                       'micro F1 score','MCC',
                       'cell')
df_null <- df_null %>% gather('metric','value',-cell)
p <- ggplot(df_null,aes(x=value,fill=cell)) + geom_histogram(color='black',alpha=0.8,bins = 50) + facet_wrap(~metric,scales = 'free')+
  ggtitle('Performance in randomly classifying correctly activation/inhibition') + theme(plot.title = element_text(hjust = 0.5))
print(p)
png('random_classification_performance.png',units = 'in',width = 12,height = 9,res=600)
print(p)
dev.off()
df_null <- df_null %>% group_by(metric,cell) %>% mutate(value=mean(value)) %>% ungroup() %>% unique()


### Read performance files
files <- list.files('Model/CVL1000_Paper/MultipleRuns/test/')
files <- files[grep('.csv',files)]
files <- files[grep('class',files)]
files_inf <- files[grep('modeltype1',files)]
files_zero <- files[grep('modeltype4',files)]
files_reg <-  list.files('Model/CVL1000_Paper/hyperParamTune/test/')
files_reg <- files_reg[grep('.csv',files_reg)]
files_reg <- files_reg[grep('class',files_reg)]
all_files  <- list(files_inf,files_zero,files_reg)

# Load train performance
all_class_train <- data.frame()
for (i in 1:length(all_files)){
  files_class <- all_files[[i]][grep('TrainPerformance',all_files[[i]])]
  #file <- files_class[1]
  for (file in files_class){
    cell <- str_split_fixed(file,'_',3)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_class_train <- rbind(all_class_train,tmp)
  }
}

# Load validation performance
all_class_val <- data.frame()
for (i in 1:length(all_files)){
  files_class <- all_files[[i]][grep('ValidationPerformance',all_files[[i]])]
  for (file in files_class){
    cell <- str_split_fixed(file,'_',3)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_class_val <- rbind(all_class_val,tmp)
  }
}

# Add random shuffle
random_files <- list.files('Model/CVL1000_Paper/RandomSuffleY/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('class',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesY <- random_files[grep('Performance',random_files)]
random_files <- list.files('Model/CVL1000_Paper/RandomSuffleX/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('class',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesX <- random_files[grep('Performance',random_files)]
all_random_files <- list(random_filesX,random_filesY)
all_class_random_val <- data.frame()
for (i in 1:length(all_random_files)){
  files_class <- all_random_files[[i]][grep('ValidationPerformance',all_random_files[[i]])]
  for (file in files_class){
    cell <- str_split_fixed(file,'_',3)[1,2]
    if (i==1) {
      file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==1){
      tmp <- tmp %>% mutate(lamda = 'shuffleX')
    }else{
      tmp <- tmp %>% mutate(lamda = 'shuffleY')
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_class_random_val <- rbind(all_class_random_val,tmp)
  }
}
all_class_random_train <- data.frame()
for (i in 1:length(all_random_files)){
  files_class <- all_random_files[[i]][grep('TrainPerformance',all_random_files[[i]])]
  for (file in files_class){
    cell <- str_split_fixed(file,'_',3)[1,2]
    if (i==1) {
      file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==1){
      tmp <- tmp %>% mutate(lamda ='shuffleX')
    }else{
      tmp <- tmp %>% mutate(lamda = 'shuffleY')
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_class_random_train <- rbind(all_class_random_train,tmp)
  }
}

df_class_train <- rbind(all_class_train,all_class_random_train)
df_class_train$lamda <- factor(df_class_train$lamda,
                              levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                         "0.001","0.005","0.01","0.05","Inf",
                                         "shuffleX","shuffleY"))
colnames(df_class_train) <- c('lamda','precision in inhibition','precision in middle bin',
                              'precision in activation','micro F1 score','MCC','cell')
df_class_train <- df_class_train %>% gather('metric','value',-lamda,-cell)
p3_1 <- ggboxplot(df_class_train, x = "lamda", y = "value",add='jitter')+
  ggtitle('Classification performance in training') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('value') + scale_y_continuous(limits = c(0,1))+
  facet_wrap(~metric) + geom_hline(data = df_null,aes(yintercept=value,group =metric,color=cell),linewidth=1,linetype='dashed')

p3_1 <- p3_1 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                   c('0','shuffleX'),
                                                   c('0','shuffleY')),
                                label = 'p.format',
                                label.y = c(0.82,0.88,0.94),
                                method = 'wilcox.test',
                                tip.length=0.05)
print(p3_1)

df_class_val <- rbind(all_class_val,all_class_random_val)
df_class_val$lamda <- factor(df_class_val$lamda,
                               levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                          "0.001","0.005","0.01","0.05","Inf",
                                          "shuffleX","shuffleY"))
colnames(df_class_val) <- c('lamda','precision in inhibition','precision in middle bin',
                              'precision in activation','micro F1 score','MCC','cell')
df_class_val <- df_class_val %>% gather('metric','value',-lamda,-cell)
p3 <- ggboxplot(df_class_val, x = "lamda", y = "value",add='jitter')+
  ggtitle('Classification performance in validation') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('value') + scale_y_continuous(limits = c(0,1))+
  facet_wrap(~metric)+ geom_hline(data = df_null,aes(yintercept=value,group =metric,color=cell),linewidth=1,linetype='dashed')

p3_2 <- p3 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                   c('0','shuffleX'),
                                                   c('0','shuffleY')),
                                label = 'p.format',
                                label.y = c(0.78,0.84,0.9),
                                method = 'wilcox.test',
                                tip.length=0.05) 
print(p3_2)

png('regularization_vs_classification.png',units = 'in',width = 18,height = 18,res = 600)
p3_1 / p3_2
dev.off()


# Performance in the first PCs--------------------------
files <- list.files('Model/CVL1000_Paper/MultipleRuns/test/')
files <- files[grep('.csv',files)]
files_inf <- files[grep('modeltype1',files)]
files_inf <- files_inf[which(!grepl('class',files_inf))]
files_inf <- files_inf[which(!grepl('quantile',files_inf))]
files_zero <- files[grep('modeltype4',files)]
files_zero <- files_zero[which(!grepl('class',files_zero))]
files_zero <- files_zero[which(!grepl('quantile',files_zero))]
files_zero <- files_zero[which(!grepl('R2',files_zero))]
files_reg <-  list.files('Model/CVL1000_Paper/hyperParamTune/test/')
files_reg <- files_reg[grep('.csv',files_reg)]
files_reg <- files_reg[which(!grepl('class',files_reg))]
files_reg <- files_reg[which(!grepl('quantile',files_reg))]
all_files  <- list(files_inf,files_zero,files_reg)
# Load train correlation data
all_corr_train <- data.frame()
for (i in 1:length(all_files)){
  files_corr <- all_files[[i]][grep('trainPerformance',all_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_train <- rbind(all_corr_train,tmp)
  }
}
# Load validation correlation data
all_corr_val <- data.frame()
for (i in 1:length(all_files)){
  files_corr <- all_files[[i]][grep('valPerformance',all_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==3) {
      file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==3){
      tmp <- tmp %>% mutate(lamda = c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02))
    }else if (i==1){
      tmp <- tmp %>% mutate(lamda=Inf)
    }else{
      tmp <- tmp %>% mutate(lamda=0)
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_val <- rbind(all_corr_val,tmp)
  }
}
# Add random shuffle
random_files <- list.files('Model/CVL1000_Paper/RandomSuffleY/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesY <- random_files[grep('Performance',random_files)]
random_filesY <- random_filesY[which(!grepl('class',random_filesY))]
random_filesY <- random_filesY[which(!grepl('quantile',random_filesY))]
random_filesY <- random_filesY[which(!grepl('R2',random_filesY))]
random_files <- list.files('Model/CVL1000_Paper/RandomSuffleX/test/')
random_files <- random_files[grep('.csv',random_files)]
random_files <- random_files[grep('modeltype4',random_files)]
random_filesX <- random_files[grep('Performance',random_files)]
random_filesX <- random_filesX[which(!grepl('class',random_filesX))]
random_filesX <- random_filesX[which(!grepl('quantile',random_filesX))]
random_filesX <- random_filesX[which(!grepl('R2',random_filesX))]
all_random_files <- list(random_filesX,random_filesY)
all_corr_random_val <- data.frame()
for (i in 1:length(all_random_files)){
  files_corr <- all_random_files[[i]][grep('valPerformance',all_random_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==1) {
      file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==1){
      tmp <- tmp %>% mutate(lamda = 'shuffleX')
    }else{
      tmp <- tmp %>% mutate(lamda = 'shuffleY')
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_random_val <- rbind(all_corr_random_val,tmp)
  }
}
all_corr_random_train <- data.frame()
for (i in 1:length(all_random_files)){
  files_corr <- all_random_files[[i]][grep('trainPerformance',all_random_files[[i]])]
  for (file in files_corr){
    cell <- str_split_fixed(file,'_',4)[1,2]
    if (i==1) {
      file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
    }else{
      file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
    }
    tmp <- data.table::fread(file)
    colnames(tmp)[1] <- 'lamda'
    if (i==1){
      tmp <- tmp %>% mutate(lamda ='shuffleX')
    }else{
      tmp <- tmp %>% mutate(lamda = 'shuffleY')
    }
    tmp <- tmp %>% mutate(cell=cell)
    all_corr_random_train <- rbind(all_corr_random_train,tmp)
  }
}
# Find PCs per cell-line and global PCA
TFoutput <- read.delim('Model/data/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv') %>% column_to_rownames('X')
cellInput <- read.delim('Model/data/L1000_lvl3-conditions_cells.tsv') %>% column_to_rownames('X')
TFoutput <- TFoutput[which(rownames(TFoutput) %in% rownames(cellInput)),]
gc()
pca_all <- prcomp(TFoutput)
png('tfs_pca_all_screeplot.png',units = 'in',width = 9,height = 6,res=600)
fviz_screeplot(pca_all)
dev.off()
fviz_contrib(pca_all, choice = "var", axes = 1:3, top = 20)
loadings <- apply(abs(pca_all$rotation[,1:5]),1,mean)
loadings <- loadings[order(-loadings)]
print(names(loadings[1:20]))

### Filter top 20 TFs in PCs and visualize performance
df_corr_train <- all_corr_train %>% gather('TF','r',-lamda,-cell)
df_corr_train <- df_corr_train %>% filter(TF %in% names(loadings[1:20])) %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_train_random <-  all_corr_random_train %>% gather('TF','r',-lamda,-cell)
df_corr_train_random <-  df_corr_train_random %>% filter(TF %in% names(loadings[1:20])) %>% 
  group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_train <- df_corr_train %>% select(-TF,-r) %>% unique()
df_corr_train_random <- df_corr_train_random %>% select(-TF,-r) %>% unique()
df_corr_train <- rbind(df_corr_train,df_corr_train_random)
df_corr_train$lamda <- factor(df_corr_train$lamda,
                              levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                         "0.001","0.005","0.01","0.05","Inf",
                                         "shuffleX","shuffleY"))
p3_1 <- ggboxplot(df_corr_train, x = "lamda", y = "mean_pear",add='jitter')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Performance in train of top 20 TFs contributing in the three first PCs') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(limits = c(-0.1,1))
p3_1 <- p3_1 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                   c('0','shuffleX'),
                                                   c('0','shuffleY')),
                                label = 'p.format',
                                label.y = c(0.75,0.8,0.85),
                                method = 'wilcox.test',
                                tip.length=0.05) 
print(p3_1)
df_corr_val <- all_corr_val %>% gather('TF','r',-lamda,-cell)
df_corr_val <- df_corr_val %>% filter(TF %in% names(loadings[1:20])) %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val_random <-  all_corr_random_val %>% gather('TF','r',-lamda,-cell)
df_corr_val_random <-  df_corr_val_random %>% filter(TF %in% names(loadings[1:20])) %>% 
  group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val <- df_corr_val %>% select(-TF,-r) %>% unique()
df_corr_val_random <- df_corr_val_random %>% select(-TF,-r) %>% unique()
df_corr_val <- rbind(df_corr_val,df_corr_val_random)
df_corr_val$lamda <- factor(df_corr_val$lamda,
                              levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                         "0.001","0.005","0.01","0.05","Inf",
                                         "shuffleX","shuffleY"))
p3_2 <- ggboxplot(df_corr_val, x = "lamda", y = "mean_pear",add='jitter')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Performance in validation of top 20 TFs contributing in the three first PCs') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(limits = c(-0.1,1))
p3_2 <- p3_2 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                     c('0','shuffleX'),
                                                     c('0','shuffleY')),
                                  label = 'p.format',
                                  label.y = c(0.5,0.55,0.6),
                                  method = 'wilcox.test',
                                  tip.length=0.05) 
print(p3_2)
png('global_pca_filtered_perTF_performance.png',units = 'in',width = 12,height = 12,res=600)
print(p3_1/p3_2)
dev.off()

### Filter PCs per cell-line
tfs_to_keep <- readRDS('top_15perc_rank_tfs_persplit.rds')
df_corr_train <- data.frame()
df_corr_val <- data.frame()
all_corr_train <- rbind(all_corr_train,all_corr_random_train)
all_corr_val <- rbind(all_corr_val,all_corr_random_val)
tfs_loadings <- data.frame()
tfs_not_loadings <- data.frame()
p.vals <- rep(1,ncol(cellInput))
for (j in 1:ncol(cellInput)){
  conds <- rownames(cellInput)[which(cellInput[,j]==1)]
  y <- TFoutput[which(rownames(TFoutput) %in% conds),]
  res_pca <- prcomp(y)
  loadings <- apply(abs(res_pca$rotation[,1:5]),1,mean)
  loadings <- loadings[order(-loadings)]
  tfs <- names(loadings[1:20])
  tfs_loadings <- rbind(tfs_loadings,
                        data.frame(cell=colnames(cellInput)[j],TF=tfs,pca = 'in_loadings'))
  tfs_not_loadings <- rbind(tfs_not_loadings,
                        data.frame(cell=colnames(cellInput)[j],
                                   TF=colnames(TFoutput)[which(!(colnames(TFoutput) %in% tfs))],
                                   pca = 'not_in_loadings'))
  tfs_contigency <- left_join(rbind(tfs_loadings,tfs_not_loadings) %>% filter(cell==colnames(cellInput)[j]) %>%
                                select(-cell) %>% unique(),
                              tfs_to_keep %>% filter(lamda==0.005) %>%
                                filter(cell==colnames(cellInput)[j]) %>% select(-lamda,-cell) %>% 
                                unique() %>% mutate(rank='in_top'),
                              by='TF')
  tfs_contigency <- tfs_contigency %>% mutate(rank = ifelse(is.na(rank),'not_in_top',rank))
  
  total_cont <- table(tfs_contigency %>% select(-TF))
  if (ncol(total_cont)<2){
    total_cont <- cbind(total_cont,c(0,0))
    colnames(total_cont)[2] <- 'in_top'
  }
  test_res <- fisher.test(total_cont)
  p.vals[j] <- test_res$p.value
  
  tmp_train <- all_corr_train %>% gather('TF','r',-lamda,-cell)
  tmp_train <- tmp_train %>% filter(cell==colnames(cellInput)[j]) %>% filter(TF %in% tfs)
  df_corr_train <- rbind(df_corr_train,tmp_train)
  
  tmp_val <- all_corr_val %>% gather('TF','r',-lamda,-cell)
  tmp_val <- tmp_val %>% filter(cell==colnames(cellInput)[j])%>% filter(TF %in% tfs)
  df_corr_val <- rbind(df_corr_val,tmp_val)
}
hist(p.vals)
hist(p.adjust(p.vals,'BH'))
tfs_contigency <- left_join(rbind(tfs_loadings,tfs_not_loadings),tfs_to_keep %>% filter(lamda==0.005) %>%
                              select(-lamda) %>% mutate(rank='in_top'),
                            by=c('cell','TF'))
tfs_contigency <- tfs_contigency %>% mutate(rank = ifelse(is.na(rank),'not_in_top',rank))
total_cont <- table(tfs_contigency %>% select(-cell,-TF))
test_res <- fisher.test(total_cont)
png('contigency_pcaloadings_toprankedTFs_total.png',units = 'in',width = 9,height = 9,res=600)
ggplot(as.data.frame(total_cont), aes(pca,rank, fill= Freq)) + 
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient2() +
  labs(x = "top 20 in the loadings of 5PCs",y = "in top 15% ranked in performace") + 
  ggtitle('Confusion matrix for overlap of PCA-important and top ranked TFs')+
  annotate('text',x=1.1,y=2.55,label=paste0('Fisher exact test: p.value=',test_res$p.value))+
  theme_minimal(base_family = "serif",base_size = 15)
dev.off()

df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_train <- df_corr_train %>% select(-TF,-r) %>% unique()
df_corr_val <- df_corr_val %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r)) %>% ungroup()
df_corr_val <- df_corr_val %>% select(-TF,-r) %>% unique()

df_corr_train$lamda <- factor(df_corr_train$lamda,
                              levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                         "0.001","0.005","0.01","0.05","Inf",
                                         "shuffleX","shuffleY"))
p3_1 <- ggboxplot(df_corr_train, x = "lamda", y = "mean_pear",add='jitter')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Performance in train of top 20 TFs contributing in the three first PCs') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(limits = c(-0.1,1))
p3_1 <- p3_1 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                     c('0','shuffleX'),
                                                     c('0','shuffleY')),
                                  label = 'p.format',
                                  label.y = c(0.75,0.8,0.85),
                                  method = 'wilcox.test',
                                  tip.length=0.05) 
print(p3_1)
df_corr_val$lamda <- factor(df_corr_val$lamda,
                            levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
                                       "0.001","0.005","0.01","0.05","Inf",
                                       "shuffleX","shuffleY"))
p3_2 <- ggboxplot(df_corr_val, x = "lamda", y = "mean_pear",add='jitter')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5)) +
  ggtitle('Performance in validation of top 20 TFs contributing in the three first PCs') +
  theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
  xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(limits = c(-0.1,1))
p3_2 <- p3_2 + stat_compare_means(comparisons = list(c('0','Inf'),
                                                     c('0','shuffleX'),
                                                     c('0','shuffleY')),
                                  label = 'p.format',
                                  label.y = c(0.5,0.55,0.6),
                                  method = 'wilcox.test',
                                  tip.length=0.05) 
print(p3_2)
png('cellLine_pca_filtered_perTF_performance.png',units = 'in',width = 12,height = 12,res=600)
print(p3_1/p3_2)
dev.off()

# # Binned performance-------------
# # Add random shuffle
# random_files <- list.files('Model/CVL1000_Paper/RandomSuffleY/test/')
# random_files <- random_files[grep('.csv',random_files)]
# random_files <- random_files[grep('modeltype4',random_files)]
# random_files <- random_files[grep('quantile',random_files)]
# random_filesY <- random_files[grep('Performance',random_files)]
# random_files <- list.files('Model/CVL1000_Paper/RandomSuffleX/test/')
# random_files <- random_files[grep('.csv',random_files)]
# random_files <- random_files[grep('modeltype4',random_files)]
# random_files <- random_files[grep('quantile',random_files)]
# random_filesX <- random_files[grep('Performance',random_files)]
# all_random_files <- list(random_filesX,random_filesY)
# all_corr_random_val <- data.frame()
# for (i in 1:length(all_random_files)){
#   files_corr <- all_random_files[[i]][grep('ValidationPerformance',all_random_files[[i]])]
#   for (file in files_corr){
#     cell <- str_split_fixed(file,'_',4)[1,2]
#     if (i==1) {
#       file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
#     }else{
#       file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
#     }
#     tmp <- data.table::fread(file)
#     colnames(tmp)[1] <- 'TF'
#     if (i==1){
#       tmp <- tmp %>% mutate(lamda = 'shuffleX')
#     }else{
#       tmp <- tmp %>% mutate(lamda = 'shuffleY')
#     }
#     tmp <- tmp %>% mutate(cell=cell)
#     all_corr_random_val <- rbind(all_corr_random_val,tmp)
#   }
# }
# all_corr_random_train <- data.frame()
# for (i in 1:length(all_random_files)){
#   files_corr <- all_random_files[[i]][grep('TrainPerformance',all_random_files[[i]])]
#   for (file in files_corr){
#     cell <- str_split_fixed(file,'_',4)[1,2]
#     if (i==1) {
#       file <- paste0('Model/CVL1000_Paper/RandomSuffleX/test/',file)
#     }else{
#       file <- paste0('Model/CVL1000_Paper/RandomSuffleY/test/',file)
#     }
#     tmp <- data.table::fread(file)
#     colnames(tmp)[1] <- 'TF'
#     if (i==1){
#       tmp <- tmp %>% mutate(lamda ='shuffleX')
#     }else{
#       tmp <- tmp %>% mutate(lamda = 'shuffleY')
#     }
#     tmp <- tmp %>% mutate(cell=cell)
#     all_corr_random_train <- rbind(all_corr_random_train,tmp)
#   }
# }
# 
# ### Load other data
# files <- list.files('Model/CVL1000_Paper/MultipleRuns/test/')
# files <- files[grep('.csv',files)]
# files_inf <- files[grep('modeltype1',files)]
# files_inf <- files_inf[grep('quantile',files_inf)]
# files_zero <- files[grep('modeltype4',files)]
# files_zero <- files_zero[grep('quantile',files_zero)]
# 
# files_reg <-  list.files('Model/CVL1000_Paper/hyperParamTune/test/')
# files_reg <- files_reg[grep('.csv',files_reg)]
# files_reg <- files_reg[grep('quantile',files_reg)]
# all_files  <- list(files_inf,files_zero,files_reg)
# 
# # Load train correlation data
# all_corr_train <- data.frame()
# lamdas <- c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02)
# l <- 1
# cell_previous <- 'nothing'
# for (i in 1:length(all_files)){
#   files_corr <- all_files[[i]][grep('TrainPerformance',all_files[[i]])]
#   for (file in files_corr){
#     cell <- str_split_fixed(file,'_',4)[1,2]
#     if (cell_previous!=cell){
#       l <- 1
#       cell_previous <- cell
#     }
#     if (i==3) {
#       file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
#     }else{
#       file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
#     }
#     tmp <- data.table::fread(file)
#     colnames(tmp)[1] <- 'TF'
#     if (i==3){
#       tmp <- tmp %>% mutate(lamda = lamdas[l])
#       l <- l+1
#     }else if (i==1){
#       tmp <- tmp %>% mutate(lamda=Inf)
#     }else{
#       tmp <- tmp %>% mutate(lamda=0)
#     }
#     tmp <- tmp %>% mutate(cell=cell)
#     all_corr_train <- rbind(all_corr_train,tmp)
#   }
# }
# # Load validation correlation data
# all_corr_val <- data.frame()
# lamdas <- c(1e-08, 1e-07, 1e-06, 1e-05, 1e-04, 1e-03, 5e-03, 1e-02, 5e-02)
# l <- 1
# cell_previous <- 'nothing'
# for (i in 1:length(all_files)){
#   files_corr <- all_files[[i]][grep('ValidationPerformance',all_files[[i]])]
#   for (file in files_corr){
#     cell <- str_split_fixed(file,'_',4)[1,2]
#     if (cell_previous!=cell){
#       l <- 1
#       cell_previous <- cell
#     }
#     if (i==3) {
#       file <- paste0('Model/CVL1000_Paper/hyperParamTune/test/',file)
#     }else{
#       file <- paste0('Model/CVL1000_Paper/MultipleRuns/test/',file)
#     }
#     tmp <- data.table::fread(file)
#     colnames(tmp)[1] <- 'TF'
#     if (i==3){
#       tmp <- tmp %>% mutate(lamda = lamdas[l])
#     }else if (i==1){
#       tmp <- tmp %>% mutate(lamda=Inf)
#     }else{
#       tmp <- tmp %>% mutate(lamda=0)
#     }
#     tmp <- tmp %>% mutate(cell=cell)
#     all_corr_val <- rbind(all_corr_val,tmp)
#   }
# }
# 
# 
# df_corr_train <- all_corr_train %>% gather('quantile','r',-lamda,-cell,-TF)
# df_corr_train <- df_corr_train %>% filter(quantile=="q0.25")
# df_corr_train_random <-  all_corr_random_train %>% gather('quantile','r',-lamda,-cell,-TF)
# df_corr_train_random <-df_corr_train_random %>% filter(quantile=="q0.25")
# df_corr_train_random <- df_corr_train_random %>% group_by(cell,TF) %>% mutate(r=mean(r,na.rm = T)) %>% ungroup() %>% unique()
# # df_corr_train_random <- df_corr_train_random %>% filter(!is.na(df_corr_train_random))
# df_corr_train_random <- df_corr_train_random %>% group_by(lamda,cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
# #### to add ranking ###
# df_corr_train <- df_corr_train %>% group_by(cell,TF) %>% 
#   mutate(r=ifelse(lamda==0 | lamda==Inf,mean(r),
#                   r)) %>% ungroup() %>% unique()
# df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>%  mutate(tf_rank=rank(-r)) %>% ungroup()
# df_corr_train <- rbind(df_corr_train,df_corr_train_random)
# df_corr_train <- df_corr_train %>% filter(tf_rank<=15)
# tfs_to_keep <- df_corr_train %>% select(lamda,cell,TF) %>% unique()
# #colnames(tfs_to_keep)[3] <- 'TF_train'
# ### Filtering rank end ###
# df_corr_train <- df_corr_train %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r,na.rm = T)) %>% ungroup()
# #df_corr_train <- df_corr_train %>% select(-TF,-r) %>% unique()
# df_corr_train$lamda <- factor(df_corr_train$lamda,
#                               levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
#                                          "0.001","0.005","0.01","0.05","Inf",
#                                          "shuffleX","shuffleY"))
# p3 <- ggboxplot(df_corr_train %>% select(lamda,mean_pear) %>% unique(), x = "lamda", y = "mean_pear",add='jitter')+
#   ggtitle('Top 15% ranked TFs in each individual training set') +
#   theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
#   xlab('λ') + ylab('average across TFs pearson`s r') + scale_y_continuous(limits = c(0,1))
# p3_1 <- p3 + stat_compare_means(comparisons = list(c('0','Inf'),
#                                                    c('0','shuffleX'),
#                                                    c('0','shuffleY')),
#                                 label = 'p.format',
#                                 label.y = c(0.88,0.92,0.96),
#                                 method = 'wilcox.test',
#                                 tip.length=0.05) 
# #p3 <- p3 + stat_compare_means(method = "kruskal.test")
# print(p3_1)
# 
# df_corr_val_random <-  all_corr_random_val %>% gather('quantile','r',-lamda,-cell,-TF)
# df_corr_val_random  <- df_corr_val_random %>% filter(quantile=="q0.25")
# df_corr_val <- all_corr_val %>% gather('quantile','r',-lamda,-cell,-TF)
# df_corr_val <- df_corr_val %>% filter(quantile=="q0.25")
# df_corr_val_random <- df_corr_val_random %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r,na.rm = T)) %>% ungroup()
# df_corr_val <- df_corr_val %>% group_by(lamda,cell) %>% mutate(mean_pear = mean(r,na.rm = T)) %>% ungroup()
# df_corr_val <- rbind(df_corr_val,df_corr_val_random)
# ### filter based on training ###
# df_corr_val <- left_join(tfs_to_keep,df_corr_val)
# #df_corr_val <- df_corr_val %>% mutate(keep= (TF==TF_train)) %>% filter(keep==TRUE) %>% 
# #  select(-TF_train,-keep) %>% unique()
# ####filtering ended ###
# df_corr_val$lamda <- factor(df_corr_val$lamda,
#                             levels = c("0","1e-08","1e-07","1e-06","1e-05","1e-04",
#                                        "0.001","0.005","0.01","0.05","Inf",
#                                        "shuffleX","shuffleY"))
# #df_corr_val <- df_corr_val %>% select(-TF,-r) %>% unique()
# p3 <- ggboxplot(df_corr_val %>% select(lamda,mean_pear) %>% unique(), x = "lamda", y = "mean_pear",add='jitter')+
#   ggtitle('Filtered TFs to keep only 15% ranked TFs from training') +
#   theme(text = element_text(size=11),plot.title = element_text(hjust = 0.5)) +
#   xlab('λ') + ylab('average across TFs pearson`s r')
# p3_2 <- p3 + stat_compare_means(comparisons = list(c('1e-04','0.005'),
#                                                    c('0.001','0.005'),
#                                                    c('0.005','0.01'),
#                                                    c('0','Inf'),
#                                                    c('0','shuffleX'),
#                                                    c('0','shuffleY')),
#                                 label = 'p.format',
#                                 method = 'wilcox.test',
#                                 tip.length=0.05) 
# #p3 <- p3 + stat_compare_means(method = "kruskal.test")
# print(p3_2)