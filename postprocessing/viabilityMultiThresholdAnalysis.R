library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(patchwork)
# library(PharmacoGx)
# library(caret)
source("viabilityModeling.R")


# Initialize
set <- 'NCI60_2021'
no_models <- 50


# Load data---------

# Load cmap data
drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]
data_cmap <- readRDS("../preprocessing/preprocessed_data/all_cmap_sigs_with_pert_info.rds")
TF_activities_2 <- read.delim('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv',row.names = 1)

data_cmap <- data_cmap %>% filter(sig_id %in% rownames(TF_activities_2))
data_cmap <- data_cmap %>% filter(canonical_smiles %in% rownames(drug_targets) | pert_iname %in% rownames(drug_targets))
gc()

# Load lethality/viability----------------------
data <- readRDS(paste0('./',set,'.rds'))

# Create initial linear models using different number of models----------------
# Do it for multiple thresholds to check. And then choose one
initial_models <- c('elasticnet')
performance_df <- data.frame()
for (md in initial_models){
  test_corr <- NULL
  test_pval <- NULL
  for (th in 1:no_models){
    drug_targets <- read.delim('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv',row.names = 1)
    merged_interactions_all <- data.table::fread('../results/A549_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
    consensus_info <- readRDS('../results/consensus_interaction_info.rds')
    colnames(consensus_info)[2] <- 'variable'
    merged_interactions_all <- left_join(consensus_info,merged_interactions_all)
    merged_interactions_all <- merged_interactions_all %>% filter(!is.na(`Prior knowledge`))
    merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = 1*(Inferred=='Interaction'))
    merged_interactions_all <- merged_interactions_all %>% mutate(`Prior knowledge` = 1*(`Prior knowledge`=='Interaction'))
    merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred",
                                                                  "model_no","consensus_inferrence","mean_frequency")
    merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = ifelse((consensus_inferrence==3) & (mean_frequency>=th/no_models),1,0))
    merged_interactions_all <- merged_interactions_all %>% select(-mean_frequency,-consensus_inferrence)
    merged_interactions_all <- distinct(merged_interactions_all)
    
    inferred_interactions <- merged_interactions_all %>% select(drug,variable,Inferred)
    inferred_interactions <- distinct(inferred_interactions)
    inferred_interactions <- inferred_interactions %>% group_by(drug) %>% mutate(interactions = sum(Inferred)) %>% ungroup()
    inferred_interactions <- inferred_interactions %>% filter(interactions>0)
    inferred_interactions <- inferred_interactions %>% select(-interactions)
    inferred_interactions <- inferred_interactions %>% spread('variable','Inferred')
    gc()
    drug_targets_inferred <- as.matrix(distinct(inferred_interactions %>% column_to_rownames('drug')))
    drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions_all$drug),
                                 which(colnames(drug_targets) %in% merged_interactions_all$variable)]
    
    data_new <- data
    data_new <- data_new %>% filter(smiles %in% merged_interactions_all$drug)
    data_new[,3:(ncol(data_new)-9)] <- apply(data_new[,3:(ncol(data_new)-9)],c(1,2),as.numeric)
    changed <- 0
    for (i in 1:nrow(data_new)){
      drug <- data_new$smiles[i]
      if (drug %in% rownames(drug_targets_inferred)){
        for (j in 3:(ncol(data_new)-9)){
          target <-   colnames(data_new)[j]
          val <- drug_targets_inferred[drug,target]
          if (val!=data_new[i,j]){
            changed <- changed+1
          }
          data_new[i,j] <- val
        }
      }
    }
    ### Call function ###
    total_results <- viability_model(data_new,
                                     colnames(drug_targets),
                                     dataset=set,
                                     lethality_data_path = "./",
                                     model = md,
                                     no_models=50)
    
    res_test_all <- total_results$res_test_all
    res_test_all <- res_test_all %>% group_by(data) %>% mutate(correlation = cor(predicted,EC50)) %>% ungroup()
    res_test_all <- res_test_all %>% mutate(frequency=th/no_models)
    res_test_all <- res_test_all %>% mutate(model=md)
    performance_df <- rbind(performance_df,
                            res_test_all)
    # saveRDS(performance_df,'../results/viability_analysis_results/viability_results_LOOCV_multi_threshold_analysis.rds')
    print(paste0('Using ',md,' finished: ',th,'/',no_models))
  }
}
# saveRDS(performance_df,'../results/viability_analysis_results/viability_results_LOOCV_multi_threshold_analysis.rds')

## plot results----------------------------------------------
performance_df <- readRDS('../results/viability_analysis_results/viability_results_LOOCV_multi_threshold_analysis.rds')
performance_df$data <- factor(performance_df$data,levels = c('test','shuffled'))
p1 <- ggplot(performance_df %>% select(-predicted,-EC50) %>% unique()%>% 
         filter(model %in% c('elasticnet','lasso','xgbTree','svmLinear','gaussprLinear')), #,'rf','neuralnet'
       aes(x=frequency,y=correlation,color=model,shape=data)) +
  geom_point(size=2) + geom_line(aes(linetype=data),size=1) +
  ylab('pearson`s r')+
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')+
  theme_pubr(base_size = 24,base_family = 'Arial') +
  theme(text = element_text(size = 24,family = 'Arial'),
        legend.text = element_text(size = 20,family = 'Arial'))
print(p1)
ggsave(paste0('../results/viability_analysis_results/performance_across_thresholds_best_lines.eps'),
       plot=p1,
       device = cairo_ps,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)
p2 <- ggscatter(performance_df %>% filter(data=='test')%>% select(-predicted,-EC50) %>% unique(),
          x='frequency',y='correlation',
          cor.coef = T,cor.coef.size = 7,cor.coef.coord = c(0,1),rug = T) +
  ylab('pearson`s r')+
  # geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  geom_smooth(se=T,color='red',size=2,linetype = 'dashed')+
  theme(text = element_text(size = 24,family = 'Arial')) +
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')
print(p2)
ggsave(paste0('../results/viability_analysis_results/performance_across_thresholds_scatter.eps'),
       plot=p2,
       device = cairo_ps,
       height = 9,
       width = 12,
       units = 'in',
       dpi = 600)
### Compare across all thresholds
p3 <- ggboxplot(performance_df %>% select(-predicted,-EC50) %>% unique(),x='model',y='correlation',color='model',add='jitter') +
  ylab('pearson`s r')+
  theme(text = element_text(size = 20,family = 'Arial'),axis.text.x = element_blank(),axis.title.x = element_blank()) +
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')+
  facet_wrap(~data)
print(wilcox.test(performance_df$correlation[which(performance_df$data=='test')],
                 performance_df$correlation[which(performance_df$data!='test')]))
# p <- (p1+p2)/p3
print(p3)
ggsave(paste0('../results/viability_analysis_results/performance_across_thresholds_boxplot.eps'),
       plot=p3,
       device = cairo_ps,
       height = 9,
       width = 12,
       units = 'in',
       dpi = 600)
### Check frequency and prior models
frequency_files <- list.files('../results/viability_analysis_results/',pattern = 'predictions_all_frequency')
frequency_performance <- data.frame()
for (file in frequency_files){
  md <- str_split_fixed(file,pattern = '.rds',n=2)[1,1]
  md <- str_split_fixed(md,pattern = '_',n=4)[1,4]
  frequency_performance <- rbind(frequency_performance,
                                 readRDS(paste0('../results/viability_analysis_results/',file)) %>% mutate(model=md))
}
frequency_performance <- frequency_performance %>% filter(!(model %in% c('LASSO_CTRPV2_2015','LASSO_PRISM_2020')))
frequency_performance <- frequency_performance %>% group_by(data,model) %>% mutate(correlation = cor(predicted,EC50)) %>% ungroup()
tmp <- frequency_performance %>% select(-predicted,-EC50) %>%
  filter(model!='LM') %>% filter(data=='test') %>% unique() 
frequency_performance$model <- factor(frequency_performance$model,levels = tmp$model[order(-tmp$correlation)])
frequency_performance$data <- factor(frequency_performance$data,levels = c('test','shuffled'))
ggplot(frequency_performance %>% select(-predicted,-EC50) %>%
         filter(model!='LM') %>% unique(),aes(x=model,y=correlation,fill=model)) +
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')+
  scale_y_continuous(breaks = c(-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25),limits = c(-0.75,1))+
  theme_pubr(base_size = 24,base_family = 'Arial') +
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 18,family = 'Arial')) +
  facet_wrap(~data)
ggsave(paste0('../results/viability_analysis_results/performance_frequency_all_models_barplot.eps'),
       device = cairo_ps,
       height = 10.5,
       width = 14,
       units = 'in',
       dpi = 600)
# prior models
prior_files <- list.files('../results/viability_analysis_results/',pattern = 'predictions_all_prior')
prior_performance <- data.frame()
for (file in prior_files){
  md <- str_split_fixed(file,pattern = '.rds',n=2)[1,1]
  md <- str_split_fixed(md,pattern = '_',n=4)[1,4]
  prior_performance <- rbind(prior_performance,
                                 readRDS(paste0('../results/viability_analysis_results/',file)) %>% mutate(model=md))
}
prior_performance <- prior_performance %>% filter(!(model %in% c('LASSO_CTRPV2_2015','LASSO_PRISM_2020')))
prior_performance <- prior_performance %>% mutate(model=ifelse(grepl('NCI60',model),'RF',model))
prior_performance <- prior_performance %>% group_by(data,model) %>% mutate(correlation = cor(predicted,EC50)) %>% ungroup()
tmp <- prior_performance %>% select(-predicted,-EC50) %>%
  filter(model!='LM') %>% filter(data=='test') %>% unique() 
prior_performance$model <- factor(prior_performance$model,levels = tmp$model[order(-tmp$correlation)])
prior_performance$data <- factor(prior_performance$data,levels = c('test','shuffled'))
ggplot(prior_performance %>% select(-predicted,-EC50) %>%
         filter(model!='LM') %>% unique(),aes(x=model,y=correlation,fill=model)) +
  geom_bar(stat = 'identity')+
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black')+
  scale_y_continuous(breaks = c(-1.25,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25),limits = c(-0.9,1))+
  theme_pubr(base_size = 24,base_family = 'Arial') +
  theme(text = element_text(size = 24,family = 'Arial'),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 18,family = 'Arial')) +
  facet_wrap(~data)
ggsave(paste0('../results/viability_analysis_results/performance_prior_all_models_barplot.eps'),
       device = cairo_ps,
       height = 10.5,
       width = 14,
       units = 'in',
       dpi = 600)

# compare prior and frequency
res <- rbind(prior_performance %>% filter(data=='test') %>% filter(model!='lm') %>% select(model,correlation) %>% mutate(type='prior') %>% unique(),
             frequency_performance %>% filter(data=='test') %>% filter(model!='lm') %>% select(model,correlation) %>% mutate(type='frequency')%>% unique())
res <- rbind(res,
             performance_df %>% filter(frequency==0.62) %>% filter(data=='test') %>% 
               filter(model!='lm') %>% select(model,correlation) %>% mutate(type='threshold:0.62')%>% unique() )
ggboxplot(res,x='type',y='correlation',color='type',add='jitter') +
  ylab('pearson`s r')+
  theme(text = element_text(size = 24,family = 'Arial'),axis.title.x = element_blank()) +
  geom_hline(yintercept = 0,size=1,linetype='dashed',color='black') + 
  stat_compare_means(comparisons = list(c('prior','frequency'),c('prior','threshold:0.62')),method='wilcox.test')
ggsave(paste0('../results/viability_analysis_results/compate_input_types_across_all_models.eps'),
       device = cairo_ps,
       height = 10.5,
       width = 14,
       units = 'in',
       dpi = 600)

# show LASSO, SVM and XGboost only
res <- rbind(prior_performance %>% filter(model=='LASSO') %>% mutate(type='prior') %>% unique(),
             frequency_performance %>% filter(model=='LASSO') %>% mutate(type='frequency')%>% unique())
res <- rbind(res,
             performance_df %>% filter(frequency==0.62) %>% filter(model=='lasso') %>% unique() %>%
               select(all_of(colnames(frequency_performance)))%>% mutate(type='threshold:0.62'))
p1 <- ggscatter(res %>% filter(data=='test'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 8) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('LASSO')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=22,family='Arial'),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~type)
# geom_point(data = res_test %>% filter(data=='test'),aes(x=predicted,y=EC50),color='orange',shape=17,size = 4)
p2 <- ggscatter(res %>% filter(data=='shuffled'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 8) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('LASSO with shuffled data')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=22,family='Arial'),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~type)
p <- p1/p2
print(p)
ggsave(paste0('../results/viability_analysis_results/LASSO_performance_all.eps'),
       device = cairo_ps,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)

#SVMs
res <- rbind(prior_performance %>% filter(model==toupper('svmLinear')) %>% mutate(type='prior') %>% unique(),
             frequency_performance %>% filter(model==toupper('svmLinear')) %>% mutate(type='frequency')%>% unique())
res <- rbind(res,
             performance_df %>% filter(frequency==0.62) %>% filter(model=='svmLinear') %>% unique() %>%
               select(all_of(colnames(frequency_performance)))%>% mutate(type='threshold:0.62'))
p1 <- ggscatter(res %>% filter(data=='test'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 8) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('SVM regression')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=22,family='Arial'),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~type)
# geom_point(data = res_test %>% filter(data=='test'),aes(x=predicted,y=EC50),color='orange',shape=17,size = 4)
p2 <- ggscatter(res %>% filter(data=='shuffled'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 8) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('SVM regression with shuffled data')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=22,family='Arial'),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~type)
p <- p1/p2
print(p)
ggsave(paste0('../results/viability_analysis_results/SVM_performance_all.eps'),
       device = cairo_ps,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)

# XGboost tree
res <- rbind(prior_performance %>% filter(model==toupper('xgbTree')) %>% mutate(type='prior') %>% unique(),
             frequency_performance %>% filter(model==toupper('xgbTree')) %>% mutate(type='frequency')%>% unique())
res <- rbind(res,
             performance_df %>% filter(frequency==0.62) %>% filter(model=='xgbTree') %>% unique() %>%
               select(all_of(colnames(frequency_performance)))%>% mutate(type='threshold:0.62'))
p1 <- ggscatter(res %>% filter(data=='test'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 8) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('XGBoost decision trees')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=22,family='Arial'),
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~type)
# geom_point(data = res_test %>% filter(data=='test'),aes(x=predicted,y=EC50),color='orange',shape=17,size = 4)
p2 <- ggscatter(res %>% filter(data=='shuffled'),x='predicted',y='EC50',cor.coef = T,rug = T,cor.coef.size = 8) +
  xlab('Predicted value') + ylab('True value')+ ggtitle('XGBoost decision trees with shuffled data')+
  geom_abline(intercept = 0,slope=1,linetype=2,color='red',linewidth=2)+
  theme(text=element_text(size=22,family='Arial'),
        plot.title = element_text(hjust = 0.5))+
  facet_wrap(~type)
p <- p1/p2
print(p)
ggsave(paste0('../results/viability_analysis_results/XGBoost_performance_all.eps'),
       device = cairo_ps,
       height = 12,
       width = 16,
       units = 'in',
       dpi = 600)
