library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Visualize Delta1 vs TF activity and with validation correlation------------------------------------------
base_cell_performance <- data.table::fread('../../results/case_study/TrainEnsemblePerformance.csv')
colnames(base_cell_performance)[1] <- 'TF'
Delta <- data.table::fread('../../results/case_study/DeltaTF1.csv',header=T)
# Load TF activities
TFoutput <- read.delim('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv') %>% filter(X %in% Delta$V1) %>%
  column_to_rownames('X') %>% rownames_to_column('sample')
gc()
TFoutput <- TFoutput %>% gather('TF','activity',-sample) 
Delta <- Delta  %>% column_to_rownames('V1') %>%  rownames_to_column('sample') %>% gather('TF','delta',-sample)
# merge everything
df <- left_join(TFoutput,Delta,by=c('sample','TF'))
df <- left_join(df,base_cell_performance)
# Load conditions to get rid of DMSO
conditions <- data.table::fread('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_A375-conditions_drugs.tsv',sep = "\t") %>% column_to_rownames('V1')
conditions <- conditions %>% rownames_to_column('sample') %>% gather('drug','value',-sample) %>% filter(value>0) %>%
  select(-value) %>% unique()
conditions <- conditions %>% filter(sample %in% df$sample) %>% filter(drug!='CS(C)=O')
annotation <- read.delim('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv') %>% dplyr::select(c('TF'='code'),name)
annotation <- annotation %>% filter(TF %in% df$TF)
df <- left_join(df,annotation)
df <- left_join(df,conditions)
df <- df %>% filter(!is.na(drug))
xlsx::write.xlsx2(df,'../../results/case_study/A375_ScoredOffTargetsCaseStudies.xlsx')
interestingTFs = df %>% filter(!is.na(drug))  %>% 
  filter(r>=0.5) %>% 
  filter((delta>0.23 | delta<(-0.23)) & (activity>0.75 | activity<0.25))
interestingTFs <- interestingTFs %>% group_by(TF) %>% mutate(max_score=max(r)) %>%
  ungroup() %>% mutate(keep=ifelse(max_score==r,TRUE,FALSE))
interestingTFs <- interestingTFs %>% filter(keep==TRUE)
interestingTFs <- interestingTFs %>% group_by(TF) %>% mutate(max_score=max(abs(delta))) %>%
  ungroup() %>% mutate(keep=ifelse(max_score==abs(delta),TRUE,FALSE))
interestingTFs <- interestingTFs %>% filter(keep==TRUE)
drugs <- unique(interestingTFs$drug)
p <- ggplot(df %>% filter(!is.na(drug)),aes(x=-delta,y=activity,color=r)) +
  geom_point() + 
  scale_colour_gradient2(low = "blue",mid="white" ,high = "red",midpoint = 0.4,limits=c(min(0,min(df$r)-0.05),1)) + xlab(expression(Delta*"TF")) + ylab('DoRothEA inferred TF activity') + 
  geom_vline(xintercept = -0.25,linetype = 'dashed') + geom_vline(xintercept = 0.25,linetype = 'dashed')+
  geom_hline(yintercept = 0.25,linetype = 'dashed') + geom_hline(yintercept = 0.75,linetype = 'dashed')+
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4),limits = c(-0.4,0.45))+
  ggtitle('Off-target effects in A375 cell line')+
  theme_pubr() + 
  theme(text = element_text(size=24),legend.position = 'right',plot.title = element_text(hjust = 0.5))+
  geom_label_repel(data=interestingTFs,
                   aes(label=name),
                   box.padding   = 0.75, 
                   point.padding = 0.5,
                   max.overlaps = 50,
                   segment.color = 'grey50',
                   size = 7)
print(p)
ggsave('../figures/figure3B_case_study.eps',
       plot = p,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)
data.table::fwrite(interestingTFs,
                   '../results/A375_ensembles/interestingSamples_A375_case_study.csv',
                   row.names = T)
