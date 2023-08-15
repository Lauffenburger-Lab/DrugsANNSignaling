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
library(ggrepel)
library(factoextra)


# Visualize Delta1 vs TF activity and with validation correlation------------------------------------------
performance <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/meanCorrPerTFEnsembleVal_lamda6.csv',
                                 header=T)
performance <- performance %>% dplyr::select(-model) %>% unique()
performance <- performance %>% group_by(TF) %>% mutate(mean_r=mean(r)) %>%
  ungroup() %>% dplyr::select(-cell,-r) %>% unique()
base_cell_performance <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/A375TrainEnsemblePerformance.csv')
colnames(base_cell_performance)[1] <- 'TF'
Delta <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/DeltaTF1.csv',header=T)
# Load TF activities
TFoutput <- read.delim('Model/data/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv') %>% filter(X %in% Delta$V1) %>%
  column_to_rownames('X') %>% rownames_to_column('sample')
gc()
TFoutput <- TFoutput %>% gather('TF','activity',-sample)  #%>% group_by(TF) %>% 
#  mutate(mean_activity = mean(activity)) %>%  select(-activity) %>% unique()
Delta <- Delta  %>% column_to_rownames('V1') %>%  rownames_to_column('sample') %>% gather('TF','delta',-sample)
#Delta <- Delta %>% group_by(TF) %>% mutate(mean_delta = mean(delta)) %>%  select(-delta) %>% unique()
# merge everything
df <- left_join(TFoutput,Delta,by=c('sample','TF'))
df <- left_join(df,performance)
df <- left_join(df,base_cell_performance)
df <- df %>% mutate(score=0.5*(mean_r+r))
# Load conditions to get rid of DMSO
conditions <- data.table::fread('Model/data/L1000_lvl3_allcells-conditions_drugs.tsv',sep = "\t") %>% column_to_rownames('V1')
conditions <- conditions %>% rownames_to_column('sample') %>% gather('drug','value',-sample) %>% filter(value>0) %>%
  select(-value) %>% unique()
conditions <- conditions %>% filter(sample %in% df$sample) %>% filter(drug!='CS(C)=O')
annotation <- read.delim('Model/data/l1000_lvl3_withsignor-Annotation.tsv') %>% dplyr::select(c('TF'='code'),name)
annotation <- annotation %>% filter(TF %in% df$TF)
df <- left_join(df,annotation)
df <- left_join(df,conditions)
df <- df %>% filter(!is.na(drug))

interestingTFs = df %>% filter(!is.na(drug))  %>% 
  filter(mean_r>0.4 & score>=0.5) %>% 
  filter((delta>0.23 | delta<(-0.23)) & (activity>0.75 | activity<0.25))
interestingTFs <- interestingTFs %>% group_by(TF) %>% mutate(max_score=max(score)) %>%
  ungroup() %>% mutate(keep=ifelse(max_score==score,TRUE,FALSE))
interestingTFs <- interestingTFs %>% filter(keep==TRUE)
interestingTFs <- interestingTFs %>% group_by(TF) %>% mutate(max_score=max(abs(delta))) %>%
  ungroup() %>% mutate(keep=ifelse(max_score==abs(delta),TRUE,FALSE))
interestingTFs <- interestingTFs %>% filter(keep==TRUE)
#samples <- NULL
#for (tf in unique(interestingTFs$TF)){
#  tmp <- interestingTFs %>% filter(TF=='tf')
#  samples <- c(samples,sample(unique(interestingTFs$sample),2))
#}
#interestingTFs <- interestingTFs %>% filter(sample %in% samples)
drugs <- unique(interestingTFs$drug)
p <- ggplot(df %>% filter(!is.na(drug)),aes(x=-delta,y=activity,color=score)) +
  geom_point() + 
  scale_colour_gradient2(low = "blue",mid="white" ,high = "red",midpoint = 0.4,limits=c(0,0.72)) + xlab(expression("TF"~Delta*"activity")) + ylab('TF activity') + 
  geom_vline(xintercept = -0.25,linetype = 'dashed') + geom_vline(xintercept = 0.25,linetype = 'dashed')+
  geom_hline(yintercept = 0.25,linetype = 'dashed') + geom_hline(yintercept = 0.75,linetype = 'dashed')+
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4),limits = c(-0.4,0.45))+
  theme_pubr() + 
  theme(text = element_text(size=20),legend.position = 'right')+
  geom_label_repel(data=interestingTFs,
            aes(label=name),
            box.padding   = 0.75, 
            point.padding = 0.5,
            max.overlaps = 50,
            segment.color = 'grey50',
            size = 7)
print(p)
png('Model/CVL1000_Paper/A375_ensembles/interestingTFs.png',units = 'in',width = 12,height = 9,res=600)
print(p)
dev.off()

ggsave('../MIT/LauffenburgerLab/drugLembasPaper/A375_interestingTFs.eps',
       plot = p,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)

data.table::fwrite(interestingTFs,
                   'Model/CVL1000_Paper/A375_ensembles/interestingSamples.csv',
                   row.names = T)

### Also plot delta vs drug-targer score
maxDelta <- df%>% group_by(drug)%>% mutate(max_delta = max(abs(delta)))%>% ungroup()
maxDelta <- maxDelta %>% select(sample,drug,max_delta) %>% unique()

grad_scores <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/InteractionScores/l1000_modeltype4_lamda6_A375_interactionScoresEnsembled.csv')
colnames(grad_scores)[1] <- 'drug'
grad_scores <- grad_scores %>% filter(drug %in% maxDelta$drug)
grad_scores <- grad_scores %>% gather('variable','InteractionScore',-drug)
merged_interactions <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/InteractionScores/l1000_modeltype4_lamda6_A375_mergedInteractions_ROC.csv',
                                         header = T)
merged_interactions <- merged_interactions %>% select(-V1)
merged_interactions <- merged_interactions %>% filter(drug %in% maxDelta$drug)
# merged_interactions <- merged_interactions %>% filter(`Prior knowledge`!=Inferred)
grad_scores <- left_join(grad_scores,merged_interactions)
grad_scores <- grad_scores %>% filter(!is.na(Inferred))
grad_scores <- grad_scores %>% group_by(drug) %>% mutate(max_score = max(abs(InteractionScore)))
grad_scores <- grad_scores %>% select(drug,max_score) %>% unique()

maxDelta <- maxDelta %>% filter(drug %in% grad_scores$drug)
maxDelta <- left_join(maxDelta,grad_scores)
maxDelta <- maxDelta %>% mutate(log_score = log10(max_score+1))

ggscatter(maxDelta,x='max_score',y='max_delta',rug = TRUE,
          alpha = 0.5,size=1,
          cor.coef=T,cor.coef.size = 5,cor.coef.coord = c(0.05, 0.9)) + 
  ggtitle('Off-target effect in activity VS gradient score of drug-target interaction') +
  xlab('Max absolute grad score') + ylab('Max absolute ΔTF')+
  theme(plot.title = element_text(hjust = 0.5,size=15))

