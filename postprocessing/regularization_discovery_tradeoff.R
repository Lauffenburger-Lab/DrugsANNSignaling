library(tidyverse)
library(gg3D)
library(ggpubr)
library(patchwork)
library(xlsx)

# Load regularization results
regularization <- read.xlsx('../results/regularizationTuneInteractionResults/hitDiscoveryResults.xlsx',
                            sheetIndex = 1)[,1:13]
regularization$threshold <- factor(regularization$threshold,
                                   levels = unique(regularization$threshold))
regularization_errorbased <- data.table::fread('../results/regularizationTuneInteractionResults/error_based_inference/all_regularizations_interactions_drugs.csv')
regularization_errorbased <- regularization_errorbased %>% select(-V1)
regularization_errorbased <- regularization_errorbased %>% mutate(tp=1*(`Prior knowledge`==Inferred & Inferred=='Interaction')) %>% 
  mutate(fp=1*(`Prior knowledge`=='No interaction' & Inferred=='Interaction')) %>%
  mutate(fn=1*(`Prior knowledge`=='Interaction' & Inferred=='No interaction')) %>%
  mutate(tn=1*(`Prior knowledge`==Inferred & Inferred=='No interaction'))
regularization_errorbased <- regularization_errorbased %>% group_by(lamda) %>% 
  mutate(tp = sum(tp)) %>%
  mutate(tn = sum(tn)) %>%
  mutate(fp = sum(fp)) %>%
  mutate(fn = sum(fn)) %>%
  ungroup() %>%
  select(lamda,tp,tn,fp,fn) %>%
  unique()
regularization_errorbased <- regularization_errorbased %>% mutate(sensitivity = tp/(tp+fn)) %>% mutate(specificity = tn/(tn+fp)) %>%
  mutate(Gmean = sqrt(sensitivity*specificity)) %>% mutate(NDR=(fp/(fp+tp)))

regularization <- regularization %>% mutate(NDR=new_hits/(new_hits+tp)) %>% group_by(lamda) %>%
  mutate(meanG=mean(G)) %>% mutate(sdG=sd(G)/sqrt(n_distinct(threshold))) %>%
  mutate(meanNDR=mean(NDR)) %>% mutate(sdNDR=sd(NDR)/sqrt(n_distinct(threshold))) %>%
  ungroup()
p2 <- ggplot(rbind(regularization %>% select(lamda,c('value'='meanG'),c('valueSE'='sdG')) %>% mutate(metric='G-mean') %>% unique(),
                   regularization %>% select(lamda,c('value'='meanNDR'),c('valueSE'='sdNDR'))%>% mutate(metric='New Discovery Rate') %>% unique()),
             aes(x=lamda,y=value,color = metric))+
  geom_point() + geom_line(linewidth = 0.75)+ 
  geom_errorbar(aes(ymin=value-valueSE, ymax=value+valueSE), width=.3)+
  scale_color_manual(values =  c("black","red"))+
  scale_x_log10() +
  xlab(expression(lambda))+ylab('value') + 
  ylim(c(0,1))+  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text=element_text(family = 'Arial',size=24),
        axis.title.x =  element_text(family = 'Arial',size=36),
        legend.position = 'right',
        legend.text = element_text(family = 'Arial',size=24),
        plot.title = element_text(hjust = 0.25,size=24),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2)
p2_error_derived <- ggplot(rbind(regularization_errorbased %>% select(lamda,c('value'='Gmean')) %>% mutate(metric='G-mean') %>% unique(),
                                 regularization_errorbased %>% select(lamda,c('value'='NDR')) %>% mutate(metric='New Discovery Rate') %>% unique()),
                           aes(x=lamda,value,color=metric)) +
  geom_point() + geom_line(linewidth = 0.75)+ 
  scale_color_manual(values =  c("black","red"))+
  ylim(c(0,1)) + scale_x_log10() +
  xlab(expression(lambda))+ylab('value') +  theme_pubr(base_family = 'Arial',base_size = 24)+ 
  theme(text=element_text(family = 'Arial',size=24),
        axis.title.x =  element_text(family = 'Arial',size=36),
        legend.text = element_text(family = 'Arial',size=24),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.25,size=24),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2_error_derived)
ggsave('../article_supplementary_info/suppl_figure6B.eps',
       plot = p2_error_derived,
       device = cairo_ps,
       scale = 1,
       width = 16,
       height = 8,
       units = "in",
       dpi = 600)

ggsave('../article_supplementary_info/suppl_figure6A.eps',
       plot = p2,
       device = cairo_ps,
       scale = 1,
       width = 16,
       height = 8,
       units = "in",
       dpi = 600)

### Get also plot for optimal G-mean
regularization_optimal <- read.xlsx('../results/regularizationTuneInteractionResults/hitDiscoveryResults.xlsx',
                            sheetIndex = 1)
regularization_optimal <- regularization_optimal %>% mutate(NDR=new_hits.1/(new_hits.1+tp.1)) %>% select(lamda,optimalG,NDR) %>% unique()
p2_optimal <- ggplot(rbind(regularization_optimal %>% select(lamda,c('value'='optimalG')) %>% mutate(metric='G-mean') %>% unique(),
                           regularization_optimal %>% select(lamda,c('value'='NDR'))%>% mutate(metric='New Discovery Rate') %>% unique()),
                     aes(x=lamda,y=value,color = metric))+
  geom_point() + geom_line(linewidth = 0.75)+ 
  scale_color_manual(values =  c("black","red"))+
  scale_x_log10() +
  xlab(expression(lambda))+ylab('value') + 
  ylim(c(0,1))+  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text=element_text(family = 'Arial',size=24),
        legend.position = 'right',
        axis.title.x =  element_text(family = 'Arial',size=36),
        legend.text = element_text(family = 'Arial',size=24),
        plot.title = element_text(hjust = 0.25,size=24),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2_optimal)
ggsave('../figures/figure1E',
       plot = p2_optimal,
       device = cairo_ps,
       scale = 1,
       width = 16,
       height = 8,
       units = "in",
       dpi = 600)