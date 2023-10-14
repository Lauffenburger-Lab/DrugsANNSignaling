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
  mutate(Gmean = sqrt(sensitivity*specificity)) %>% mutate(NDR=(fp/(fp+tp))) %>% mutate(`True Positive Rate`=sensitivity)

regularization <- regularization %>% mutate(NDR=new_hits/(new_hits+tp)) %>% 
  mutate(`True Positive Rate`=sensitivity) %>% group_by(lamda) %>%
  mutate(meanG=mean(G)) %>% mutate(sdG=sd(G)/sqrt(n_distinct(threshold))) %>%
  mutate(meanNDR=mean(NDR)) %>% mutate(sdNDR=sd(NDR)/sqrt(n_distinct(threshold))) %>%
  mutate(meanTPR=mean(`True Positive Rate`)) %>% mutate(sdTPR=sd(`True Positive Rate`)/sqrt(n_distinct(threshold))) %>%
  ungroup()
p2 <- ggplot(rbind(regularization %>% select(lamda,c('value'='meanG'),c('valueSE'='sdG')) %>% mutate(metric='G-mean = sqrt(sensitivity*specificity)') %>% unique(),
                   regularization %>% select(lamda,c('value'='meanNDR'),c('valueSE'='sdNDR'))%>% mutate(metric='New Discovery Rate') %>% unique(),
                   regularization %>% select(lamda,c('value'='meanTPR'),c('valueSE'='sdTPR'))%>% mutate(metric='True Positive Rate') %>% unique()),
             aes(x=lamda,y=value,color = metric))+
  geom_point() + geom_line(linewidth = 0.75)+ 
  geom_errorbar(aes(ymin=value-valueSE, ymax=value+valueSE), width=.3)+
  scale_color_manual(values =  c("black","red","blue"))+
  scale_x_log10() +
  xlab(expression(lambda))+ylab('value') + 
  ylim(c(0,1))+  theme_pubr(base_family = 'Arial',base_size = 28)+
  theme(text=element_text(family = 'Arial',size=28),
        axis.title.x =  element_text(family = 'Arial',size=40),
        legend.position = 'right',
        legend.text = element_text(family = 'Arial',size=26),
        plot.title = element_text(hjust = 0.25,size=28),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2)
p2_error_derived <- ggplot(rbind(regularization_errorbased %>% select(lamda,c('value'='Gmean')) %>% mutate(metric='G-mean = sqrt(sensitivity*specificity)') %>% unique(),
                                 regularization_errorbased %>% select(lamda,c('value'='NDR')) %>% mutate(metric='New Discovery Rate') %>% unique(),
                                 regularization_errorbased %>% select(lamda,c('value'='True Positive Rate')) %>% mutate(metric='True Positive Rate') %>% unique()),
                           aes(x=lamda,value,color=metric)) +
  geom_point() + geom_line(linewidth = 0.75)+ 
  scale_color_manual(values =  c("black","red","blue"))+
  ylim(c(0,1)) + scale_x_log10() +
  xlab(expression(lambda))+ylab('value') +  theme_pubr(base_family = 'Arial',base_size = 28)+ 
  theme(text=element_text(family = 'Arial',size=28),
        axis.title.x =  element_text(family = 'Arial',size=40),
        legend.text = element_text(family = 'Arial',size=26),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.25,size=28),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2_error_derived)
ggsave('../article_supplementary_info/suppl_figure7B.eps',
       plot = p2_error_derived,
       device = cairo_ps,
       scale = 1,
       width = 16,
       height = 8,
       units = "in",
       dpi = 600)

ggsave('../article_supplementary_info/suppl_figure7A.eps',
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
regularization_optimal <- regularization_optimal %>% mutate(NDR=new_hits.1/(new_hits.1+tp.1)) %>% mutate(TPR=(tp.1/(tp.1+missed_hits.1))) %>% 
  select(lamda,optimalG,NDR,TPR) %>% 
  unique()
p2_optimal <- ggplot(rbind(regularization_optimal %>% select(lamda,c('value'='optimalG')) %>% mutate(metric='G-mean = sqrt(sensitivity*specificity)') %>% unique(),
                           regularization_optimal %>% select(lamda,c('value'='NDR'))%>% mutate(metric='New Discovery Rate') %>% unique(),
                           regularization_optimal %>% select(lamda,c('value'='TPR'))%>% mutate(metric='True Positive Rate') %>% unique()),
                     aes(x=lamda,y=value,color = metric))+
  geom_point() + geom_line(linewidth = 0.75)+ 
  scale_color_manual(values =  c("black","red",'blue'))+
  scale_x_log10() +
  xlab(expression(lambda))+ylab('value') + 
  ylim(c(0,1))+  theme_pubr(base_family = 'Arial',base_size = 28)+
  theme(text=element_text(family = 'Arial',size=28),
        legend.position = 'right',
        axis.title.x =  element_text(family = 'Arial',size=40),
        legend.text = element_text(family = 'Arial',size=26),
        plot.title = element_text(hjust = 0.25,size=28),
        panel.grid.major = element_line(colour="black",
                                        linetype = 'dashed', 
                                        linewidth=0.25))
print(p2_optimal)
ggsave('../figures/figure1E.eps',
       plot = p2_optimal,
       device = cairo_ps,
       scale = 1,
       width = 18,
       height = 9,
       units = "in",
       dpi = 600)
