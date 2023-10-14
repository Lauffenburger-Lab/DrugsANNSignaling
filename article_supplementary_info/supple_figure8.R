library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data for the plot--------
a375_scores <-  data.table::fread('../results/A375_ensembles/grads_vs_weights_A375.csv') %>% select(-V1)
a549_scores <-  data.table::fread('../results/A549_ensembles/grads_vs_weights_A549.csv')%>% select(-V1)
vcap_scores <- data.table::fread('../results/FinalEnsemble/grads_vs_weights_VCAP.csv')%>% select(-V1)
all_scores <- rbind(a375_scores,a549_scores,vcap_scores)
all_scores <- all_scores %>% group_by(cell,model_no) %>% 
  mutate(R2=MLmetrics::R2_Score(grad_scores,weight_scores)) %>% ungroup()

# Create plots-------------------
ggscatter(all_scores,x='grad_scores',y='weight_scores',cor.coef = T,cor.coef.size = 10,
          shape='cell',color = 'model_no') + labs(color='model #')+
  xlab('integrated gradient score') + ylab('weight-derived score')+
  ggtitle('Integrated gradient scores VS weight-derived scores')+
  theme(text=element_text(family = 'Arial',size = 24),
        legend.text = element_text(size=20),
        legend.position = 'right',
        title = element_text(size=22))

ggsave('grads_vs_weights.png',
       scale = 1,
       width = 12,
       height = 12,
       units = "in",
       dpi = 600)
