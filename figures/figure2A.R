library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(viridis)
library(patchwork)

### Load global errors of all drugs
global_errors <- data.table::fread('../results/A375_ensembles/all_drugs_global_errors.csv') %>% select(-V1)
global_errors <- global_errors %>% gather('drug','error',-grad_threshold)
no_thresh <- length(unique(global_errors$grad_threshold))
max_thresh <- max(global_errors$grad_threshold)
global_errors <- global_errors %>% group_by(grad_threshold) %>% mutate(mean_error = mean(error)) %>%
  mutate(sd_error = sd(error)) %>% mutate(max_error = max(error)) %>% mutate(min_error = min(error)) %>% ungroup()
max_global <- global_errors$mean_error[which(global_errors$grad_threshold==max_thresh)[1]]

unmasked_global_errors <- data.table::fread('../results/A375_ensembles/all_drugs_unmasked_global_errors.csv',skip=1)
colnames(unmasked_global_errors) <- c('drug','unmasked_error')

global_errors <- left_join(global_errors,unmasked_global_errors)
global_errors <- global_errors %>% group_by(grad_threshold) %>% 
  mutate(mean_unmasked_error = mean(unmasked_error)) %>% ungroup()

### Select mean threshold
global_errors_mean <- global_errors %>% select(grad_threshold,mean_error,mean_unmasked_error) %>% unique()
global_errors_mean <- global_errors_mean %>% mutate(percentage_change = (mean_error - mean_unmasked_error)/(max_global-mean_unmasked_error))
if (any(global_errors_mean$percentage_change <= 0.25)) {
  inds <- which(global_errors_mean$percentage_change <= 0.25 & global_errors_mean$percentage_change>=0)
  global_thresh_selection <- inds[length(inds)]
} else{
  global_thresh_selection = 0
}

## Visualize global error vs grad threshold
p <- ggplot(global_errors,aes(x=grad_threshold,y=mean_error)) +
  geom_line(color='black',linewidth=1)+
  geom_ribbon(aes(ymin = mean_error - sd_error,ymax=mean_error + sd_error),alpha=0.1)+
  geom_segment(aes(x = min(global_errors_mean$grad_threshold), xend = global_errors_mean$grad_threshold[global_thresh_selection],
               y = global_errors_mean$mean_error[global_thresh_selection],yend = global_errors_mean$mean_error[global_thresh_selection]),
               linetype='dashed',linewidth=0.75,color='red',arrow = arrow(type = "closed",length = unit(2,"mm")))+
  geom_segment(aes(x = global_errors_mean$grad_threshold[global_thresh_selection], xend = global_errors_mean$grad_threshold[global_thresh_selection],
                   y = min(mean_error- sd_error),yend = global_errors_mean$mean_error[global_thresh_selection]),
               linetype='dashed',linewidth=0.75,color='red',arrow = arrow(type = "closed",length = unit(2,"mm"),ends = 'first'))+
  annotate('text',x=2.5e-02,y=1.1*global_errors_mean$mean_error[global_thresh_selection],
           label='25% error increase threshold',size=10)+
  scale_x_log10(breaks = c(1e-03,1e-02,1e-01,1,1e+01,1e+02,1e+03))+
  ggtitle('Global error of the model when masking interactions')+
  xlab('absolute integrated gradient score cut-off') +
  ylab('error')+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(text = element_text(family = 'Arial',size=24),
        plot.title = element_text(family = 'Arial',size=26,hjust=0.5))
print(p)
png('../figures/figure2A.png',units = 'in',width = 12,height = 8,res = 600)
print(p)
dev.off()

ggsave('../figures/figure2A.eps',
       plot = p,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)
