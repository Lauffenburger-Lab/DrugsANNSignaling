library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(viridis)
library(patchwork)
library(caret)

### Give as input drug of interest and TF of interest
TF = "Q08050"
TF_gene = "FOXM1"
drug_smile = "C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13"
drug_name = "lestaurtinib"
sample = "CPC014_A375_6H:BRD-K23192422-001-01-1:10"

### Load global errors of all drugs
global_errors <- data.table::fread('../results/A375_ensembles/allmodels_all_drugs_global_errors.csv',header = T) 
colnames(global_errors)[1] <- 'model_no'
global_errors <- global_errors %>% filter(drug == drug_smile)
global_errors <- global_errors %>% select(-drug)
global_errors <- global_errors %>% gather('grad_threshold','error',-model_no)
global_errors$grad_threshold <- as.numeric(global_errors$grad_threshold)
unmasked_global_errors <- data.table::fread('../results/A375_ensembles/allmodels_all_drugs_unmasked_global_errors.csv',header=T)
colnames(unmasked_global_errors)[1] <- 'drug'
unmasked_global_errors <- unmasked_global_errors %>% gather('model_no','unmasked_error',-drug)
unmasked_global_errors <- unmasked_global_errors %>% filter(drug == drug_smile)
unmasked_global_errors <- unmasked_global_errors %>% select(-drug)
unmasked_global_errors$model_no <- as.numeric(unmasked_global_errors$model_no)
global_errors <- left_join(global_errors,unmasked_global_errors)
global_errors <- global_errors %>% group_by(grad_threshold) %>% mutate(mean_umnasked = mean(unmasked_error)) %>%
  mutate(mean_error = mean(error)) %>% mutate(sd_error = sd(error)) %>%
  mutate(max_error = max(error)) %>% mutate(min_error = min(error)) %>% ungroup() %>%
  group_by(model_no) %>% mutate(model_max_err = max(error)) %>% ungroup()
no_thresh <- length(unique(global_errors$grad_threshold))
max_thresh <- max(global_errors$grad_threshold)
max_global <- global_errors$mean_error[which(global_errors$grad_threshold==max_thresh)[1]]

### Load TF error of all drugs
TF_errors <- data.table::fread(paste0('../results/A375_ensembles/allmodels_all_drugs_',TF_gene,'_errors.csv'),header = T) 
colnames(TF_errors)[1] <- 'model_no'
TF_errors <- TF_errors %>% filter(drug == drug_smile)
TF_errors <- TF_errors %>% select(-drug)
TF_errors <- TF_errors %>% gather('grad_threshold','error',-model_no)
TF_errors$grad_threshold <- as.numeric(TF_errors$grad_threshold)
unmasked_TF_errors <- data.table::fread(paste0('../results/A375_ensembles/allmodels_all_drugs_unmasked_',TF_gene,'_errors.csv'),header=T)
colnames(unmasked_TF_errors)[1] <- 'drug'
unmasked_TF_errors <- unmasked_TF_errors %>% gather('model_no','unmasked_error',-drug)
unmasked_TF_errors <- unmasked_TF_errors %>% filter(drug == drug_smile)
unmasked_TF_errors <- unmasked_TF_errors %>% select(-drug)
unmasked_TF_errors$model_no <- as.numeric(unmasked_TF_errors$model_no)
TF_errors <- left_join(TF_errors,unmasked_TF_errors)
TF_errors <- TF_errors %>% group_by(grad_threshold) %>% mutate(mean_umnasked = mean(unmasked_error)) %>%
  mutate(mean_error = mean(error)) %>% mutate(sd_error = sd(error)) %>%
  mutate(max_error = max(error)) %>% mutate(min_error = min(error)) %>% ungroup() %>%
  group_by(model_no) %>% mutate(model_max_err = max(error)) %>% ungroup()
no_thresh <- length(unique(TF_errors$grad_threshold))
max_TF_thresh <- max(TF_errors$grad_threshold)
max_TF <- TF_errors$mean_error[which(TF_errors$grad_threshold==max_TF_thresh)[1]]


### Select mean threshold
global_errors_mean <- global_errors %>% select(grad_threshold,mean_error,mean_umnasked) %>% unique()
global_errors_mean <- global_errors_mean %>% mutate(percentage_change = (mean_error - mean_umnasked)/(max_global-mean_umnasked))
if (any(global_errors_mean$percentage_change <= 0.25)) {
  inds <- which(global_errors_mean$percentage_change <= 0.25 & global_errors_mean$percentage_change>=0)
  global_thresh_selection <- inds[length(inds)]
} else{
  global_thresh_selection = 0
}
TF_errors_mean <- TF_errors %>% select(grad_threshold,mean_error,mean_umnasked) %>% unique()
TF_errors_mean <- TF_errors_mean %>% mutate(percentage_change = (mean_error - mean_umnasked)/(max_TF-mean_umnasked))
if (any(TF_errors_mean$percentage_change <= 0.25)) {
  inds <- which(TF_errors_mean$percentage_change <= 0.25 & TF_errors_mean$percentage_change>=0)
  TF_thresh_selection <- inds[length(inds)]
} else{
  TF_thresh_selection = 0
}

# Combine global and local into one graph
all_errors <- rbind(TF_errors %>% mutate('error type'=TF_gene),
                    global_errors %>% mutate('error type'='Global'))

## Visualize global error vs grad threshold
p <- ggplot(all_errors,aes(x=grad_threshold,y=mean_error,color=`error type`,fill=`error type`)) +
  geom_line(linewidth=1)+
  geom_ribbon(aes(ymin = mean_error - sd_error,ymax=mean_error + sd_error),colour = NA,alpha=0.1)+
  scale_color_manual(values = c('#d97b38','#4878CF'))+
  scale_fill_manual(values = c('#d97b38','#4878CF'))+
  geom_segment(aes(x = min(TF_errors_mean$grad_threshold), xend = TF_errors_mean$grad_threshold[TF_thresh_selection],
                   y = TF_errors_mean$mean_error[TF_thresh_selection],yend = TF_errors_mean$mean_error[TF_thresh_selection]),
               linetype='dashed',linewidth=1,color='#d97b38',arrow = arrow(type = "closed",length = unit(2,"mm")))+
  geom_segment(aes(x = TF_errors_mean$grad_threshold[TF_thresh_selection], xend = TF_errors_mean$grad_threshold[TF_thresh_selection],
                   y = min(mean_error- sd_error),yend = TF_errors_mean$mean_error[TF_thresh_selection]),
               linetype='dashed',linewidth=1,color='#d97b38',arrow = arrow(type = "closed",length = unit(2,"mm"),ends = 'first'))+
  geom_segment(aes(x = min(global_errors_mean$grad_threshold), xend = global_errors_mean$grad_threshold[global_thresh_selection],
                   y = global_errors_mean$mean_error[global_thresh_selection],yend = global_errors_mean$mean_error[global_thresh_selection]),
               linetype='dashed',linewidth=1,color='#4878CF',arrow = arrow(type = "closed",length = unit(2,"mm")))+
  geom_segment(aes(x = global_errors_mean$grad_threshold[global_thresh_selection], xend = global_errors_mean$grad_threshold[global_thresh_selection],
                   y = min(mean_error- sd_error),yend = global_errors_mean$mean_error[global_thresh_selection]),
               linetype='dashed',linewidth=1,color='#4878CF',arrow = arrow(type = "closed",length = unit(2,"mm"),ends = 'first'))+
  # annotate('text',x=3e-03,y=1.2*global_errors_mean$mean_error[global_thresh_selection],label='25% error increase threshold',size=6)+
  scale_x_log10(breaks = c(1e-03,1e-02,1e-01,1,1e+01,1e+02,1e+03))+
  ggtitle(paste0('Removing interactions for ',
                 paste0(toupper(substr(drug_name, 1, 1)), substr(drug_name, 2, nchar(drug_name))),
                 ' in A375 cell line'))+
  xlab('absolute integrated gradient score cut-off') +
  ylab('error')+
  theme_pubr(base_family = 'Arial',base_size = 26)+
  theme(text = element_text(family = 'Arial',size=26),
        plot.title = element_text(family = 'Arial',size=24,hjust=1))
print(p)

png('../figures/figure3C.png',units = 'in',width = 12,height = 8,res = 600)
print(p)
dev.off()

ggsave('../figures/figure3C.eps',
       plot = p,
       device = cairo_ps,
       scale = 1,
       width = 9,
       height = 6,
       units = "in",
       dpi = 600)

### Build average confusion-matrix and bar-plot for Lestaurtinib-------------------------------------------------

no_models <- 50
merged_interactions_all <- data.table::fread('../results/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions_all <- merged_interactions_all %>% filter(drug==drug_smile)
ensemble_table <- NULL
for (i in 1:no_models){
  merged_interactions <- merged_interactions_all %>% filter(model_no == i-1) %>% select(-model_no)
  confusion <- confusionMatrix(data=as.factor(merged_interactions$Inferred), 
                               reference = as.factor(merged_interactions$`Prior knowledge`))
  ensemble_table[[i]] <-  confusion$table
  message(paste0('Finished model ',i))
}
ensemble_mat <- do.call(cbind,ensemble_table)
ensemble_mat <- array(ensemble_mat,c(dim=dim(ensemble_table[[1]]),length(ensemble_table)))
ensemble_mat_mean <- apply(ensemble_mat, c(1,2), mean, na.rm = TRUE)
colnames(ensemble_mat_mean) <- colnames(ensemble_table[[1]])
rownames(ensemble_mat_mean) <- rownames(ensemble_table[[1]])
ensemble_mat_sd <- apply(ensemble_mat, c(1,2), sd, na.rm = TRUE)
colnames(ensemble_mat_sd) <- colnames(ensemble_table[[1]])
rownames(ensemble_mat_sd) <- rownames(ensemble_table[[1]])
ensemble_mat_mean <- melt(ensemble_mat_mean, variable.name = c("Reference", "Prediction"), value.name = "Count")
ensemble_mat_sd <- melt(ensemble_mat_sd, variable.name = c("Reference", "Prediction"), value.name = "Count_sd")
ensemble_final <- left_join(ensemble_mat_mean,ensemble_mat_sd)
plot_confusion <- ggplot(ensemble_final %>% mutate(Count = ifelse(Count>=1000,round(Count,-3),ifelse(Count>=10,round(Count,-1),round(Count,0)))) %>%
                           mutate(Count_se = Count_sd/sqrt(50)) %>%
                           mutate(Count_se = ifelse(Count_se>=1000,round(Count_se,-3),ifelse(Count_se>=10,round(Count_se,-1),round(Count_se,0)))) , 
                         aes(x = Var1, y = Var2, fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue")+
  # scale_fill_gradientn(colours = colorspace::heat_hcl(7),limits = c(0,round(max(ensemble_final$Count),-2)))+
  labs(fill='Count')+
  geom_text(aes(label = paste0('~',paste(Count, Count_se, sep = "\u00B1"))), size = 14) +
  labs(title = paste0("Inferred drug-target interactions for ",
       paste0(toupper(substr(drug_name, 1, 1)), substr(drug_name, 2, nchar(drug_name)))),
       x = "Inferred",
       y = "Prior knowledge") +
  theme_minimal() +
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24, face = "bold"),
        legend.text = element_text(family = 'Arial',size = 20),
        plot.title = element_text(family = 'Arial',size = 24, hjust = 0.5))
print(plot_confusion)

ggsave('../figures/A375_ensembled_drugtarget_lestrautinib.eps',
       plot = plot_confusion,
       device = cairo_ps,
       scale = 1,
       width = 9,
       height = 6,
       units = "in",
       dpi = 600)
