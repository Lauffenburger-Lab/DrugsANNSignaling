library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load data for the plot--------
files_weights <- list.files('../../results/case_study/InteractionScores/',full.names = TRUE)
files_weights <- files_weights[grepl('.csv',files_weights)]
files_weights <- files_weights[grepl('interactionScores_',files_weights)]
vals <- sapply(files_weights,FUN=function(x){strsplit(x,'_')[[1]][3]})
vals <- unname(vals)
vals <- as.numeric(unname(sapply(vals,FUN=function(x){substr(x,start = 1,stop = nchar(x)-4)[[1]][1]})))
files_weights <- files_weights[order(vals)]
# files_weights <- files_weights[1:9]
df_weights <- data.frame()
i <- 0
for (file in files_weights){
  df <- data.table::fread(file)
  colnames(df)[1] <- 'drug'
  df <- df %>% gather('target','weight',-drug) %>% mutate(model_no = i)
  df_weights <- rbind(df_weights,df)
  i <- i+1
}

files_grads <- list.files('../../results/case_study/interactions_gpu/',full.names = FALSE)
files_grads <- files_grads[grepl('.csv',files_grads)]
files_grads <- files_grads[grepl('interactionScores_',files_grads)]
vals <- sapply(files_grads,FUN=function(x){strsplit(x,'_')[[1]][2]})
vals <- unname(vals)
vals <- as.numeric(unname(sapply(vals,FUN=function(x){substr(x,start = 1,stop = nchar(x)-4)[[1]][1]})))
files_grads <- files_grads[order(vals)]
# files_grads <- files_grads[1:9]
df_score <- data.frame()
i <- 0
for (file in files_grads){
  df <- data.table::fread(paste0('../../results/case_study/interactions_gpu/',file))
  colnames(df)[1] <- 'drug'
  df <- df %>% gather('target','score',-drug)%>% mutate(model_no = i)
  df_score <- rbind(df_score,df)
  i <- i+1
}

res <- left_join(df_score,df_weights)


# ### The weights were calculated by passing an identity matrix a input
# files_weights_v2 <- list.files('../../results/case_study/InteractionScores_v2/',full.names = TRUE)
# files_weights_v2 <- files_weights_v2[grepl('.csv',files_weights_v2)]
# files_weights_v2 <- files_weights_v2[grepl('interactionScores_',files_weights_v2)]
# vals_v2 <- sapply(files_weights_v2,FUN=function(x){strsplit(x,'_')[[1]][5]})
# vals_v2 <- unname(vals_v2)
# vals_v2 <- as.numeric(unname(sapply(vals_v2,FUN=function(x){substr(x,start = 1,stop = nchar(x)-4)[[1]][1]})))
# files_weights_v2 <- files_weights_v2[order(vals_v2)]
# df_weights_v2 <- data.frame()
# i <- 0
# for (file in files_weights_v2){
#   df <- data.table::fread(file)
#   colnames(df)[1] <- 'drug'
#   df <- df %>% gather('target','weight',-drug) %>% mutate(model_no = i)
#   df_weights_v2 <- rbind(df_weights_v2,df)
#   i <- i+1
# }
# 
# res <- left_join(df_weights,df_weights_v2,by=c('drug','target','model_no'))





# Create plots-------------------
ggscatter(res,x='score',y='weight',color = 'model_no',cor.coef = T,cor.coef.size = 10) + labs(color='model #')+
  xlab('integrated gradient score') + ylab('weight-derived score')+
  ggtitle('Integrated gradient scores VS weight-derived scores')+
  theme(text=element_text(family = 'Arial',size = 24),
        legend.text = element_text(size=20),
        legend.position = 'right',
        title = element_text(size=22))

ggsave('grads_vs_weights_gpu_10steps.png',
       scale = 1,
       width = 12,
       height = 12,
       units = "in",
       dpi = 600)
