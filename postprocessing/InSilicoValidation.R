# Use Series GSE31534 from GEO, with siRNAs experiments.
# Compare CDK2 KOs with FOXM1 KOs.
# Compare also with the know lestrautinib Kos targets: none of them in the dataset we use 
# Calculate the gene expression and also the activity of FOXM1

library(tidyverse)
library(hgu133a.db)
library(limma)
library(GEOquery)
library(affy)
library(dorothea)
library(ggplot2)
library(ggpubr)

# Load data
geo<-  getGEO(GEO = 'GSE31534', AnnotGPL = TRUE)
raw <- geo[["GSE31534_series_matrix.txt.gz"]]@assayData[['exprs']]
pheno <- geo[["GSE31534_series_matrix.txt.gz"]]@phenoData@data
list_files <- list.files('GSE31534/')
raw <- ReadAffy(celfile.path = "GSE31534",filenames = list_files,
                phenoData = pheno)
exprs_rma  <- affy::rma(raw)
cols <- data.frame(samples=colnames(exprs(exprs_rma)))
df <- left_join(cols,pheno %>% rownames_to_column('samples'))
df <- df %>%  separate(`knockdown:ch1`,into = c('pert_type','knockdown'),sep = '-')
df <- df %>% mutate(knockdown = ifelse(is.na(knockdown),'untreated',knockdown))
group<-as.factor(df$knockdown)
group <- relevel(group, ref = "untreated")
design<-model.matrix(~group)

# get the normalized gene expression and aggregate probes into genes
gex_rma <- exprs(exprs_rma)
anno <- AnnotationDbi::select(hgu133a.db,
                              keys = (rownames(gex_rma)),
                              columns = c("SYMBOL"),
                              keytype = "PROBEID")
anno <- subset(anno , !is.na(SYMBOL))
annog <- group_by(anno, PROBEID)
anno_summarized <-
  dplyr::summarize(annog, no_of_matches = n_distinct(SYMBOL))
anno_filtered <- filter(anno_summarized, no_of_matches == 1)
ids <- (rownames(gex_rma) %in% anno_filtered$PROBEID)
gex_rma <- subset(gex_rma,ids)
gex_rma <- as.data.frame(gex_rma) %>%
  rownames_to_column("PROBEID")
gex_rma <- left_join(gex_rma,anno)
gex_rma <- gex_rma  %>% dplyr::select(-PROBEID)
gex_rma <- aggregate(gex_rma[,1:(ncol(gex_rma)-1)],by=list(gex_rma$SYMBOL),FUN=median)
gex_rma <- gex_rma %>% column_to_rownames('Group.1')

# Calculate TF activity
dorotheaData = read.table('TF activities/annotation/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]

#Standarize gene expression
gex_standarized <- scale(t(gex_rma),scale = F)
hist(gex_standarized)

# Estimate TF activities
settings = list(verbose = TRUE, minsize = 5)
TF_activities = run_viper(t(gex_standarized), dorotheaData, options =  settings)

# Get nice data frame for visualization
TF_activities <- as.data.frame(TF_activities) %>% rownames_to_column('TF') %>% gather('samples','value',-TF)
TF_activities <- left_join(TF_activities,df %>% dplyr::select(samples,knockdown) %>% unique())

# plot results
results <- TF_activities %>% filter(TF=='FOXM1') %>% filter(knockdown %in% c('untreated','CONTROL','FOXM1','CDK2',
                                                                             'EGFR','CDK6','CDK1'))
ggdotplot(results,x='knockdown',y='value',fill='knockdown') +
  geom_hline(yintercept = 0,linetype='dashed',color='black',linewidth=1) +
  ylab('FOXM1 TF activity') + 
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle=90),
        legend.position='none')
ggsave('Model/CVL1000_Paper/ExperimentalValidation/geo_kos_foxm1.eps',
       device=cairo_ps,
       scale = 1,
       width = 6,
       height = 9,
       units = "in",
       dpi = 600)
