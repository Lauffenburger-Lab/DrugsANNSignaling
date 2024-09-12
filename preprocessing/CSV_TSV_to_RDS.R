library(tidyverse)
inputFile <- 'preprocessed_data/l1000_all_genes_lvl3_drugs_with_targets_exemplar.tsv'
outputFile <- 'preprocessed_data/l1000_all_genes_lvl3_drugs_with_targets_exemplar.rds'
type <- 'TSV'

### Load to file
if (type=='TSV'){
  df <- data.table::fread(inputFile,sep = "\t") %>% column_to_rownames('V1')
}else{
  df <- data.table::fread(inputFile) %>% column_to_rownames('V1')
}

### Load file
saveRDS(df,outputFile)