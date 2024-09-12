inputFile <- 'preprocessed_data/l1000_all_genes_lvl3_drugs_with_targets_exemplar.rds'
outputFile <- 'preprocessed_data/l1000_all_genes_lvl3_drugs_with_targets_exemplar.tsv'
type <- 'TSV'

### Load file
df <- readRDS(inputFile)
df <- as.data.frame(df)

### Save to file
if (type=='TSV'){
  data.table::fwrite(df,outputFile,sep = "\t",row.names = TRUE)
}else{
  data.table::fwrite(df,outputFile,row.names = TRUE)
}