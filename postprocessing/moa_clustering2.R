library(tidyverse)
library(ggsignif)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(viridis)
library(patchwork)
library(caret)
library(scales)
library(doFuture)
#options(future.globals.maxSize= 891289600)
#options(future.globals.maxSize= 4831838208)
# parallel set number of workers
cores <- 16
registerDoFuture()
plan(multisession,workers = cores)
library(doRNG)

df_drugbank <- readRDS('df_drugbank.rds')
drug_targets <- readRDS('drug_targets_space.rds')
#df_km <- readRDS('moa_cluster_elbow.rds')
#km <- readRDS('moa_clusters.rds')

### Load drug target data and build the MOA space---------
cmap_drugs <- readRDS('l1000_drugs_with_targets_all.rds')
annot <- read.delim('Network Construction/annotation/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab')
annot <- annot %>% select(Entry,c('Target'=Gene.names...primary..)) %>% unique()
cmap_drugs <-  cmap_drugs %>% mutate(Target=strsplit(Target,", ")) %>% unnest(Target)
cmap_drugs <- left_join(cmap_drugs,annot)
cmap_drugs <- cmap_drugs %>% filter(!is.na(Entry))
data <- readRDS("all_cmap_sigs_with_pert_info.rds")
data <- data %>% filter(is_exemplar==1)  %>%
  filter(pert_itime %in% c("6 h"))
data <- data %>% dplyr::select(sig_id,pert_iname,pert_type,pert_idose,pert_itime,cell_id,pert_id,canonical_smiles) %>% unique()
cmap_ctrls <- data %>% filter(pert_type=='ctl_vehicle' | pert_type=='ctl_untrt') %>% unique()
cmap_ctrls <- cmap_ctrls %>% mutate(canonical_smiles = ifelse(pert_id=='DMSO','CS(C)=O',
                                                              ifelse(pert_id=='H2O','O',
                                                                     ifelse(pert_id=='PBS','Ctrl_PBS','Ctrl_UnTrt'))))
cmap_ctrls <- cmap_ctrls %>% mutate(Entry=ifelse(pert_id=='DMSO',c('Q01344', 'Q8WXI7', 'P01106'),'None')) %>% unnest(Entry)
cmap_ctrls <- cmap_ctrls  %>%  mutate(Target=ifelse(Entry=='Q01344','IL5RA',ifelse(Entry=='P01106','MYC',ifelse(Entry=='Q8WXI7','MUC16',Entry)))) 
cmap_ctrls <- cmap_ctrls %>% mutate(MOA=NA,Phase=NA,Disease.Area='control',pubchem_cid=NA)
cmap_drugs <- left_join(cmap_drugs,data)
cmap_ctrls <- cmap_ctrls %>% select(colnames(cmap_drugs))
cmap_drugs <- rbind(cmap_drugs,cmap_ctrls)
cmap_drugs <- cmap_drugs %>% filter(!is.na(Entry))
cmap_drugs <- cmap_drugs %>% filter(Entry!='')
cmap_drugs <- cmap_drugs %>% filter(Entry!=' ')

cmap_drugs_unnested <-  cmap_drugs %>% mutate(Target=strsplit(Target,", ")) %>% unnest(Target)
pkn <- data.table::fread('Network Construction/model/pkn.tsv')
annot <- read.csv('TF activities/annotation/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t')
annot <- annot %>% select(Gene.names...primary..,Entry)
colnames(annot)[1] <- 'Target'
cmap_drugs_unnested <- left_join(cmap_drugs_unnested,annot)
cmap_drugs_unnested <- cmap_drugs_unnested %>% mutate(Entry=ifelse(Target=='None','None',Entry))
cmap_drugs_unnested <- cmap_drugs_unnested %>% filter(!is.na(Entry))
cmap_drugs_unnested <- cmap_drugs_unnested %>% select(canonical_smiles,Entry) %>% unique()
# Read training data (it has pbs and dmso)
long_drug_target <-  read.delim('Model/data/L1000_lvl3_allcells-drugs_targets.tsv')
colnames(long_drug_target)[1] <- 'canonical_smiles'
long_drug_target <- long_drug_target %>% gather('Entry','value',-canonical_smiles)
long_drug_target <- long_drug_target %>% filter(value==1) %>% select(-value) %>% unique()
long_drug_target <- rbind(long_drug_target,data.frame(canonical_smiles='Ctrl_PBS',Entry='None'))
cmap_drugs_unnested <- rbind(cmap_drugs_unnested,long_drug_target)
cmap_drugs_unnested <- cmap_drugs_unnested %>% unique()

## Manually add dmso from drugbank and also add more targets for other drugs too
library(dbparser)
library(XML)
read_drugbank_xml_db("full database.xml")
drugs <- drugs()
snp_effects <- drugs$snp_effects
drug_interactions <- drugs$interactions
snp_effects <- snp_effects %>% select(`protein-name`,`gene-symbol`,`uniprot-id`,`pubmed-id`,parent_key)
drug_interactions <- drug_interactions %>% select(-description)
drug_interactions <- left_join(drug_interactions,snp_effects)
df_drugbank <- drugs$general_information
df_drugbank <- left_join(df_drugbank,drug_interactions,by=c('primary_key'='drugbank-id','name'='name'))
df_drugbank$name <- tolower(df_drugbank$name)
df_drugbank <- df_drugbank %>%
  select(c('drugbank_id'='primary_key'),'name','cas_number',
         c('target_symbol'='gene-symbol'),c('target_uniprot'='uniprot-id')) %>% 
  filter(!is.na(target_symbol) | !is.na(target_uniprot)) %>% unique()
df_drugbank_cmap <- left_join(cmap_drugs %>% select(-Target) %>% mutate(cmap_name=tolower(pert_iname)) %>% unique(),
                              df_drugbank,
                              by=c('pert_iname'='name')) %>%
  select(-target_symbol,-cas_number) %>% filter(!is.na(target_uniprot)) %>% 
  filter(target_uniprot!='') %>% filter(target_uniprot!=' ')%>%
  unique()
df_drugbank_cmap <-  df_drugbank_cmap %>% select("pert_iname","canonical_smiles",
                                                 "Entry"='target_uniprot') %>% unique()
df_drugbank_cmap <-  df_drugbank_cmap %>% filter(!is.na(Entry)) %>% filter(Entry!='') %>% filter(Entry!=' ')
df_drugbank_cmap$Entry <- toupper(df_drugbank_cmap$Entry) 
df_drugbank_cmap <- df_drugbank_cmap %>% select(-pert_iname) %>% unique()
df_drugbank_not_in_cmap <- df_drugbank %>% filter(!(name %in% tolower(cmap_drugs$pert_iname))) %>%
  select(c('canonical_smiles'='name'),c('Entry'='target_uniprot')) %>% unique()

drug_targets_space <- rbind(cmap_drugs_unnested,df_drugbank_cmap,df_drugbank_not_in_cmap) %>% unique()


uniq_targets <- unique(drug_targets_space$Entry)

binary_classes <- function(sig,df,targets){
  #df <- df %>% select(sig_id,Target)
  #df <- df %>% filter(sig_id==sig)
  df <- df %>% select(canonical_smiles,Entry) %>% unique()
  df <- df %>% filter(canonical_smiles==sig)
  m <- as.matrix(df[,2])
  ind <- which(targets %in% m)
  n <- targets
  targets <- rep(0,length(targets))
  targets[ind] <- 1
  names(targets) <- n
  return(targets)
}
drug_targets <- t(sapply(unique(drug_targets_space$canonical_smiles),binary_classes,drug_targets_space,uniq_targets))
ind <- which(colnames(drug_targets)=='None')
drug_targets <- drug_targets[,-ind]
#drug_targets <- read.delim('Model/data/L1000_lvl3_allcells-drugs_targets.tsv')
#colnames(drug_targets)[1] <- 'drug'
#drug_targets <- drug_targets %>% column_to_rownames('drug')

### Network analysis based on target similarity-------------------------------
jaccard <- function(M=NULL,x=NULL,y=NULL){
  if (is.null(M)){
    intersection  <-  sum(x*y)
    union  <-  sum(x) + sum(y) - intersection
  }else{
    intersection <- M %*% t(M)
    union <- apply(as.matrix(rowSums(M)),1,'+',as.matrix(rowSums(M))) - intersection
  }
  J <- intersection/union
  return(J)
}
### Perform spectral clustering
# Create null distribution of observing a similarity score randomly
NullJaccard <- function(X,max_iter=100){
  #r <- nrow(X)
  #c <- ncol(X)
  Snull <- NULL
  for (i in 1:max_iter){
    #Y <- matrix(rbinom(r*c,1,0.5),r,c)
    Y <- Rfast::colShuffle(X)
    SIM <-  jaccard(Y)
    names(SIM) <- NULL
    SIM <- SIM[lower.tri(SIM,diag = F)]
    SIM <- SIM[which(!is.na(SIM))]
    Snull <- c(Snull,SIM)
    message(paste0('Finished iteration ',i,'/',max_iter))
  }
  
  return(Snull)
}

Snull <- NullJaccard(drug_targets,max_iter = 100)
#Snull <- Snull[which(!is.na(Snull))]
Snull <- readRDS('Snull_new_100iters.rds')
hist(Snull,main='Histogram of null distribution of jaccard similarities',xlab='jaccard similarity')
#saveRDS(Snull,'Snull_new_100iters.rds')
S <- jaccard(drug_targets)
ind <- grep('_',rownames(S))
S[ind[1],ind[1]] <- 1 # it gives NaN cause PBS has no targets at all but we can say that with itself it has similarity 1
S[ind[2],ind[2]] <- 1
hist(S)

## Create edge list for node2vec
name_map <- data.frame(name=rownames(S),id=seq(1,length(rownames(S))))
S[lower.tri(S,diag = T)] <- -10
df_sim <- gather(as.data.frame(S) %>% rownames_to_column('name1'),key='name2',value = 'tanimoto',-name1)
df_sim <- df_sim %>% filter(tanimoto!=-10)
df_sim <- df_sim %>% filter(tanimoto!=0)
df_sim <-  distinct(df_sim)
df_sim <- left_join(df_sim,name_map,by=c('name1'='name'))
df_sim <- left_join(df_sim,name_map,by=c('name2'='name'))
#df_sim <- distinct(df_sim %>% select(c('node_1'='id.x'),c('node_2'='id.y'),tanimoto))
data.table::fwrite(distinct(df_sim %>% select(c('node_1'='id.x'),c('node_2'='id.y'),tanimoto)),
                   'similarity_graph_all_weighted_newdata.edgelist',sep=' ',col.names = FALSE)

df_sim <- distinct(df_sim %>% filter(tanimoto>=0.477))
data.table::fwrite(distinct(df_sim %>% select(c('node_1'='id.x'),c('node_2'='id.y'),tanimoto)),
                   'similarity_graph_filtered_weighted_newdata.edgelist',sep=' ',col.names = FALSE)

data.table::fwrite(distinct(df_sim %>% select(c('node_1'='id.x'),c('node_2'='id.y'))),
                   'similarity_graph_filtered_newdata.edgelist',sep=' ',col.names = FALSE)

make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}

A <- make.affinity(S, 3)  # use 3 neighboors (includes self)
A[1:8,1:8]

#S[diag(S)] <- 0
#S[S>0.4] <- 1
#S[S<=0.4] <- 0
hist(S)
hist(A)
#S = 1-S # make it distance
D <- rowSums(A)
L <- diag(D)-A
L_norm <- Rfast::mat.mult(diag(D^(-1/2)),L) # random walk normalized if D^-1
L_norm <- Rfast::mat.mult(L_norm,diag(D^(-1/2))) # symmetric normalized

# First make a nice correlation matrix
#library("RColorBrewer")
#library(ComplexHeatmap)
#col <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
#od =  hclust(dist(S))$order
#m2 = S[od, od]
# Heatmap(S, rect_gp = gpar(type = "none"), column_dend_side = "bottom",
#         cell_fun = function(j, i, x, y, w, h, fill) {
#           if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
#             grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
#           }
#         })
#png('similarity_heatmap_target_space.png',width=12,height=8,units = "in",res=300)
#Heatmap(S,scale = "none",col=col)
#dev.off()
saveRDS(S,'Affinity.rds')
saveRDS(L,'L.rds')
saveRDS(L_norm,'L_norm.rds')

##plot graph
library(igraph)
#Adj <- A
Adj <- 1*(A>=0.5)
colnames(Adj) <- colnames(S)
rownames(Adj) <- rownames(S)
network <- graph_from_adjacency_matrix(Adj,mode='undirected', diag=F)
plot(network,vertex.label=NA,layout=layout.sphere,main="sphere",
     vertex.size=2)


ev = eigen(L_norm,symmetric = TRUE)
plot(1:80,rev(ev$values)[1:80])
k <-20
Z <- ev$vectors[,(ncol(ev$vectors)-k+1):ncol(ev$vectors)]
km_eigen <- kmeans(Z,k,iter.max = 20,nstart = 30)
hist(km_eigen[["cluster"]])
saveRDS(km_eigen,'km_eigen.rds')

df <- Z[,1:2]
df <- data.frame(df)
df$cluster <- factor(km_eigen$cluster,levels = seq(1,k,1))

p <- ggplot(df, aes(X1, X2))+
  geom_point(aes(col =cluster))  + labs(title="Eigen vector plot of drug-target interactions space") +
  xlab('X1') + ylab('X2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(p)

# t-SNE plot
library(umap)
umap <- umap(Z,n_components=2)
df_all <- data.frame(V1 = umap$layout[,1], V2 =umap$layout[,2])
df_all$cluster <- factor(km_eigen$cluster,levels = seq(1,k,1))


gumap <- ggplot(df_all, aes(V1, V2),alpha=0.2)+
  geom_point(aes(col =cluster)) + labs(title="u-MAP of drug-target interactions space") + 
  xlab('uMAP 1') + ylab('uMAP 2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(gumap)

## Clustering analysis of mechanism of action-------------

### Read drug-network and node2vec embs
S <- readRDS('S.rds')
name_map <- data.frame(name=rownames(S),id=seq(1,length(rownames(S))))
embs <- data.table::fread('../node2vec/emb/filtered_weighted_similarity_graph_60walks.emb',skip=1,sep = ' ')
colnames(embs)[1] <- 'id'
colnames(embs)[2:ncol(embs)] <- paste0('z',seq(0,ncol(embs)-2))

ks <- c(2,5,seq(10,400,10))
dunn <- NULL
ParallelKmeansCluster <- function(kcenters,embeddings){
  km <- kmeans(embeddings,kcenters,iter.max = 20, nstart = 30)
  dunn <- (km$tot.withinss)/(km$tot.withinss+km$betweenss)
  return(data.frame(k=kcenters,dunn=dunn))
}

library(doRNG)
df_km <- NULL
df_km <- foreach(k = ks) %dorng% {
  ParallelKmeansCluster(kcenters = k ,embeddings = embs[,2:ncol(embs)])
}
df_km <- do.call(rbind,df_km)
saveRDS(df_km,'moa_cluster_elbow_60walks.rds')
# for (i in 1:length(ks)){
#   km <- kmeans(as.matrix(drug_targets), ks[i],iter.max = 15, nstart = 30)
#   dunn[i] <- (km$tot.withinss)/(km$tot.withinss+km$betweenss)
# }
# df_km <- data.frame(dunn,k=ks)
## Elbow plot
library(ggplot2)
# library(KneeArrower)
#xy <- findCutoff(df_km$k,df_km$dunn,method='first')
x <- 20
y <-  0.32621732
png('kmeans_elbow_plot_node2vec_60walks.png',width=8,height=6,units = "in",res=600)
ggplot(df_km,aes(k,dunn)) +geom_line(size=1)+ geom_point() +ylim(c(0,1))+
  scale_x_continuous(breaks = c(2,seq(20,400,20))) + 
  geom_segment(aes(x=x,xend =x,y=0, yend=y),size=0.75,color="red",
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text",x=65,y=0.1,label="Elbow point",size=5)+
  xlab('k-number of clusters') + ylab('Dunn Index') +ggtitle('Elbow plot for optimal number of clusters')+
  theme(text = element_text(size=13),legend.position = "none",plot.title = element_text(hjust = 0.5))
dev.off()
k <- 20
km <- kmeans(embs[,2:ncol(embs)], k,iter.max = 20, nstart = 30)
#saveRDS(km,'moa_clusters_nod2vec_60walks.rds')

## PCA plot
library(factoextra)
pca.samples <- prcomp(embs[,2:ncol(embs)],scale = F)
fviz_eig(pca.samples,ncp=20)

df_pca <- pca.samples$x[,1:2]
df_pca <- data.frame(df_pca)
df_pca$cluster <- factor(km$cluster,levels = seq(1,k,1))

pca_plot <- ggplot(df_pca, aes(PC1, PC2))+
  geom_point(aes(col =cluster))  + labs(title="PCA plot of drug-target interactions space") +
  xlab('PC1') + ylab('PC2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(pca_plot)
png('pca_plot_drugtargets_space_node2vec_60walks.png',width=8,height=6,units = "in",res=600)
print(pca_plot)
dev.off()

# t-SNE plot
library(Rtsne)
perpl = DescTools::RoundTo(sqrt(nrow(embs)), multiple = 5, FUN = round)
#perpl= 50
init_dim = 8
iter = 1000
emb_size = ncol(embs[,2:ncol(embs)])
set.seed(42)
tsne_all <- Rtsne(embs[,2:ncol(embs)], 
                  dims = 2, perplexity=perpl, 
                  verbose=TRUE, 
                  max_iter = iter,
                  initial_dims = init_dim,
                  check_duplicates = T,
                  normalize = T,pca_scale = F,
                  num_threads = 15)
df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2])
df_all$cluster <- factor(km$cluster,levels = seq(1,k,1))
# df_all$cluster <- km$cluster


gtsne <- ggplot(df_all, aes(V1, V2),alpha=0.2)+
  geom_point(aes(col =cluster)) + labs(title="t-SNE plot of drug-target interactions space") + 
  xlab('Dim 1') + ylab('Dim 2')+
  xlab('t-SNE 1') + ylab('t-SNE 2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(gtsne)
png('tsne_plot_drugtargets_space_node2vec_60walks.png',width=8,height=6,units = "in",res=600)
print(gtsne)
dev.off()

# UMAP plot
library(umap)
set.seed(42)
map <- umap(embs[,2:ncol(embs)], n_components = 2,metric='euclidean',n_neighbors=10,
            init=pca.samples$x[,1:2],knn_repeats=10, random_state = 42)
df_map <- data.frame(V1 = map$layout[,1], V2 = map$layout[,2])
df_map$cluster <- factor(km$cluster,levels = seq(1,k,1))

gg_map <- ggplot(df_map, aes(V1, V2))+
  geom_point(aes(col =cluster))  + labs(title="UMAP plot of drug-target interactions space") +
  xlab('uMAP-1') + ylab('uMAP-2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(gg_map)

png('umap_plot_drugtargets_space_node2vec_60walks.png',width=8,height=6,units = "in",res=600)
print(gtsne)
dev.off()


### Assing super-MOA to each drug
embs <- left_join(embs,name_map) %>% select(-id) %>% column_to_rownames('name')
drug_moa <- data.frame(canonical_smiles=rownames(embs),moa=km$cluster)
drug_targets_moa <- left_join(drug_targets_space,drug_moa)
moa_targets_summary <- drug_targets_moa %>% group_by(canonical_smiles) %>% mutate(n_targets=n_distinct(Entry)) %>% ungroup() %>%
  group_by(moa) %>% mutate(median_targetrs_per_moa=median(n_targets)) %>% ungroup()
moa_targets_summary <-  moa_targets_summary %>% select(-Entry,-n_targets) %>% unique()
moa_targets_summary <- moa_targets_summary %>% group_by(moa) %>% mutate(n_drugs=n_distinct(canonical_smiles)) %>% ungroup()  
moa_targets_summary$moa <- factor(moa_targets_summary$moa,levels = seq(1,k,1)) 
cluster_barplot <- ggplot(moa_targets_summary %>% filter(!is.na(moa)),aes(x=moa,y=median_targetrs_per_moa,fill=n_drugs)) +
  geom_col(position="identity") +scale_fill_gradient(low="blue", high="red") + xlab('Cluster') + 
  ylab('Median number of targets per drug') + labs(fill='Cluster size') + 
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(cluster_barplot)
png('cluster_barplot_node2vec_60walks.png',width=8,height=6,units = "in",res=600)
print(cluster_barplot)
dev.off()

cmap_drugs_moa <- left_join(cmap_drugs_unnested,drug_moa) %>% select(-Entry) %>% unique()
# Now keep only those in in my model training space
drug_targets_model <- read.delim('Model/data/L1000_lvl3_allcells-drugs_targets.tsv')
colnames(drug_targets_model)[1] <- 'drug'
drug_targets_model <- drug_targets_model %>% column_to_rownames('drug')
cmap_drugs_moa_inmodel <- cmap_drugs_moa %>% filter(canonical_smiles %in% rownames(drug_targets_model))
drugs_in_model <- cmap_drugs_moa_inmodel$canonical_smiles[!is.na(cmap_drugs_moa_inmodel$moa)]

# library(Rtsne)
# #perpl = DescTools::RoundTo(sqrt(nrow(drug_targets)), multiple = 5, FUN = round)
# perpl= 50
# init_dim = 5
# iter = 1000
# emb_size = ncol(embs)
# set.seed(42)
# tsne_all <- Rtsne(embs, 
#                   dims = 2, perplexity=perpl, 
#                   verbose=TRUE, 
#                   max_iter = iter,
#                   initial_dims = init_dim,
#                   check_duplicates = T,
#                   normalize = T,pca_scale = F,
#                   num_threads = 15)
# df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2])
# df_all$cluster <- factor(km$cluster,levels = seq(1,k,1))
# df_all <- df_all[which(rownames(drug_targets) %in% cmap_drugs_moa_inmodel$canonical_smiles),] 
# df_all$cluster <- factor(df_all$cluster,levels = unique(df_all$cluster)[order(unique(df_all$cluster))])
df_all <- df_all[which(rownames(embs) %in% drugs_in_model),]

gtsne <- ggplot(df_all, aes(V1, V2),alpha=0.2)+
  geom_point(aes(col =cluster)) + labs(title="t-SNE plot of drug-target interactions in our model") + 
  xlab('Dim 1') + ylab('Dim 2')+
  xlab('t-SNE 1') + ylab('t-SNE 2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(gtsne)
png('tsne_plot_drugtargets_in_model_node2vec_60walks.png',width=8,height=6,units = "in",res=600)
print(gtsne)
dev.off()

pca.samples <- prcomp(embs[,2:ncol(embs)],scale = F)
df_pca <- pca.samples$x[,1:2]
df_pca <- data.frame(df_pca)
df_pca$cluster <- factor(km$cluster,levels = seq(1,k,1))
df_pca <- df_pca[which(rownames(embs) %in% drugs_in_model),]
df_pca$cluster <- factor(df_pca$cluster,levels = unique(df_pca$cluster)[order(unique(df_pca$cluster))])
pca_plot <- ggplot(df_pca, aes(PC1, PC2))+
  geom_point(aes(col =cluster))  + labs(title="PCA plot of drug-target interactions space") +
  xlab('PC1') + ylab('PC2')+
  theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
print(pca_plot)

png('pca_plot_drugtargets_in_model_node2vec_6walks.png',width=8,height=6,units = "in",res=600)
print(pca_plot)
dev.off()


# Those disconnected from the graph S assign them their self as MoA
cmap_drugs_moa_inmodel <- cmap_drugs_moa_inmodel %>% mutate(moa=ifelse(is.na(moa),canonical_smiles,moa))
# write.csv(cmap_drugs_moa_inmodel,'models_drugs_moa_node2vec.csv')

## Similarities scatter plot----
ecfp4_sims <- data.table::fread('all_drugbank_plus_cmap_ecfp4_similarities.csv',header = T) %>% column_to_rownames('V1')
ecfp4_sims <- as.matrix(ecfp4_sims)
S <- readRDS('S.rds')
cols <- which(colnames(S) %in% colnames(ecfp4_sims))
rows <- which(rownames(S) %in% rownames(ecfp4_sims))
S <- S[rows,cols]

cols2 <- which(colnames(ecfp4_sims) %in% colnames(S))
rows2 <- which(rownames(ecfp4_sims) %in% rownames(S))
ecfp4_sims <- ecfp4_sims[rows2,cols2]

ecfp4_sims[lower.tri(ecfp4_sims,diag = T)] <- 10
ecfp4_sims <- reshape2::melt(ecfp4_sims)
ecfp4_sims <- ecfp4_sims %>% filter(value != 10)
ecfp4_sims <- ecfp4_sims %>% filter(!is.na(value))
  
S[lower.tri(S,diag = T)] <- 10
S <- reshape2::melt(S)
S <- S %>% filter(value != 10)

df <- left_join(S,ecfp4_sims,by=c("Var1","Var2"))
df <- distinct(df %>% filter(!is.na(value.x))  %>% filter(!is.na(value.y)))

library(ggplot2)
library(ggpubr)
p <- ggplot(df,aes(x=value.x,y=value.y)) +geom_point(size=1) +
  xlab('Target tanimoto similarity') + ylab('ECFP4 chemical similarity') + 
  ggtitle('Correlation between chemical similarity and target similarity of drugs')+
  stat_cor(method = "pearson")+
  theme(text = element_text(size=13))
png('ecfp4_vs_target_smi.png',width=10,height = 10,units = "in",res=600)
print(p)
dev.off()

### Load inferred tagets and evaluate similarity in MoA---------------
drug_targets <- readRDS('drug_targets_space.rds')
no_models <- 51
threshold <- 39/no_models
merged_interactions <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
drug_targets_prior <- as.data.frame(drug_targets) %>% rownames_to_column('drug') %>% 
  gather('variable','Prior knowledge',-drug)
merged_interactions <- merged_interactions %>% select(-`Prior knowledge`)
merged_interactions <- left_join(merged_interactions,drug_targets_prior)
merged_interactions <- merged_interactions %>% filter(!is.na(`Prior knowledge`))
merged_interactions <- merged_interactions %>% mutate(Inferred = 1*(Inferred=='Interaction'))
merged_interactions <- merged_interactions %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions <- merged_interactions %>% group_by(drug,variable) %>% 
  mutate(model_counts = sum(Inferred==1)) %>% ungroup()
merged_interactions <- merged_interactions %>% mutate(model_counts=model_counts/no_models)
merged_interactions <- distinct(merged_interactions %>% mutate(Inferred = ifelse(model_counts>=threshold,1,0)) %>% 
           select("drug","variable","Prior knowledge","Inferred"))
### APO EDO KAI META PREPEI NA TO DO ###
embs_original <- data.table::fread('../node2vec/emb/filtered_weighted_similarity_graph_60walks.emb',skip=1,sep = ' ')
colnames(embs_original)[1] <- 'id'
colnames(embs_original)[2:ncol(embs_original)] <- paste0('z',seq(0,ncol(embs_original)-2))
S_original <- readRDS('S.rds')
#embs_original$cluster <- km_original$cluster 
name_map_original <- data.frame(name=rownames(S_original),id=seq(1,length(rownames(S_original))))
embs_original <- left_join(embs_original,name_map_original)
embs_original <- embs_original %>% dplyr::select(-id)

### Find new embeddings and try to then perform clustering and see if a drug is in the same cluster both times.
drug_targets_original <- readRDS('drug_targets_space.rds')
rows <- which(rownames(drug_targets_original) %in% merged_interactions$drug)
cols <- which(colnames(drug_targets_original) %in% merged_interactions$variable)
drug_targets_new <- drug_targets_original
drug_targets_new[rows,cols] <- 0

new_interactions <- merged_interactions %>% select(drug,variable,Inferred) %>% unique()
new_interactions <- new_interactions %>% spread('variable','Inferred')
new_interactions <- as.matrix(new_interactions %>% column_to_rownames('drug'))
#new_interactions <- 1*(new_interactions=='Interaction')
drug_targets_new[rownames(new_interactions),colnames(new_interactions)] <- new_interactions

### Check jaccard similarity of new target drugs and old
### Network analysis based on target similarity
jaccard <- function(M=NULL,x=NULL,y=NULL){
  if (is.null(M)){
    intersection  <-  sum(x*y)
    union  <-  sum(x) + sum(y) - intersection
  }else{
    intersection <- M %*% t(M)
    union <- apply(as.matrix(rowSums(M)),1,'+',as.matrix(rowSums(M))) - intersection
  }
  J <- intersection/union
  return(J)
}
merged_drug_targets <- rbind(drug_targets_original[rownames(new_interactions),colnames(new_interactions)],
                             drug_targets_new[rownames(new_interactions),colnames(new_interactions)])
id <- paste0('v',seq(1:nrow(merged_drug_targets)))
all_rows_map <- c(rownames(new_interactions),rownames(new_interactions))
all_rows_map <- data.frame(drug_name=all_rows_map,id=id)
rownames(merged_drug_targets) <- id  
S_merged <- jaccard(merged_drug_targets)
self_sim_merged <- S_merged[paste0('v',seq(1,length(unique(merged_interactions$drug)))),
                            paste0('v',seq(length(unique(merged_interactions$drug))+1,2*length(unique(merged_interactions$drug))))]
self_sim_merged <- diag(self_sim_merged)
png('self_drug_target_smililarity_after_inference_withsometargets_onlycmap.png',width = 9,height = 6,units = 'in',res = 600)
hist(self_sim_merged[which(rowSums(drug_targets_new[rownames(new_interactions),])!=0)],breaks = 20,
     main='Target similarity of the same drug if some target is inferred',
     xlab = 'Jaccard similarity',
     ylab= 'Counts')
abline(v=0.477,col="red",lwd=2,lty='dashed')
text(0.34, 12, 'Similarity threshold')
dev.off()
png('self_drug_target_smililarity_after_inference_withzeros_onlycmap.png',width = 9,height = 6,units = 'in',res = 600)
hist(self_sim_merged,breaks = 20,
     main='Target similarity of the same drug',
     xlab = 'Jaccard similarity',
     ylab= 'Counts')
abline(v=0.477,col="red",lwd=2,lty='dashed')
text(0.34, 80, 'Similarity threshold')
text(0.34, 40, paste0(100*round(sum(1*(self_sim_merged<0.477))/length(self_sim_merged),4),"%"))
text(0.614, 40, paste0(100*round(sum(1*(self_sim_merged>=0.477))/length(self_sim_merged),4),"%"))
dev.off()

S_merged[which((merged_drug_targets %*% t(merged_drug_targets))==0)] <- 0
pbs_ids <- all_rows_map$id[grep('Ctrl',all_rows_map$drug_name)]
print(pbs_ids)
if (length(pbs_ids)>0){
  S_merged[pbs_ids[1],pbs_ids[1]] <- 1
  S_merged[pbs_ids[2],pbs_ids[2]] <- 1
}
S_merged[lower.tri(S_merged,diag = T)] <- -10
df_sim_merged <- gather(as.data.frame(S_merged) %>% rownames_to_column('name1'),key='name2',value = 'tanimoto',-name1)
df_sim_merged <- df_sim_merged %>% filter(tanimoto!=-10)
df_sim_merged <- df_sim_merged %>% 
  filter((name1 %in% paste0('v',seq(length(unique(merged_interactions$drug))+1,2*length(unique(merged_interactions$drug)))))&
           (name2 %in% paste0('v',seq(length(unique(merged_interactions$drug))+1,2*length(unique(merged_interactions$drug))))))
df_sim_merged <-  distinct(df_sim_merged)
#df_sim_merged <- df_sim_merged %>% filter(tanimoto!=0)
df_sim_merged <- left_join(df_sim_merged,all_rows_map,by=c('name1'='id'))
df_sim_merged <- left_join(df_sim_merged,all_rows_map,by=c('name2'='id'))
df_sim_merged <- df_sim_merged %>% select(drug_name.x,drug_name.y,tanimoto) %>% unique()
df_sim_merged <- df_sim_merged %>% mutate(interaction = ifelse(tanimoto>=0.477,'Similar','Not similar'))
#df_sim_merged <- distinct(df_sim_merged %>% filter(tanimoto>=0.4))
S_true <- S_original[rownames(new_interactions),rownames(new_interactions)] 
#S_true <- 1*(S_true>=0.4)
df_sim_true <- gather(as.data.frame(S_true) %>% rownames_to_column('drug_name.x'),key='drug_name.y',value = 'true_tanimoto',-drug_name.x)
df_sim_true <- df_sim_true %>% mutate(true_interaction=ifelse(true_tanimoto>=0.477,'Similar','Not similar'))
df_sim_merged <- left_join(df_sim_merged,df_sim_true)

p <- ggscatter(df_sim_merged,x='true_tanimoto',y='tanimoto',
               rug = TRUE,
               alpha = 0.5,size=1,
               cor.coef=T,cor.coef.size = 5,cor.coef.coord = c(0.05, 0.9)) + 
  geom_abline(slope=1,intercept=0,color='red',lty=2,size=1)+
  #stat_cor(method="pearson", label.x=0.1, label.y=0.75,size=5)+
  ggtitle('Inferred target similarity') + xlab('reference jaccard similarity') + ylab('predicted jaccard similarity') +
  xlim(c(0,1))+ylim(c(0,1))+
  theme_minimal(base_family = "serif",base_size = 15)+
  theme(plot.title = element_text(hjust = 0.5,size=15))
print(p)
png('jaccard_correlation_60walks_onlycmap.png',width=9,height=9,units = "in",res = 600)
print(p)
#ggMarginal(p, type = "boxplot",)
dev.off()
#Creating confusion matrix
preds <- as.factor(df_sim_merged$interaction)
true <-   as.factor(df_sim_merged$true_interaction)
confusion <- confusionMatrix(data=factor(preds, levels=rev(levels(preds))), 
                             reference = factor(true, levels=rev(levels(true))))
confusion
plt <- as.data.frame(confusion$table)
#plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
plt$Reference <- factor(plt$Reference, levels=rev(levels(plt$Reference)))
ggplot(plt, aes(Prediction,Reference, fill= Freq)) + 
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient2() +
  labs(x = "Reference",y = "Prediction") + ggtitle('Confusion matrix for similarity in the drug-target space')+
  theme_minimal(base_family = "serif",base_size = 15)
png('moa_performance_60walks_onlycmap.png',width=9,height=9,units = "in",res = 600)
ggplot(plt, aes(Prediction,Reference, fill= Freq)) + 
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient2() +
  labs(x = "Reference",y = "Prediction") + ggtitle('Confusion matrix for similarity in the drug-target space')+
  theme_minimal(base_family = "serif",base_size = 15)
dev.off()

Snew <- jaccard(drug_targets_new)
#ind <- which(rownames(Snew)=='Ctrl_PBS')
ind <- grep('Ctrl',rownames(Snew))
Snew[ind,ind] <- 1 # it gives NaN cause PBS has no targets at all but we can say that with itself it has similarity 1
saveRDS(Snew,'Snew_60walks.rds')
## Creta edge list for node2vec
name_map <- data.frame(name=rownames(Snew),id=seq(1,length(rownames(Snew))))
Snew[lower.tri(Snew,diag = T)] <- -10
df_sim <- gather(as.data.frame(Snew) %>% rownames_to_column('name1'),key='name2',value = 'tanimoto',-name1)
df_sim <- df_sim %>% filter(tanimoto!=-10)
df_sim <- df_sim %>% filter(tanimoto!=0)
df_sim <-  distinct(df_sim)
df_sim <- left_join(df_sim,name_map,by=c('name1'='name'))
df_sim <- left_join(df_sim,name_map,by=c('name2'='name'))
df_sim <- distinct(df_sim %>% filter(tanimoto>=0.4))
data.table::fwrite(distinct(df_sim %>% select(c('node_1'='id.x'),c('node_2'='id.y'),tanimoto)),
                   'similarity_graph_filtered_weighted_newinferred_60walks.edgelist',sep=' ',col.names = FALSE)

# Snew <- readRDS('Snew.rds')
# name_map <- data.frame(name=rownames(Snew),id=seq(1,length(rownames(Snew))))
# embs <- data.table::fread('../node2vec/emb/filtered_weighted_similarity_graph_newinferred.emb',skip=1,sep = ' ')
# colnames(embs)[1] <- 'id'
# colnames(embs)[2:ncol(embs)] <- paste0('z',seq(0,ncol(embs)-2))
# embs <- left_join(embs,name_map) %>% select(-id)
# 
# embs_merged <- rbind(embs,embs_original)
# dim <- 128
# k <- 20
# km <- kmeans(embs_merged[,1:dim], k,iter.max = 15, nstart = 30)
# saveRDS(km,'moa_clusters_nod2vec_merged_withnew.rds')
# 
# ## PCA plot
# library(factoextra)
# pca.samples <- prcomp(embs_merged[,1:dim],scale = F)
# fviz_eig(pca.samples,ncp=20)
# 
# df_pca <- pca.samples$x[,1:2]
# df_pca <- data.frame(df_pca)
# df_pca$cluster <- factor(km$cluster,levels = seq(1,k,1))
# 
# pca_plot <- ggplot(df_pca, aes(PC1, PC2))+
#   geom_point(aes(col =cluster))  + labs(title="PCA plot of drug-target interactions space") +
#   xlab('PC1') + ylab('PC2')+
#   theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
# print(pca_plot)
# 
# # t-SNE plot
# library(Rtsne)
# perpl = DescTools::RoundTo(sqrt(nrow(embs_merged)), multiple = 5, FUN = round)
# #perpl= 50
# init_dim = 5
# iter = 1000
# emb_size = dim
# set.seed(42)
# tsne_all <- Rtsne(embs_merged[,1:dim], 
#                   dims = 2, perplexity=perpl, 
#                   verbose=TRUE, 
#                   max_iter = iter,
#                   initial_dims = init_dim,
#                   check_duplicates = T,
#                   normalize = T,pca_scale = F,
#                   num_threads = 15)
# df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2])
# df_all$cluster <- factor(km$cluster,levels = seq(1,k,1))
# 
# 
# gtsne <- ggplot(df_all, aes(V1, V2),alpha=0.2)+
#   geom_point(aes(col =cluster)) + labs(title="t-SNE plot of drug-target interactions space") + 
#   xlab('Dim 1') + ylab('Dim 2')+
#   xlab('t-SNE 1') + ylab('t-SNE 2')+
#   theme(text = element_text(size=13),plot.title = element_text(hjust = 0.5))
# print(gtsne)
# png('tsne_plot_drugtargets_space_node2vec_withnew_infered.png',width=8,height=6,units = "in",res=600)
# print(gtsne)
# dev.off()
# 
# df_cluster <- left_join(embs_merged[1:nrow(embs),] %>% select(name) %>% 
#                           mutate(cluster=km$cluster[1:nrow(embs)]),
#                         embs_merged[nrow(embs)+1:nrow(embs_merged),] %>% select(name) %>% 
#                           mutate(cluster=km$cluster[nrow(embs)+1:nrow(embs_merged)]),
#                         by='name')
# df_cluster <- df_cluster %>% filter(name %in% unique(merged_interactions$drug))
# df_cluster <- df_cluster %>% filter(!is.na(cluster.x) & !is.na(cluster.y))
# 
# embs_merged <- embs_merged %>% mutate(id=paste0('v',seq(1,nrow(embs_merged))))
# ## Cosine similarity
# library(lsa)
# mat <- t(embs_merged[,1:dim])
# sim <- cosine(mat)
# colnames(sim) <- embs_merged$id
# rownames(sim) <- embs_merged$id
# # Conver to long format data frame
# # Keep only unique (non-self) pairs
# sim[lower.tri(sim,diag = T)] <- NA
# sim <- reshape2::melt(sim)
# sim <- sim %>% filter(!is.na(value))
# # Merge meta-data info and distances values
# sim <- left_join(sim,embs_merged %>% select(id,name),by = c("Var1"="id"))
# sim <- left_join(sim,embs_merged %>% select(id,name),by = c("Var2"="id"))
# sim <- sim %>% filter(!is.na(value))
# sim <- sim %>% mutate(pair = ifelse(name.x==name.y,'Same drugs','Random pair'))
# sim <- sim %>% mutate(metric = 'cosine')
# 
# ## Eucledian distance
# dist <- as.matrix(dist(embs_merged[,1:dim], method = 'euclidean'))
# colnames(dist) <- embs_merged$id
# rownames(dist) <- embs_merged$id
# # Conver to long format data frame
# # Keep only unique (non-self) pairs
# dist[lower.tri(dist,diag = T)] <- NA
# dist <- reshape2::melt(dist)
# dist <- dist %>% filter(!is.na(value))
# # Merge meta-data info and distances values
# dist <- left_join(dist,embs_merged %>% select(id,name),by = c("Var1"="id"))
# dist <- left_join(dist,embs_merged %>% select(id,name),by = c("Var2"="id"))
# dist <- dist %>% filter(!is.na(value))
# dist <- dist %>% mutate(pair = ifelse(name.x==name.y,'Same drugs','Random pair'))
# dist <- dist %>% mutate(metric = 'euclidean')
# 
# ### Visualize
# # sim$value <- 1 - sim$value#
# #rbind(dist,sim)
# violin_separation <- ggplot(sim, aes(x=pair, y=value, fill = pair)) + 
#   geom_violin(position = position_dodge(width = 1),width = 1)+geom_boxplot(position = position_dodge(width = 1),width = 0.05,
#                                                                            outlier.shape = NA)+
#   scale_fill_discrete(name="Embeddings distance distribution",
#                       labels=c("Random pair","Same drugs"))+
#   ylim(-1,1)+ #ylim(0,max(dist$value))
#   xlab("")+ylab("Cosine Similarity")+ 
#   theme(axis.ticks.x=element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(),text = element_text(family = "serif",size = 17),legend.position = "bottom")+
#   theme_minimal(base_family = "serif",base_size = 17)#+facet_wrap(~ metric)
# violin_separation <- violin_separation + theme(legend.position = "bottom")
# print(violin_separation)

#### Compare inferred with drugbank------------------------------
library(tidyverse)
library(caret)
library(reshape2)
# ensemble_table <- readRDS('Model/CVL1000_Paper/A375_ensembles/ensemble_table.rds')
# new_predictions <- readRDS('Model/CVL1000_Paper/A375_ensembles/new_predictions.rds')
# pvals <- readRDS('Model/CVL1000_Paper/A375_ensembles/pvals.rds')
# drug_bank_tp <- readRDS('Model/CVL1000_Paper/A375_ensembles/drug_bank_tp.rds')
# ensemble_table_null <- readRDS('Model/CVL1000_Paper/A375_ensembles/ensemble_table_null.rds')
no_models <- 50
ensemble_table <- NULL
drug_targets <- readRDS('drug_targets_space.rds')
merged_interactions_all <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
i <- 1
# merged_interactions <- data.table::fread(paste0('Model/CVL1000_Paper/FinalEnsemble/InteractionScores/l1000_modeltype4_lamda6_VCAP_mergedInteractions_ROC_',
#                                                 i-1,
#                                                 '.csv'),header=T) %>%
#   column_to_rownames('V1')
merged_interactions <- merged_interactions_all %>% filter(model_no==i-1) %>% dplyr::select(-model_no)
# drugs_of_interest <- c('C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13',
#                        'CCc1nn(CCCN2CCN(CC2)c2cccc(Cl)c2)c(=O)n1CCOc1ccccc1',
#                        'COc1ccc2n(C(=O)c3cccc(Cl)c3Cl)c(C)c(CCN3CCOCC3)c2c1',
#                        'COC(=O)C1=C(C)NC(C)=C(C1C)C(=O)OCCSc1ccccc1',
#                        '[O-][N+](=O)c1ccc(Cl)c(c1)C(=O)Nc1ccncc1')
# merged_interactions <- merged_interactions %>% filter(drug %in% drugs_of_interest)
drug_targets <- drug_targets[which(rownames(drug_targets) %in% merged_interactions$drug),which(colnames(drug_targets) %in% merged_interactions$variable)]
new_predictions <- NULL
ensemble_table_null <- NULL
drug_bank_tp <- NULL 
pvals <- NULL
for (i in 1:no_models){
  # merged_interactions <- data.table::fread(paste0('Model/CVL1000_Paper/FinalEnsemble/InteractionScores/l1000_modeltype4_lamda6_VCAP_mergedInteractions_ROC_',
  #                                                 i-1,
  #                                                 '.csv'),header=T) %>%
  #   column_to_rownames('V1')
  merged_interactions <- merged_interactions_all %>% filter(model_no == i-1) %>% select(-model_no)
  # merged_interactions <- merged_interactions %>% filter(drug %in% drugs_of_interest)
  drug_targets_prior <- as.data.frame(drug_targets) %>% rownames_to_column('drug') %>% 
    gather('variable','Prior knowledge',-drug)
  # training_prior <- merged_interactions %>% filter(`Prior knowledge`=='Interaction')
  # training_prior <- training_prior %>% select(-`Prior knowledge`,-Inferred) %>% mutate(mode='training_prior')
  # drug_targets_prior <- anti_join(drug_targets_prior,training_prior)
  potential_new_pairs <- merged_interactions %>% filter(`Prior knowledge`=='No interaction') %>% 
    select(-`Prior knowledge`) %>%
    mutate(Inferred='No interaction')
  merged_interactions <- merged_interactions %>% filter(!(`Prior knowledge`==Inferred & `Prior knowledge`=='Interaction'))
  merged_interactions <- merged_interactions %>% select(-`Prior knowledge`)
  merged_interactions <- suppressMessages(left_join(merged_interactions,drug_targets_prior))
  merged_interactions <- merged_interactions %>% filter(!is.na(`Prior knowledge`))
  merged_interactions <- merged_interactions %>%
    mutate(`Prior knowledge`=ifelse(`Prior knowledge`==0,'No interaction','Interaction'))
  confusion <- confusionMatrix(data=as.factor(merged_interactions$Inferred), 
                               reference = as.factor(merged_interactions$`Prior knowledge`))
  ensemble_table[[i]] <-  confusion$table
  new_predictions[i] <- sum(merged_interactions$Inferred == 'Interaction')
  
  # sample a number of predictions 10000 times
  table_null <- NULL
  drug_bank_tp_null <- NULL
  for (j in 1:1000){
    tmp <- potential_new_pairs
    inds <- sample(seq(1:nrow(tmp)),new_predictions[i])
    tmp$Inferred[inds] <- 'Interaction'
    tmp <- suppressMessages(left_join(tmp,drug_targets_prior))
    tmp <- tmp %>% filter(!is.na(`Prior knowledge`))
    tmp <- tmp %>%
      mutate(`Prior knowledge`=ifelse(`Prior knowledge`==0,'No interaction','Interaction'))
    confusion_null <- confusionMatrix(data=as.factor(tmp$Inferred), 
                                 reference = as.factor(tmp$`Prior knowledge`))
    table_null[[j]] <- confusion_null$table
    drug_bank_tp_null[j] <- confusion_null$table['Interaction','Interaction']
    if (j %% 250 == 0){
      print(paste0('Finished null iteration ',j, ' of model ',i))
    }
  }
  ensemble_mat <- do.call(cbind,table_null)
  ensemble_mat <- array(ensemble_mat,c(dim=dim(table_null[[1]]),length(table_null)))
  #ensemble_mat_mean <- apply(ensemble_mat, c(1,2), mean, na.rm = TRUE)
  ensemble_table_null[[i]] <- apply(ensemble_mat, c(1,2), mean, na.rm = TRUE)
  drug_bank_tp[[i]] <- drug_bank_tp_null
  pvals[i] <- sum(drug_bank_tp_null>=confusion$table['Interaction','Interaction'])/length(drug_bank_tp_null)
  message(paste0('Finished model ',i))
}
# saveRDS(ensemble_table,'Model/CVL1000_Paper/A375_ensembles/ensemble_table_interest.rds')
# saveRDS(new_predictions,'Model/CVL1000_Paper/A375_ensembles/new_predictions_interest.rds')
# saveRDS(pvals,'Model/CVL1000_Paper/A375_ensembles/pvals_interest.rds')
# saveRDS(drug_bank_tp,'Model/CVL1000_Paper/A375_ensembles/drug_bank_tp_interest.rds')
# saveRDS(ensemble_table_null,'Model/CVL1000_Paper/A375_ensembles/ensemble_table_null_interest.rds')
ensemble_table <- readRDS('Model/CVL1000_Paper/A375_ensembles/ensemble_table_interest.rds')
new_predictions <- readRDS('Model/CVL1000_Paper/A375_ensembles/new_predictions_interest.rds')
pvals <- readRDS('Model/CVL1000_Paper/A375_ensembles/pvals_interest.rds')
drug_bank_tp <- readRDS('Model/CVL1000_Paper/A375_ensembles/drug_bank_tp_interest.rds')
ensemble_table_null <- readRDS('Model/CVL1000_Paper/A375_ensembles/ensemble_table_null_interest.rds')

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
plot_confusion <- ggplot(ensemble_final, 
                         aes(x = Var1, y = Var2, fill = Count)) +
  geom_tile() +
  scale_fill_viridis(trans='log10', limits = c(1, max(ensemble_final$Count)),breaks=c(1,10,100,1000,10000)) +
  labs(fill='Count')+
  geom_text(aes(label = paste(round(Count, 3), round(Count_sd/sqrt(10), 2), sep = "\u00B1")), size = 8) +
  # annotate('text',x = 1.5,y=1.5,
  #          label=paste0(paste0('F1 score = ',paste(round(mean(f1_mouse)*100,2),round(sd(f1_mouse)*100/sqrt(10),2),sep="\u00B1"),'%'),
  #                       "\n",
  #                       paste0('Accuracy = ',paste(round(mean(acc_mouse)*100,2),round(sd(acc_mouse)*100/sqrt(10),2),sep="\u00B1"),'%')),
  #          size=8,
  #          fontface =2)+
  labs(title = "Confusion matrix of off-target and drugbank prior not in cmap space",
       x = "Inferred",
       y = "Prior knowledge") +
  theme_minimal() +
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24, face = "bold"),
        legend.text = element_text(family = 'Arial',size = 20),
        plot.title = element_text(family = 'Arial',size = 20, face = "bold", hjust = 0.5))
print(plot_confusion)

## Null distribution characteristics
# hist(new_predictions)
ensemble_mat_null <- do.call(cbind,ensemble_table_null)
ensemble_mat_null <- array(ensemble_mat_null,c(dim=dim(ensemble_table_null[[1]]),length(ensemble_table_null)))
ensemble_mat_mean_null <- apply(ensemble_mat_null, c(1,2), mean, na.rm = TRUE)
colnames(ensemble_mat_mean_null) <- colnames(ensemble_table_null[[1]])
rownames(ensemble_mat_mean_null) <- rownames(ensemble_table_null[[1]])
ensemble_mat_sd_null <- apply(ensemble_mat_null, c(1,2), sd, na.rm = TRUE)
colnames(ensemble_mat_sd_null) <- colnames(ensemble_table_null[[1]])
rownames(ensemble_mat_sd_null) <- rownames(ensemble_table_null[[1]])
ensemble_mat_mean_null <- melt(ensemble_mat_mean_null, variable.name = c("Reference", "Prediction"), value.name = "Count")
ensemble_mat_sd_null <- melt(ensemble_mat_sd_null, variable.name = c("Reference", "Prediction"), value.name = "Count_sd")
ensemble_final_null <- left_join(ensemble_mat_mean_null,ensemble_mat_sd_null)
ensemble_final_null <- ensemble_final_null %>% mutate(Var1=ifelse(Var1==1,'Interaction','No interaction')) %>%
  mutate(Var2=ifelse(Var2==1,'Interaction','No interaction'))
plot_confusion_null <- ggplot(ensemble_final_null, 
                         aes(x = Var2, y = Var1, fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, max(ensemble_final_null$Count))) +
  labs(fill='Count')+
  geom_text(aes(label = paste(round(Count, 3), round(Count_sd/sqrt(10), 2), sep = "\u00B1")), size = 6) +
  labs(title = "Confusion matrix of random off-target prediction",
       x = "True Class",
       y = "Predicted Class") +
  theme_minimal() +
  theme(text = element_text(family = 'Arial',size=16),
        axis.text = element_text(family = 'Arial',size = 16),
        axis.title = element_text(family = 'Arial',size = 16, face = "bold"),
        legend.text = element_text(family = 'Arial',size = 16),
        plot.title = element_text(family = 'Arial',size = 16, face = "bold", hjust = 0.5))
print(plot_confusion_null)
# merged_interactions <- data.table::fread('Model/CVL1000_Paper/FinalEnsemble/InteractionScores/l1000_modeltype4_lamda6_VCAP_mergedInteractions_ROC_0.csv',header=T) %>% 
#   column_to_rownames('V1')
# drug_targets <- readRDS('drug_targets_space.rds')
# drug_targets_prior <- as.data.frame(drug_targets) %>% rownames_to_column('drug') %>% 
#   gather('variable','Prior knowledge',-drug)
# training_prior <- merged_interactions %>% filter(`Prior knowledge`=='Interaction')
# training_prior <- training_prior %>% select(-`Prior knowledge`,-Inferred) %>% mutate(mode='training_prior')
# drug_targets_prior <- anti_join(drug_targets_prior,training_prior)
# merged_interactions <- merged_interactions %>% select(-`Prior knowledge`)
# merged_interactions <- left_join(merged_interactions,drug_targets_prior)
# merged_interactions <- merged_interactions %>% filter(!is.na(`Prior knowledge`))
# merged_interactions <- merged_interactions %>%
#   mutate(`Prior knowledge`=ifelse(`Prior knowledge`==0,'No interaction','Interaction'))
# confusion <- confusionMatrix(data=as.factor(merged_interactions$Inferred), 
#                              reference = as.factor(merged_interactions$`Prior knowledge`))
# confusion
# plt <- as.data.frame(confusion$table)
# #plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))
# png('drugbank_confusion.png',width=9,height=9,units = "in",res = 600)
# ggplot(plt, aes(Reference,Prediction, fill= Freq)) + 
#   geom_tile() + geom_text(aes(label=Freq)) +
#   scale_fill_gradient2() +
#   labs(x = "Reference",y = "Prediction") + ggtitle('Confusion matrix for similarity in the drug-target space')+
#   theme_minimal(base_family = "serif",base_size = 15)
# dev.off()

### See if I consider as ground-truth cmap plus drugbank
### how precision and recall change with respect to including predictions
### based on how many model they appear
no_models <- 50
merged_interactions_all <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
cmap_prior <- merged_interactions_all %>% select(drug,variable,`Prior knowledge`) %>% unique()
cmap_prior <- cmap_prior %>% filter(`Prior knowledge`=='No interaction') %>% select(-`Prior knowledge`) %>% unique()
drug_targets_prior <- as.data.frame(drug_targets) %>% rownames_to_column('drug') %>% 
  gather('variable','Prior knowledge',-drug)
merged_interactions_all <- merged_interactions_all %>% select(-`Prior knowledge`)
merged_interactions_all <- left_join(merged_interactions_all,drug_targets_prior)
merged_interactions_all <- merged_interactions_all %>% filter(!is.na(`Prior knowledge`))
merged_interactions_all <- merged_interactions_all %>% mutate(Inferred = 1*(Inferred=='Interaction'))
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
  mutate(model_counts = sum(Inferred==1)) %>% ungroup()
merged_interactions_all <- merged_interactions_all %>% mutate(model_counts=model_counts/no_models)
print(all(merged_interactions_all %>% filter(model_counts==0) %>% mutate(logic = (Inferred==`Prior knowledge`)) %>% select(logic)))
# # Get rid of things already in the cmap
# cmap_prior <- cmap_prior %>% mutate(in_camp=1)
# merged_interactions_all <- left_join(merged_interactions_all,cmap_prior)
# merged_interactions_all <- merged_interactions_all %>% filter(is.na(in_camp)) %>% select(-in_camp) %>% unique()


thresholds <- seq(1,no_models)/no_models
prec <- NULL
rec <- NULL
acc <- NULL
F1 <- NULL
pvals <- NULL
for (i in 1:length(thresholds)){
  # For those inferred keep everything above the threshold
  # Those never inferred in the model are indeed never in prior knowledge
  th <- thresholds[i]
  #thresholded_predictions <- distinct(merged_interactions_all %>% filter(model_counts>=th | model_counts==0) %>% 
  #  select("drug","variable","Prior knowledge","Inferred"))
  thresholded_predictions <- distinct(merged_interactions_all %>% mutate(Inferred = ifelse(model_counts>=th,1,0)) %>% 
                                        select("drug","variable","Prior knowledge","Inferred"))
  confusion <- confusionMatrix(data=factor(thresholded_predictions$Inferred,levels = c(1,0)), 
                               reference = factor(thresholded_predictions$`Prior knowledge`,levels = c(1,0)))
  #print(confusion$table)
  #print(confusion$byClass)
  acc[i] <- confusion$overall['Accuracy']
  prec[i] <- confusion$byClass['Precision']
  rec[i] <- confusion$byClass['Recall']
  F1[i] <- confusion$byClass['F1']
  pvals[i] <- confusion$overall['AccuracyPValue']
}
df_res <- data.frame(thresholds=thresholds,accuracy = acc,precision=prec,recall=rec,F1=F1,pvalue = pvals)
df_res <- df_res %>% gather('metric','value',-thresholds)
# visualize
p1 <- ggplot(df_res %>% filter(metric!='pvalue'),
       aes(x=thresholds,y=value,color=metric)) + 
  geom_point() +geom_line(linewidth=1) + 
  xlab('minimum frequency to consider an interaction') + 
  ggtitle('Discovery of new interactions as the frequency of appearing in multiple models increases')+
  theme_minimal() +
  #geom_hline(yintercept = 0.5,linetype='dashed',linewidth=1.5,color='black')+
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24),
        legend.text = element_text(family = 'Arial',size = 24),
        plot.title = element_text(family = 'Arial',size = 20, hjust = 0.5))
print(p1)
png('A375model_interaction_discovery_thresholding.png',width=12,height=12,units = "in",res=600)
print(p1)
dev.off()

ggsave('../MIT/LauffenburgerLab/drugLembasPaper/A375model_interaction_discovery_thresholding.eps',
       plot = p1,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 8,
       units = "in",
       dpi = 600)

## Calculate p.values of individual models
p_vals_one <- NULL
for (i in 1:no_models){
  model_predictions <- distinct(merged_interactions_all %>% filter(model_no==i-1) %>% 
                                        select("drug","variable","Prior knowledge","Inferred"))
  confusion <- confusionMatrix(data=factor(model_predictions$Inferred,levels = c(1,0)), 
                               reference = factor(model_predictions$`Prior knowledge`,levels = c(1,0)))
  p_vals_one[i] <- confusion$overall['AccuracyPValue']
}
hist(p_vals_one)

p2 <- ggplot(df_res %>% filter(metric=='pvalue'),
       aes(x=thresholds,y=value)) + 
  #ylim(c(0,1))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  #ylim(c(0,1))+
  geom_point() +geom_line(linewidth=1,linetype='dashed') + 
  xlab('minimum frequency to consider prediction') + 
  ylab('p-value') +
  ggtitle('P-value for accuracy to be greater than no-information rate')+
  theme_minimal() +
  geom_hline(yintercept = 0.01,linetype='dashed',linewidth=1.5,color='red')+
  annotate('text',x=0.25,y=1e-07,
           label='p-value = 0.01',size=10)+
  geom_hline(yintercept = mean(p_vals_one),linetype='dashed',linewidth=1.5,color='orange')+
  annotate('text',x=0.25,y=10000,
           label='p-value of individual models',size=10)+
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24),
        legend.text = element_text(family = 'Arial',size = 24),
        plot.title = element_text(family = 'Arial',size = 22, hjust = 0.5))
print(p2)
png('A375model_interaction_discovery_thresholding_NIR_pvalues.png',width=12,height=12,units = "in",res=600)
print(p2)
dev.off()

ggsave('../MIT/LauffenburgerLab/drugLembasPaper/A375model_interaction_discovery_thresholding_NIR_pvalues.eps',
       plot = p2,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 8,
       units = "in",
       dpi = 600)

# plot confusion matrix for all ensembles A375-----------------
no_models <- 50
#th <- 40/no_models
merged_interactions_all <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
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
  scale_fill_viridis(trans='log10', limits = c(1,round(max(ensemble_final$Count),-3)),breaks=c(1,10,100,1000,10000)) +
  labs(fill='Count')+
  geom_text(aes(label = paste0('~',paste(Count, Count_se, sep = "\u00B1"))), size = 14) +
  labs(title = "Confusion matrix of inferred and prior knowledge drug-target interactions",
       x = "Inferred",
       y = "Prior knowledge") +
  theme_minimal() +
  theme(text = element_text(family = 'Arial',size=24),
        axis.text = element_text(family = 'Arial',size = 24),
        axis.title = element_text(family = 'Arial',size = 24, face = "bold"),
        legend.text = element_text(family = 'Arial',size = 20),
        plot.title = element_text(family = 'Arial',size = 24, hjust = 0.5))
print(plot_confusion)


png('../MIT/LauffenburgerLab/drugLembasPaper/A375_ensembled_drugtarget_confusion.matrix.png',units = 'in',width = 12,height = 9,res = 600)
print(plot_confusion)
dev.off()

ggsave('../MIT/LauffenburgerLab/drugLembasPaper/A375_ensembled_drugtarget_confusion.matrix.eps',
       plot = plot_confusion,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)
## Plot drug-target matrix filled by frequency
merged_interactions_all <- data.table::fread('Model/CVL1000_Paper/A375_ensembles/all_merged_interactions_drugs.csv',header=T) %>% column_to_rownames('V1')
merged_interactions_all <- merged_interactions_all %>% select("drug","variable","Prior knowledge","Inferred","model_no")
merged_interactions_all <- merged_interactions_all %>% group_by(drug,variable) %>% 
  mutate(model_counts = sum(Inferred=='Interaction')) %>% ungroup()

ann_text <- ensemble_final %>% filter(Var1=='Interaction') %>% select(c('Prior knowledge'='Var2'),Count,Count_sd)
ann_text <- ann_text %>% mutate(Count_se = Count_sd/sqrt(no_models)) %>% select(`Prior knowledge`,Count,Count_se)
total_prior <- distinct(merged_interactions_all %>% select(drug,variable,`Prior knowledge`))
total_prior_negatives <-  sum(total_prior$`Prior knowledge`!='Interaction')
total_prior <- sum(total_prior$`Prior knowledge`=='Interaction')
ann_text <- ann_text %>% mutate(total_prior = ifelse(`Prior knowledge`=='Interaction',total_prior,total_prior_negatives))
ann_text <- ann_text %>% mutate(Count=Count/total_prior) %>% mutate(Count_se=Count_se/total_prior)
ann_text <- ann_text %>% mutate(lab=paste0('~',paste(round(100*Count,2),round(100*Count_se,2), sep = "\u00B1"),'%'))
ann_text <- ann_text %>% mutate(`Prior knowledge`=paste0('Prior : ',`Prior knowledge`))
ann_text <- ann_text %>% select(`Prior knowledge`,lab)
ann_text <- ann_text %>% mutate(lab = ifelse(`Prior knowledge`=='Prior : Interaction',
                                             paste0('retrieved',"\n",lab,"\n",' of prior knowledge'),
                                             paste0(lab,"\n",' of prior negatives')))

p <- ggplot(distinct(merged_interactions_all %>% filter(Inferred=='Interaction') %>%
                  select(drug,variable,`Prior knowledge`,model_counts) %>% mutate(`Prior knowledge`=paste0('Prior : ',`Prior knowledge`))),
       aes(x=model_counts)) +
  geom_histogram(aes(y = ..count..),fill='#0F52BA') +
  scale_x_continuous(limits = c(0,NA))+
  facet_wrap(~`Prior knowledge`,scales = "free") + xlab('number of models')+
  ggtitle('Inferred interactions by the model') +
  #scale_y_log10()+
  #scale_y_continuous(labels = scales::percent_format(accuracy = 1L))+
  theme_pubr(base_family = 'Arial',base_size = 24)+
  theme(panel.grid.major.y = element_line(linetype = 'dashed',colour = 'lightgrey'),
        panel.grid.minor.y = element_line(linetype = 'dashed',colour = 'lightgrey'),
        panel.grid.major.x = element_line(linetype = 'dashed',colour = 'lightgrey'),
        panel.grid.minor.x = element_line(linetype = 'dashed',colour = 'lightgrey'),
        plot.title = element_text(size=20,hjust = 0.5))
p <- p + geom_text(data = ann_text,mapping = aes(x = c(20,30), y = c(280,5200), label = lab),
                   size=9
)
print(p)

png('../MIT/LauffenburgerLab/drugLembasPaper/A375_barplot_confusion_2.png',units = 'in',width = 12,height = 9,res = 600)
print(p)
dev.off()

ggsave('../MIT/LauffenburgerLab/drugLembasPaper/A375_barplot_confusion_2.eps',
       plot = p,
       device = cairo_ps,
       scale = 1,
       width = 12,
       height = 9,
       units = "in",
       dpi = 600)
