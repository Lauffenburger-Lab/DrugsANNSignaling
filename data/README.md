## Data: Here you should put your raw data.
Most of the raw data are too big to be uploaded but they can be accessed in the:
* [L1000 dataset:GSE92742](https://www-ncbi-nlm-nih-gov.libproxy.mit.edu/geo/query/acc.cgi?acc=GSE92742)
* [siRNA dataset:GSE31534](https://www-ncbi-nlm-nih-gov.libproxy.mit.edu/geo/query/acc.cgi?acc=GSE31534)
* [Broad’s Institute Repurposing Hub](https://clue.io/repurposing)
* [DrugBank](https://go.drugbank.com/)

Note: In the [CLUE Glossary](https://clue.io/connectopedia/glossary) you can find definitions about many variables and identifiers in the L1000 data.

## Data description
1. GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz : The level 3 normalized gene expression data of experiments of drugs tested on specific cancer cell lines. Replicates are NOT aggregated at this level.
2. GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz : The level 5 transformed (z-scored) gene expression data of experiments of drugs tested on specific cancer cell lines. Replicates are aggregated at this level
3. GSE92742_Broad_LINCS_gene_info.txt.gz : Meta-data information about the measured and inferred (in the platform) genes.
4. GSE92742_Broad_LINCS_inst_info.txt.gz : Meta-data information about each sample at the individual technical replicates level. Each sample is identified by a unique inst_id.
5. GSE92742_Broad_LINCS_sig_info.txt.gz : Meta-data information about each experiment/sample. Each experiment/sample is identified by a unique sig_id after replicate aggregation.
6. GSE92742_Broad_LINCS_pert_info.txt.gz : Meta-data information about each perturbation, e.g. the drugs name, SMILE, dose, time, etc.
7. GSE92742_Broad_LINCS_sig_metrics.txt.gz : Quality metrics' information about each individual experiment.
