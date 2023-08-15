## TF_activities
 It contains files related to the TF activity inference process.

## Files in this folder 
Some files are in .rds format for saving storage space, and they can be read fast and easily in R with the readRDS() command.
If a file is not present, it means it is too big to be uploaded in GitHub and can be provided upon reasonable request.

1. TFs_random_variances_all.rds : Random replicates variance, derived from random permutations of rows 100 times. It is a file useful for reproducability of the study. 
2. filtered_tfs_from_repls.rds : TFs that would be filtered out in this study. All TFs survived this process in this study.
3. difference_means_of_variances.rds : Difference in means of the variacne distributions of TFs, for the random (NULL) distribution and our data.
4. ks_tests_repls.rds : The result of KS test, comparing our data with the NULL distribution.
5. random_correlations.rds :
6. TFs_cor_repls.rds :
7. l1000_allgenes_lvl3_tfs.tsv :
8. tfs_targetd_alls_genes_lvl3.tsv :
9. Trimmed_l1000_allgenes_lvl3_tfs.tsv :
10. TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv : 
