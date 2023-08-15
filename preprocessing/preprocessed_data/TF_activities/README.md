## TF_activities
 It contains files related to the TF activity inference process.

## Files in this folder 
Some files are in .rds format for saving storage space, and they can be read fast and easily in R with the readRDS() command.
If a file is not present, it means it is too big to be uploaded to GitHub and can be provided upon reasonable request.

1. TFs_random_variances_all.rds : Random replicates variance, derived from random permutations of rows 100 times. It is a file useful for the reproducibility of the study. 
2. filtered_tfs_from_repls.rds : TFs that would be filtered out in this study. All TFs survived this process in this study.
3. difference_means_of_variances.rds : Difference in means of the variance distributions of TFs, for the random (NULL) distribution and our data.
4. ks_tests_repls.rds : The result of the KS test, comparing our data with the NULL distribution.
5. random_correlations.rds : Null distribution of correlations between random pairs of samples in our data.
6. TFs_cor_repls.rds : Correlation between technical replicates, together with the probability of observing an equally high or higher correlation in our data by chance.
7. l1000_allgenes_lvl3_tfs.tsv : Inferred TF activity after filtering and merging replicates, but before the trimming of the PKN and before finding cell lines with a sufficiently high number of drugs tested on them.
8. tfs_targetd_alls_genes_lvl3.tsv : TFs that are directly targeted by drugs in our data.
9. Trimmed_l1000_allgenes_lvl3_tfs.tsv : TF activity data after trimming the network.
10. TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv : TF activity data after trimming the network and final cleaning and filtering of conditions.
