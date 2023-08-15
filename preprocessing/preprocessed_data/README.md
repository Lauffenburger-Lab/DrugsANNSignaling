## Preprocessed data
 It contains the split and pre-processed data used to train models and evaluate them. 
 Additionally, it contains the retrieved:
 1. Drug-target interaction data.
 2. Protein-protein interactions and signaling network
 3. And all the training data aggregated.

## Folder structure
1. **ChemicalSims** : Script and data to calculate the chemical similarity of drugs.
2. **PKN** : Prior knowledge signaling network and each trimmed version.
3. **RL** : Receptor-ligand interactions data.
4. **TF_activities** : Inferred TF activity data.
5. **TrainingValidationData** : Training conditions and drugs that will be used for validation.

## Files in this folder 
These files are in .rds format for saving storage space, and they can be read fast and easily in R with the readRDS() command.
If a file is not present, it means it is too big to be uploaded in GitHub and can be provided upon reasonable request.

1. sig_info_5_with_exempler.rds : All available metadata of drug perturbations and control conditions in the L1000 (sample information).
2. all_cmap_sigs_with_pert_info.rds : All available sample information combined with perturbation information, i.e. drug name, SMILE and other chemical IDs etc. 
3. all_cmapdrugs_moa+target_v1.rds : All drugs in the Connectivity Map (CMap) dataset, together with their known targets retrieved from the Broad's Repurposing Hub.
4. l1000_drugs_with_targets_all.rds : All 1223 unique drugs available in the L1000, that have some known target and their chemical structure is documented.
5. l1000_drugs_with_targets_exemplar.rds : The 1224 L1000 drugs that are present in exemplar perturbations and have known targets. 
6. l1000_drugs_with_targets_unnested_exemplar.rds : The 1015 drugs with targets in the prior knowledge signaling network. The drug-target interaction information is also given in a long format (unnested).
7. l1000_all_genes_lvl3_drugs_with_targets_exemplar.rds : The level 3 normalized gene expression of the 10,174 measured and imputed genes of samples perturbed by the 1015 drugs mentioned before.
8. TrimmedFinal_lvl3_allgenes_all_conditions.rds : The conditions used in the L1000, with drugs having a target in the trimmed prior knowledge network.
9. TrimmedFinal_allgenes_lvl3_all_conditions_with_targets.rds : The conditions used in the L1000, with drugs having a target in the trimmed prior knowledge network, after manually adding DMSO information and controls. Also, we kept perturbations were cell lines were treated for 6 hours.
10. TrimmedFinal_lvl3_conditions_for_pairedcells.rds : The final conditions used for the 9 cell lines that have in common 233 drugs.
