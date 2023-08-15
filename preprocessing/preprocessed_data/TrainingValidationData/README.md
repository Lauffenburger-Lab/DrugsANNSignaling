## TrainingValidationData
Folder containing the final data/conditions used for training and validating the models.

## Files in this folder
1. L1000_lvl3_allcells-conditions_drugs.tsv : Final conditions, in terms of normalized dose of drug tested on a cell line. 
2. L1000_lvl3_allcells-drugs_targets.tsv : Final prior drug-target interactions knowledge that will be used for training models.
3. L1000_lvl3-conditions_cells.tsv : Binary matrix mapping each perturbation to its corresponding cell line.
4. smiles_and_sigs_to_drop.csv : SMILEs of drugs that were held out for valdiation purposes and identifiers of the corresponding samples (sig_id).
