## Algorithms for the pre-processing of the raw data (specialized for the L1000 dataset)
1. retrieveDrugTargetInfo.R: Script to retrieve drug-target interactions, merge with the L1000 dataset, and retrieve also level-3 gene expression data.
2. drug_sigs_for_cell.R: Script with function to retrieve samples metadat and adds quality score to each signature and for each drug/cell line combination selects the signature with the highest quality.
3. inferTFactivity.R: Script to infer TF activity and pre-process to remove samples and TFs of bad quality.
4. preprocessing_L1000.R: Script to pre-process the remaining data and get the final conditions, as well as split drugs into training and test drugs.
5. extractRL.py: Script to extract receptor-ligand (RL) interactions.
6. extractPKN.py: Script to extract prior knowledge network (PKN) of protein-protein interactions.
7. trimPKN.py: Script to trim the prior knowledge network (PKN).

## Scripts to be used (**with appropriate input arguments**) in a user case study:
1. inferTFactivityCaseStudy.py: Run this to infer TF activity for some samples for your own case study.
2. PreprocessTFactivityCaseStudy.R: Run this (optionally) to pre-process and remove samples and TFs of bad quality.

## Folder structure
1. preprocessed_data: It contains the pre-processed data from the above scripts to be used later for training models and other downstream analyses of the study.
	* **ChemicalSims** : Script and data to calculate the chemical simililarity of drugs.
	* **PKN** : Prior knowledge signaling network and each trimmed version.
	* **RL** : Receptor-ligand interactions data.
	* **TF_activities** : Inferred TF activity data.
	* **TrainingValidationData** : Training conditions and drugs that will be used for validation

The produced data are used to create the figures stored in the figures folder and the results in the results folder.
Also these data were used to train machine learning models.
