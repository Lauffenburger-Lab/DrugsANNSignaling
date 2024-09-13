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
3. extractRL.py: Script to extract receptor-ligand (RL) interactions.
4. extractPKN.py: Script to extract prior knowledge network (PKN) of protein-protein interactions.
5. ConnectData2PKN.py: Script to keep TFs, drugs, and targets that can be connected to the PKN, and save parts of the PKN to forcefully keep even after trimming.
6. trimPKN.py: Script to trim the prior knowledge network (PKN).
7. CSV_TSV_to_RDS.R : **Simple** R script that converts CSV or TSV files to an RDS file. All of them need to have rownames and column names. The first column will be considered that it contains the rownames. **The saved RDS file is going to contain a data frame, not a matrix/array.**
8. RDS_to_CSV_TSV.R : **Simple** R script that converts RDS file to CSV or TSV file.
9. MatrixLongFormatConversion.py: Python script to transform CSV/TSV file that contains a data frame in matrix format to a CSV/TSV file in long format (and vice versa).

## Folder structure
1. preprocessed_data: It contains the pre-processed data from the above scripts to be used later for training models and other downstream analyses of the study.
	* **ChemicalSims** : Script and data to calculate the chemical simililarity of drugs.
	* **PKN** : Prior knowledge signaling network and each trimmed version.
	* **RL** : Receptor-ligand interactions data.
	* **TF_activities** : Inferred TF activity data.
	* **TrainingValidationData** : Training conditions and drugs that will be used for validation

The produced data are used to create the figures stored in the figures folder and the results in the results folder.
Also these data were used to train machine learning models.
