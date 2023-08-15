## Data and Code for chemical similarity calculations 

## Scripts
1. smiles_similarity_ecfp4.py : Python script to calculate the chemical similarity between drugs. It has the option to get a list of SMILES and calculate all pairwise similarities in the list (**what was done**), or 2 separate lists and calculate similarities between these lists.
2. proc.sh : Bash script to execute the Python script using a file to pass automatically the appropriate arguments required by the script. You can also run the Python script by itself and follow the instructions that will be printed on your screen.

## Files in this folder
1. input_all : Text file having the arguments that will be passed in the Python script.
2. lvl3_smiles.csv : A CSV file containing in one column the list of SMILEs of drugs that will be used in chemical similarity calculation. **The column needs to be named: smiles**.
3. out_lvl3_similaritiess.csv : The output CSV file of the Python script. It is a matrix in CSV format containing the chemical similarity between drugs, and the columns and rows are named with the SMILE of each drug, which can be used as an identifier later on.
4. allcells_test_smiles.csv : This is an output of the preprocessing_L1000.R script. It contains the SMILEs of the drugs that are candidates for test drugs as their maximum chemical similarity with any other drug is less than 0.6.
