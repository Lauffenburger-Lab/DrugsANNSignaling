## Scripts for exploring the MoA of off-target effects of drugs.
Here we deposit code to perform mechanism of action MoA exploration due to off-target effects of a drug on a transcription factor in a specific sample.

Scripts list:
1. inferEnsembleTargetScores.py: This script contains code to calculate integrated gradient scores and infer interactions based on them, using an ensemble of multiple models.
2. inferTFerrorBasedThreshold.py: This script contains code to calculate the average(global) error (or specifically for a TF of interest like FOXM1) and the corresponding gradient score threshold for inferring drug-target interactions for each drug.
3. inferDrugTargetsBasedonTFpredictionsALL.py: This script contains code to infer all drug-target interactions of all drugs in a cell line and save them.
4. inferMoA.py: This script contains code to infer and save the mechanism of action networks that explain the off-target effect of a drug on a specific TF.
5. InSilicoSimulation.py:  This script contains code to perform the in-silico KOs of signaling nodes and save the results (the results are used for Figures 4C-4D and supplementary Figure 10)
4. bionetwork.py: Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the functions necessary to build and use a LEMBAS signaling model.
7. bionetworkWithDrugs.py: Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the same functions as before, but we have added a function to build and construct a Drug Layer.