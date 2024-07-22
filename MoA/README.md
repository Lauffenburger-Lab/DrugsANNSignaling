## Scripts for exploring the MoA of off-target effects of drugs.
Here we deposit code to perform mechanism of action MoA exploration due to off-target effects of a drug on a transcription factor in a specific sample.

### Scripts to be used (**with appropriate input arguments**) in a user case study:
**Run all these after the having run the scripts in the `learning` and `postprocessing` folders:**
1. inferEnsembleScoreCaseStudy.py: **Generally for any model and drug module you trained and used** you may run this script first to calculate integrated gradient scores and infer interactions based on them, using an ensemble of multiple models, and make TF activity predictions for masking out interactions for different thresholds of absolute gradient scores (and save them).
2. InferInteractionScoresLinAlg.py: **ONLY for the specific linear drug module used in the manuscript**, you may run this script first to instead of 1. to calculate drug-target interaction scores.
3. InferDTICaseStudy.py: After inferring drug-target interaction scores you may run this script second to infer interactions and calculate the error of the model as you remove interactions based on their interaction score, previously calculated.
4. inferMoACaseStudy.py: Finally using the files generated from the scripts above, you may run this to infer and save the mechanism of action networks that explain the off-target effect of a drug on a specific TF.
5. DrugTargetInteractionSignCaseStudy.py: You can examine the sing of the interaction between a target in the constructed MoA network and the drug of interest using this script. It will print out and plot information about the activity of the target node in the initial perturbation, the pseudo-steady state, as well as the sign of the interaction score.

### Scripts used in the original [manuscript](https://doi.org/10.1016/j.isci.2024.109509):
1. inferEnsembleTargetScores.py: This script contains code to calculate integrated gradient scores and infer interactions based on them, using an ensemble of multiple models.
2. inferTFerrorBasedThreshold.py: This script contains code to calculate the average(global) error (or specifically for a TF of interest like FOXM1) and the corresponding gradient score threshold for inferring drug-target interactions for each drug.
3. inferDrugTargetsBasedonTFpredictionsALL.py: This script contains code to infer all drug-target interactions of all drugs in a cell line and save them.
4. inferMoA.py: This script contains code to infer and save the mechanism of action networks that explain the off-target effect of a drug on a specific TF.
5. InSilicoSimulation.py:  This script contains code to perform the in-silico KOs of signaling nodes and save the results (the results are used for Figures 4C-4D and supplementary Figure 10)
4. bionetwork.py: Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the functions necessary to build and use a LEMBAS signaling model.
7. bionetworkWithDrugs.py: Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the same functions as before, but we have added a function to build and construct a Drug Layer.