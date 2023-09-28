## Post-processing of model's predictions results.

This folder contains code for analyzing models' results and producing some of the article's figures and supplementary figures.
This folder does not contain code regarding inferring the network describing the MoA of a frug due to off-target effects (scripts can be found in the MoA folder).
While it contains some Python scripts for inferring dug-target interactions the rest can be found again in the MoA folder.

Scripts list:
1. EstimateRequiredEnsembles.R: Script containing code to analyze how the performance of the model changes with an increasing number of ensembles. (produces figures to estimate the number of ensembles VS performance)
2. EnsembleEval.R: Script containing code to evaluate ensembles of models and compare their performance with that of individual models (produces supplementary figures 1,5)
3. CompareResultsVanillasEnsembles.R: Script containing code to compare our model with other machine learning approaches and generally analyze its performance (produces figures 1B-1C)
4. regularizationEval.R: Script containing code to evaluate performance as a function of the Drug Module's regularization (produces figure 1F and supplementary figure 13)
5. regularization_discovery_tradeoff.R: Script containing code to evaluate the trade-off between regularization and the rate of discoveries of new interactions (produces figure 1E and supplementary figures 6A,6B)
6. DTI_evaluation.R: Script containing code to augment prior knowledge of drug-target interactions and evaluate model's performance in inferring interactions (produces Figure 2B-2E and can be used for supplementary figures 8,9)
7. chooseTFsWithOffTargets.R: Script containing code to choose a TF affected by a specific drug in a sample due to large off-target effects (produces figure 3B)
8. ExternalGEOStudy.R: Script containing code to retrieve data for the analysis done using an external GEO dataset of siRNA experiments (produces figure 4A)
9. InSilicoValidation.R: Script containing code to validate in-silico whether the model inferred Lestautrinib-CDK2 interaction and perform in-silico KO experiments (produces figure 4B-4D and supplementary figures 10A-10B)
10. viabilityTargets.R: Script containing code to train Radom Forest models using known and inferred targets' information of each drug as input (produces supplementary figure 11,12)
11. inferOffTargetEffectDelta.py: Script containing code to calculate the metric that acts as a proxy for having a large off-target effect and calculating training ensemble performance (results will be used for figure 3 and for inferring MoA)
12. inferGradScoresVSRegularizations.py: Script containing code to produce gradient scores for regularization VS inference analysis for the error-based method. Additionally, it saves ensembled predictions and masked predictions for multiple gradient thresholds (will be used to produce supplementary figure 6)
13. inferDrugTargetsRegularizationVSInference.py: Script containing code to finally produce results for the regularization VS inference analysis (will be used to produce supplementary figure 6)
14. figure2_error_vs_gradscore.py: Script containing code to calculate global gradient thresholds to consider interactions for each drug (produces results to be used for figure 2A)
15. bionetwork.py: Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the functions necessary to build and train a LEMBAS signaling model.
16. bionetworkWithDrugs.py: Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the same functions as before, but we have added a function to build and construct a Drug Layer.