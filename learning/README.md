## Machine and Deep learning algorithms of the project
This folder contains scripts to train machine learning models to predict the TF activity of drug perturbations, as well as some scripts related to tuning some parameters of the models.

### Scripts to be used (**with appropriate input arguments**) in a user case study:
1. macrophageNetTrainOne.py : First run this script with the data of your choice to train a whole cell line-specific model (drug layer + LEMBAS signaling). With the argument 'no' a numeric id is given to the trained models and it can be used to parallelize the training of the ensemble of models.
2. CellLineSpecificEvalEnsembleALL.py: Then run this script to evaluate how well the model fitted the activity of each TF (the training performance) 

### Scripts used in the original [manuscript](https://doi.org/10.1016/j.isci.2024.109509):
1. bionetwork.py : Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the functions necessary to build and train a LEMBAS signaling model.
2. bionetworkWithDrugs.py : Adapted from https://github.com/Lauffenburger-Lab/LEMBAS , this script contains the same functions as before, but we have added a function to build and construct a Drug Layer.
3. macrophageNetTrainOne.py : This script contains code to train a whole model (drug layer + LEMBAS signaling) using the data of only one cell line. With the argument 'no' a numeric id is given to the trained models and it can be used to parallelize the training of the ensemble of models.
4. macrophageNet.py : This script contains code to re-train only the signaling model part on the rest of the cell lines we are going to use for validation. With the argument 'pretrained' we specify which drug layer from which model of the ensemble to use.
5. macrophageNetEvalEnsemble.py : This script contains code to evaluate an ensemble of models (individual model performance per TF and ensemble performance).
6. EstimateRequiresNumberOfEnsembles.py : This script contains code to estimate how the performance of the model change as we increase the number of models in the ensemble.
7. LearningParamsTuning.py : Contains code to train a model for a few different learning rates to observe how learning rate affects training and validation performance.
8. activationFunctions.py : Contains a couple of activation functions that can be used in LEMBAS/
9. saveSimulations.py : Contains code to save progress and results of training a model.
10. plotting.py : Contains code to plot the results of training and validation.
11. vanillaANNSignal.py : Contains code to train a simple ANN to predict TF activity from the input signal of the drug layer, without any constraint using some PKN. 
12. vanillaANNRandomSignal.py : Contains code to perform the previous task, but the input signal is randomly shuffled to create a null/random model.
13. vanillaModels.py : Contains code to train standard machine learning models such as ANNs, GCNNs, KNN, and Support Vector Regression (SVRs).
14. vanillaModelsEnsembles.py : Contains code to train ensembled of standard machine learning models such as ANNs, GCNNs, KNN, and Support Vector Regression (SVRs).
15. features.py : Contains code to create molecular graphs from the SMILEs of drugs.
16. graph_layers.py : Contains code necessary to build Graph Convolutional Layers for the GCNN models.
