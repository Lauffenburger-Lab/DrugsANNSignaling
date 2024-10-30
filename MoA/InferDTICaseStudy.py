import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Drug-Targets for drugs of interest')
parser.add_argument('--ensembles_path', action='store', default="../../results/case_study/models/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_model_modeltype4_a375_case_study")
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
parser.add_argument('--drugInputFile', action='store',default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_A375-conditions_drugs.tsv')
parser.add_argument('--drugTargetsFile', action='store',default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_A375-drugs_targets.tsv')
parser.add_argument('--TFOutFile', action='store',default='../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv')
parser.add_argument('--drugSimilarityFile', action='store',default='../preprocessing/preprocessed_data/ChemicalSims/lvl3_similarities_A375.csv')
parser.add_argument('--res_dir', action='store',default='../results/A375_ensembles/')
parser.add_argument('--Y_ALL_path', action='store',default='../results/A375_ensembles/Y_ALL_testing.pt')
parser.add_argument('--Y_ALL_masked_path', action='store',default='../results/A375_ensembles/Y_ALL_masked_testing.pt')
parser.add_argument('--interactionsPath', action='store',default='../../results/case_study/InteractionScores/')
parser.add_argument('--error_threshold', action='store',default=0.25)

args = parser.parse_args()
ensembles_path = args.ensembles_path
numberOfModels = int(args.numberOfModels)
ConvertToEmpProb = args.ConvertToEmpProb
if type(ConvertToEmpProb) == str:
    ConvertToEmpProb = eval(ConvertToEmpProb)
inputPattern = args.inputPattern
inputPath = ensembles_path + inputPattern
drugInputFile = args.drugInputFile
drugTargetsFile = args.drugTargetsFile
TFOutFile = args.TFOutFile
drugSimilarityFile = args.drugSimilarityFile
res_dir = args.res_dir
Y_ALL_masked_path = args.Y_ALL_masked_path
Y_ALL_path = args.Y_ALL_path
interactionsPath = args.interactionsPath

def inferDrugTarget(interactions, global_thresholds, prior_mask, drugInput, drugTargets, model_no,
                    grad_thresholds=list(np.logspace(-3.5, 3.5, num=45)), thresh=0.25):
    global_thresholds = global_thresholds.detach().numpy()
    scores = torch.abs(torch.tensor(interactions.values))
    grad_thresh = global_thresholds[:,model_no:model_no+1]
    mask = 1.0 * (scores.detach().numpy() >= grad_thresh)
    
    df_mask = pd.DataFrame(prior_mask.numpy())
    df_mask.index = drugInput.columns
    df_mask.columns = drugTargets.columns
    df_mask.reset_index(inplace=True)
    df_mask = df_mask.rename(columns={'index': 'drug'})
    df_mask = pd.melt(df_mask, id_vars=['drug'])
    df_mask['value'] = df_mask['value'].transform(lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
    df_mask = df_mask.rename(columns={'value': 'Prior knowledge'})
    
    df_masked_scores = pd.DataFrame(mask)
    df_masked_scores.index = drugInput.columns
    df_masked_scores.columns = drugTargets.columns
    df_masked_scores.reset_index(inplace=True)
    df_masked_scores = df_masked_scores.rename(columns={'index': 'drug'})
    df_masked_scores = pd.melt(df_masked_scores, id_vars=['drug'])
    df_masked_scores['value'] = df_masked_scores['value'].transform(
        lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
    df_masked_scores = df_masked_scores.rename(columns={'value': 'Inferred'})
    merged_df = pd.merge(df_mask, df_masked_scores, how='left', left_on=['drug', 'variable'],
                         right_on=['drug', 'variable'])
    return merged_df

networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations=120, clipping=1, targetPrecision=1e-6, leak=0.01)

drugInput = pd.read_csv(drugInputFile, sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pd.read_csv(drugTargetsFile, sep='\t', low_memory=False, index_col=0)
TFOutput = pd.read_csv(TFOutFile, sep='\t', low_memory=False, index_col=0)
TFOutput = TFOutput.loc[drugInput.index,:]

druginName = drugInput.columns.values
inName = drugTargets.columns.values
outconds = drugInput.index.values
outName = TFOutput.columns.values
inName = np.intersect1d(nodeNames, inName)
outName = np.intersect1d(nodeNames, outName)
inNameGene = [uniprot2gene[x] for x in inName]
outNameGene = [uniprot2gene[x] for x in outName]
TFOutput = TFOutput.loc[outconds, outName]
drugTargets = drugTargets.loc[:, inName]
if ConvertToEmpProb:
    print2log('Converted to Empirical probabilities')
    TFOutput = 1 / (1 + np.exp(-TFOutput))

drugTargets = drugTargets.loc[drugInput.columns.values, :]
drugSim = pd.read_csv(drugSimilarityFile, index_col=0)
drugSim = drugSim.loc[drugInput.columns.values, drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

Y = torch.tensor(TFOutput.values, dtype=torch.double)

thresholds = list(np.logspace(-3.5, 3.5, num=50))
thresh = args.error_threshold
Y_ALL = torch.load(Y_ALL_path)
Y_ALL_masked = torch.load(Y_ALL_masked_path)
global_threshold_all = np.zeros((drugInput.shape[1], numberOfModels))

print2log('Begin calculating thresholds for each drug')
drug_counts = np.count_nonzero(drugInput, axis=0)

for i in range(numberOfModels):
    for j in range(drugInput.shape[0]):
        drug_ind = np.where(drugInput.iloc[j, :] != 0)[0]
        global_error = torch.mean(torch.abs(Y_ALL_masked[i, :, j, :] - Y[j, :]).squeeze(), 1).detach().squeeze().numpy()
        unmasked_global_error = torch.mean(torch.abs(Y_ALL[i, j, :] - Y[j, :])).detach().squeeze().numpy()
        global_perc_err = (global_error - unmasked_global_error) / (global_error[-1] - unmasked_global_error)
        if np.where(global_perc_err <= thresh)[0].shape[0] > 0:
            global_thresh_selection = np.where(global_perc_err <= thresh)[0][-1]
        else:
            global_thresh_selection = 0
        global_threshold_all[drug_ind, i] += thresholds[global_thresh_selection]
    print2log('Finished calculating thresholds for model %s' % i)

global_threshold_all = global_threshold_all / drug_counts.reshape(drugInput.shape[1], 1)
torch.save(torch.tensor(global_threshold_all), res_dir + 'all_drugs_global_thresholds_testing.pt')
global_thresholds_df = pd.DataFrame(global_threshold_all)
global_thresholds_df.index = drugInput.columns.values
global_thresholds_df.to_csv(res_dir + 'all_drugs_global_thresholds_testing.csv')

print2log('Start inferring interaction with the error-based method')
drug_ratioMatrix = torch.arange(0, 1.01, 0.01).T.repeat(drugInput.shape[1], 1).double()

for i in range(numberOfModels):
    model = torch.load(inputPath + str(i) + ".pt")
    model.eval()
    interactions = pd.read_csv(interactionsPath + 'interactionScores_%s.csv' % i, index_col=0)
    merged_interactions = inferDrugTarget(interactions, torch.tensor(global_threshold_all), model.drugLayer.mask.T.detach(), drugInput, drugTargets, i, thresh=thresh)
    merged_interactions['model_no'] = i
    if i == 0:
        final_interactions = merged_interactions
    else:
        final_interactions = pd.concat([final_interactions, merged_interactions], axis=0)

final_interactions.to_csv(res_dir + "A375_masked_threshold_all_drugs_global_thresh_final_interactions.csv")
print2log('Completed all thresholding')