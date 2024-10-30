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
parser.add_argument('--ensembles_path', action='store', default="../results/A375_ensembles/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_model_modeltype4_model")
parser.add_argument('--cell_line', action='store', default="A375")
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
args = parser.parse_args()
ensembles_path = args.ensembles_path
cell = args.cell_line
numberOfModels = int(args.numberOfModels)
ConvertToEmpProb = args.ConvertToEmpProb
if type(ConvertToEmpProb) == str :
    ConvertToEmpProb = eval(ConvertToEmpProb)
inputPattern = args.inputPattern
inputPattern = cell + "_" + inputPattern
inputPath = ensembles_path + inputPattern

def inferDrugTarget(interactions, global_thresholds, prior_mask, drugInput, drugTargets, model_no,
                    grad_thresholds=list(np.logspace(-3.5, 3.5, num=45)), thresh=0.25):
    global_thresholds = global_thresholds.detach().numpy()
    # interactions = interactions[interactions.index.values!='CS(C)=O']
    scores = torch.abs(torch.tensor(interactions.values))
    grad_thresh = global_thresholds[:,model_no:model_no+1]
    # Get new mask with drug-targets interactions
    mask = 1.0 * (scores.detach().numpy() >= grad_thresh)
    #mask = mask.detach().numpy()

    # Create merged dataframe
    df_mask = pd.DataFrame(prior_mask.numpy())
    df_mask.index = drugInput.columns
    df_mask.columns = drugTargets.columns
    df_mask.reset_index(inplace=True)
    df_mask = df_mask.rename(columns={'index': 'drug'})
    df_mask = pd.melt(df_mask, id_vars=['drug'])
    df_mask['value'] = df_mask['value'].transform(lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
    df_mask = df_mask.rename(columns={'value': 'Prior knowledge'})
    # New interactios
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
    return (merged_df)

#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
spectralCapacity = np.exp(np.log(1e-2)/bionetParams['iterations'])
### Load the data
drugInput = pd.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pd.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput = pd.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
cellInput = pd.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput=TFOutput.loc[cellInput[cellInput[cell]==1].index,:]
drugInput=drugInput.loc[cellInput[cellInput[cell]==1].index,:]
#Subset input and output to intersecting nodes
druginName = drugInput.columns.values
inName = drugTargets.columns.values
outconds = drugInput.index.values
outName = TFOutput.columns.values
inName = np.intersect1d(nodeNames, inName)
outName = np.intersect1d(nodeNames, outName)
inNameGene = [uniprot2gene[x] for x in inName]
outNameGene = [uniprot2gene[x] for x in outName]
TFOutput = TFOutput.loc[outconds,outName]
drugTargets = drugTargets.loc[:,inName]
if ConvertToEmpProb==True:
    print2log('Convereted to Emprirical probabilities')
    TFOutput = 1/(1+np.exp(-TFOutput))
#make sure they are on the same order
drugTargets = drugTargets.loc[drugInput.columns.values,:]
drugSim = pd.read_csv('../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv',index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

#drugInput = drugInput[drugInput.loc[:,'CS(C)=O']==0]
dmso_ind = np.where(drugInput.columns.values=='CS(C)=O')[0][0]
#all_drugs = list(drugInput.columns.values)

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

### Define gradient thresholds
thresholds = list(np.logspace(-3.5, 3.5, num=50))
# thresholds = list(np.logspace(-5, 5, num=50))
thresh = 0.25
Y_ALL = torch.load(ensembles_path+'preds/Y_'+cell+'_ALL.pt')
Y_ALL_masked = torch.load(ensembles_path+'preds/Y_'+cell+'_ALL_MASKED.pt')
#no_dmso_samples =  drugInput[drugInput.loc[:,'CS(C)=O']==0]
#global_threshold_all = np.zeros((no_dmso_samples.shape[0],numberOfModels))
global_threshold_all = np.zeros((drugInput.shape[1],numberOfModels))
all_drugs = []
ind = 0

print2log('Begin caclulating thresholds for each drug')
dmso_thresh_array = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),numberOfModels))
dmso_counter = 0
for j in range(X.shape[0]):
    if drugInput.iloc[j,dmso_ind] == 0:
        all_drugs.append(drugInput.columns.values[np.where(drugInput.iloc[j,:]!=0)[0]][0])
        l = []
        for i in range(numberOfModels):
            # Calculate error seperetaly
            global_error = torch.mean(torch.abs(Y_ALL_masked[i,:,j,:]-Y[j,:]).squeeze(),1).detach().squeeze().numpy()
            unmasked_global_error = torch.mean(torch.abs(Y_ALL[i,j, :] - Y[j, :])).detach().squeeze().numpy()
            global_perc_err = (global_error-unmasked_global_error)/(global_error[-1]-unmasked_global_error)
            if np.where(global_perc_err <= thresh)[0].shape[0]>0:
                global_thresh_selection = np.where(global_perc_err <= thresh)[0][-1]
            else:
                global_thresh_selection = 0
            l.append(thresholds[global_thresh_selection])
        global_threshold_all[ind,:] = np.array(l)
        ind = ind +1
    else:
        dmso_thresh = []
        for i in range(numberOfModels):
            # Calculate error seperetaly
            global_error = torch.mean(torch.abs(Y_ALL_masked[i,:,j,:]-Y[j,:]).squeeze(),1).detach().squeeze().numpy()
            unmasked_global_error = torch.mean(torch.abs(Y_ALL[i,j, :] - Y[j, :])).detach().squeeze().numpy()
            global_perc_err = (global_error-unmasked_global_error)/(global_error[-1]-unmasked_global_error)
            if np.where(global_perc_err <= thresh)[0].shape[0]>0:
                global_thresh_selection = np.where(global_perc_err <= thresh)[0][-1]
            else:
                global_thresh_selection = 0
            dmso_thresh.append(thresholds[global_thresh_selection])
        dmso_thresh_array[dmso_counter,:] = np.array(dmso_thresh)
        if dmso_counter==0:
            final_dmso_ind = ind
            global_threshold_all[ind,:] = 0.
            all_drugs.append(drugInput.columns.values[np.where(drugInput.iloc[j,:]!=0)[0]][0])
            ind = ind+1
        dmso_counter = dmso_counter + 1   
    print2log('Finished drug %s'%j)
global_threshold_all[final_dmso_ind,:] = np.mean(dmso_thresh_array,0)
torch.save(torch.tensor(global_threshold_all),ensembles_path +'all_drugs_global_thresholds.pt')
global_thresholds_df = pd.DataFrame(global_threshold_all)
global_thresholds_df.index = all_drugs
global_thresholds_df.to_csv(ensembles_path +'all_drugs_global_thresholds.csv')

#%%
print2log('Start infering interaction with the error-based method')
drug_ratioMatrix = torch.arange(0, 1.01, 0.01)
drug_ratioMatrix = drug_ratioMatrix.T.repeat(drugInput.shape[1],1)
drug_ratioMatrix = drug_ratioMatrix.double()
### Get drug-target interactions
for i in range(numberOfModels):#range(1):
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()
    # # Get range of input activity for different concetrations of drugs
    # X_binned =  drug_ratioMatrix * X.unsqueeze(2)
    # Xin_range = torch.zeros(X_binned.shape[1],drugTargets.shape[1]).double()
    # drug_sum = torch.zeros(X_binned.shape[1]).double()
    # for k in range(X_binned.shape[0]):
    #     Xin_binned =  model.drugLayer(X_binned[k,:,:].T)
    #     drug_ind = torch.where(X[k,:]!=0)[0]
    #     Xin_range[drug_ind,:] = Xin_range[drug_ind,:] + torch.abs(torch.max(Xin_binned,0)[0] - torch.min(Xin_binned,0)[0]).unsqueeze(0)
    #     drug_sum[drug_ind] = drug_sum[drug_ind] + 1.0
    # Xin_range = Xin_range.T / drug_sum
    # Xin_range = Xin_range.T
    
    ### Get drug-target interactions
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/l1000_modeltype4_lamda6_' + cell + '_interactionScores_%s.csv' % i,index_col=0)
    drugs = interactions.index.values
    target_space = interactions.columns.values
    # interactions = torch.tensor(interactions.values) * Xin_range
    # interactions = pd.DataFrame(interactions.detach().numpy())
    # interactions.index =drugs
    # interactions = interactions.T
    # interactions.index = target_space
    # interactions = interactions.T
    # merged_interactions = inferDrugTarget(interactions,
    #                                       torch.tensor(global_threshold_all),
    #                                       model.drugLayer.mask.T.detach()[no_dmso_samples.columns != 'CS(C)=O',:], 
    #                                       no_dmso_samples.loc[:,no_dmso_samples.columns != 'CS(C)=O'],
    #                                       drugTargets,
    #                                       i)
    merged_interactions = inferDrugTarget(interactions,
                                          torch.tensor(global_threshold_all),
                                          model.drugLayer.mask.T.detach(), 
                                          drugInput,
                                          drugTargets,
                                          i)
    merged_interactions['model_no'] = i
    if i == 0 :
        all_merged_interactions = merged_interactions.copy()
    else:
        all_merged_interactions = all_merged_interactions.append(merged_interactions)
        all_merged_interactions = all_merged_interactions.reset_index(drop=True)
    print2log('Finished model %s'%i)
all_merged_interactions.to_csv(ensembles_path +'all_merged_interactions_drugs.csv')