import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
import logging
sns.set()

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Drug-Targets for drugs of interest')
parser.add_argument('--ensembles_path', action='store', default="CVL1000_Paper/A375_ensembles/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_model_modeltype4_model")
parser.add_argument('--cell_line', action='store', default="A375")
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
parser.add_argument('--TF_of_interest', action='store',default="FOXM1")
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
TF_gene =  args.TF_of_interest

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
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('data/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('data/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
TF = annotation.code.values[np.where(annotation.name==TF_gene)[0]][0]
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
spectralCapacity = np.exp(np.log(1e-2)/bionetParams['iterations'])
### Load the data
drugInput = pd.read_csv('data/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pd.read_csv('data/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput = pd.read_csv('data/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
cellInput = pd.read_csv('data/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
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
drugSim = pd.read_csv('../out_lvl3_similaritiess.csv',index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

#drugInput = drugInput[drugInput.loc[:,'CS(C)=O']==0]
dmso_ind = np.where(drugInput.columns.values=='CS(C)=O')[0][0]
#all_drugs = list(drugInput.columns.values)
TF_ind = np.where(TFOutput.columns==TF)[0]

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

### Define gradient thresholds
thresholds = list(np.logspace(-3.5, 3.5, num=50))
# thresholds = list(np.logspace(-5, 5, num=50))
thresh = 0.25
Y_ALL = torch.load(ensembles_path+'preds/Y_'+cell+'_ALL.pt')
Y_ALL_masked = torch.load(ensembles_path+'preds/Y_'+cell+'_ALL_MASKED.pt')

# Yhat = torch.mean(Y_ALL,0)
# Yhat_masked = torch.mean(Y_ALL_masked,0)

### TF errors initialization
dmso_TF_errors = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),numberOfModels,len(thresholds)))
dmso_unmasked_TF_errors = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),numberOfModels,1))
all_TF_errors = np.zeros((drugInput.shape[1],numberOfModels,len(thresholds)))
all_unmasked_TF_errors = np.zeros((drugInput.shape[1],numberOfModels,1))

### Global errors initialization
dmso_global_errors = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),numberOfModels,len(thresholds)))
dmso_unmasked_global_errors = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),numberOfModels,1))
all_global_errors = np.zeros((drugInput.shape[1],numberOfModels,len(thresholds)))
all_unmasked_global_errors = np.zeros((drugInput.shape[1],numberOfModels,1))

all_drugs = []
ind = 0

print2log('Begin caclulating thresholds for each drug')
dmso_counter = 0

for j in range(X.shape[0]):
    if drugInput.iloc[j,dmso_ind] == 0:
        all_drugs.append(drugInput.columns.values[np.where(drugInput.iloc[j,:]!=0)[0]][0])
        # Calculate for specific TF
        TF_error = torch.abs(Y_ALL_masked[:,:,j,TF_ind]-Y[j,TF_ind]).detach().squeeze().numpy()
        unmasked_TF_error = torch.abs(Y_ALL[:,j, TF_ind] - Y[j, TF_ind]).detach().squeeze().numpy()
        all_TF_errors[ind,:,:] = TF_error
        all_unmasked_TF_errors[ind,:,0] = unmasked_TF_error
        # Calculate for globally
        global_error = torch.mean(torch.abs(Y_ALL_masked[:,:,j,:]-Y[j,:]),2).detach().squeeze().numpy()
        unmasked_global_error = torch.mean(torch.abs(Y_ALL[:,j, :] - Y[j, :]),1).detach().squeeze().numpy()
        all_global_errors[ind,:,:] = global_error
        all_unmasked_global_errors[ind,:,0] = unmasked_global_error
        ind = ind +1
    else:
        # Calculate for specific TF
        TF_error = torch.abs(Y_ALL_masked[:,:,j,TF_ind]-Y[j,TF_ind]).detach().squeeze().numpy()
        unmasked_TF_error = torch.abs(Y_ALL[:,j, TF_ind] - Y[j, TF_ind]).detach().squeeze().numpy()
        dmso_TF_errors[dmso_counter,:,:] = TF_error
        dmso_unmasked_TF_errors[dmso_counter,:,0] = unmasked_TF_error
        # Calculate for globally
        global_error = torch.mean(torch.abs(Y_ALL_masked[:,:,j,:]-Y[j,:]),2).detach().squeeze().numpy()
        unmasked_global_error = torch.mean(torch.abs(Y_ALL[:,j, :] - Y[j, :]),1).detach().squeeze().numpy()
        dmso_global_errors[dmso_counter,:,:] = global_error
        dmso_unmasked_global_errors[dmso_counter,:,0] = unmasked_global_error
        if dmso_counter==0:
            final_dmso_ind = ind
            all_TF_errors[ind,:,:] = 0.
            all_unmasked_TF_errors[ind,:,:] = 0.
            all_drugs.append(drugInput.columns.values[np.where(drugInput.iloc[j,:]!=0)[0]][0])
            ind = ind+1
        dmso_counter = dmso_counter + 1   
    print2log('Finished drug %s'%j)

# save for specific TF
all_TF_errors[final_dmso_ind,:,:] = np.mean(dmso_TF_errors,0)
TF_error_df = pd.DataFrame([list(l) for l in np.rollaxis(all_TF_errors,1,0)]).stack().apply(pd.Series).reset_index(1, drop=True)
TF_error_df.columns = thresholds
TF_error_df['drug'] = all_drugs*numberOfModels
TF_error_df.to_csv(ensembles_path +'allmodels_all_drugs_'+TF_gene+'_errors.csv')
all_unmasked_TF_errors[final_dmso_ind,:,:] = np.mean(dmso_unmasked_TF_errors,0)
all_unmasked_TF_errors = all_unmasked_TF_errors.squeeze()
all_unmasked_TF_errors_df = pd.DataFrame(all_unmasked_TF_errors)
all_unmasked_TF_errors_df.index = all_drugs
all_unmasked_TF_errors_df.to_csv(ensembles_path +'allmodels_all_drugs_unmasked_'+TF_gene+'_errors.csv')

# save globally for specific TF
all_global_errors[final_dmso_ind,:,:] = np.mean(dmso_global_errors,0)
global_error_df = pd.DataFrame([list(l) for l in np.rollaxis(all_global_errors,1,0)]).stack().apply(pd.Series).reset_index(1, drop=True)
global_error_df.columns = thresholds
global_error_df['drug'] = all_drugs*numberOfModels
global_error_df.to_csv(ensembles_path +'allmodels_all_drugs_global_errors.csv')
all_unmasked_global_errors[final_dmso_ind,:,:] = np.mean(dmso_unmasked_global_errors,0)
all_unmasked_global_errors = all_unmasked_global_errors.squeeze()
all_unmasked_global_errors_df = pd.DataFrame(all_unmasked_global_errors)
all_unmasked_global_errors_df.index = all_drugs
all_unmasked_global_errors_df.to_csv(ensembles_path +'allmodels_all_drugs_unmasked_global_errors.csv')
