import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
import seaborn as sns
import argparse
import logging
sns.set()

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

dmso_global_errors = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),len(thresholds)))
dmso_unmasked_global_errors = np.zeros((len(drugInput[drugInput.loc[:,'CS(C)=O']!=0]),1))
all_global_errors = np.zeros((drugInput.shape[1],len(thresholds)))
all_unmasked_global_errors = np.zeros((drugInput.shape[1],1))

Yhat = torch.mean(Y_ALL,0)
Yhat_masked = torch.mean(Y_ALL_masked,0)

dmso_counter = 0
ind = 0

for j in range(X.shape[0]):
    if drugInput.iloc[j,dmso_ind] == 0:
        all_drugs.append(drugInput.columns.values[np.where(drugInput.iloc[j,:]!=0)[0]][0])
        global_error = torch.mean(torch.abs(Yhat_masked[:,j,:]-Y[j,:]).squeeze(),1).detach().squeeze().numpy()
        unmasked_global_error = torch.mean(torch.abs(Yhat[j, :] - Y[j, :])).detach().squeeze().numpy()
        all_global_errors[ind,:] = global_error
        all_unmasked_global_errors[ind,:] = unmasked_global_error
        ind = ind +1
    else:
        global_error = torch.mean(torch.abs(Yhat_masked[:,j,:]-Y[j,:]).squeeze(),1).detach().squeeze().numpy()
        unmasked_global_error = torch.mean(torch.abs(Yhat[j, :] - Y[j, :])).detach().squeeze().numpy()
        dmso_global_errors[dmso_counter,:] = global_error
        dmso_unmasked_global_errors[dmso_counter,:] = unmasked_global_error
        if dmso_counter==0:
            final_dmso_ind = ind
            all_global_errors[ind,:] = 0.
            all_unmasked_global_errors[ind,:] = 0.
            all_drugs.append(drugInput.columns.values[np.where(drugInput.iloc[j,:]!=0)[0]][0])
            ind = ind+1
        dmso_counter = dmso_counter + 1   
    print2log('Finished drug %s'%j)

all_global_errors[final_dmso_ind,:] = np.mean(dmso_global_errors,0)
all_global_errors_df = pd.DataFrame(all_global_errors)
all_global_errors_df.index = all_drugs
all_global_errors_df = all_global_errors_df.T
all_global_errors_df["grad_threshold"] = thresholds
all_global_errors_df.to_csv(ensembles_path +'all_drugs_global_errors.csv')

all_unmasked_global_errors[final_dmso_ind,:] = np.mean(dmso_unmasked_global_errors,0)
all_unmasked_global_errors = pd.DataFrame(all_unmasked_global_errors)
all_unmasked_global_errors.index = all_drugs
all_unmasked_global_errors.to_csv(ensembles_path +'all_drugs_unmasked_global_errors.csv')
