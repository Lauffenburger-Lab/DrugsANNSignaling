import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
from captum.attr import IntegratedGradients
from matplotlib import pyplot as plt
import seaborn as sns
import argparse
import logging
sns.set()

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Drug-Targets for drugs of interest')
parser.add_argument('--ensembles_path', action='store', default="../results/regularizationTune/")
parser.add_argument('--save_path', action='store', default="../results/regularizationTuneInteractionResults/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_modeltype4_lamda")
parser.add_argument('--cell_line', action='store', default="VCAP")
parser.add_argument('--numberOfModels', action='store', default=9)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
args = parser.parse_args()
ensembles_path = args.ensembles_path
save_path = args.save_path
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

#%% Data loaded  now begin analysis
### Define gradient thresholds
thresholds = list(np.logspace(-3.5, 3.5, num=50))

Y_ALL = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
Y_ALL_masked = torch.zeros(numberOfModels,len(thresholds),Y.shape[0],Y.shape[1])
filtered_all = torch.zeros(numberOfModels,len(thresholds))
drug_ratioMatrix = torch.arange(0, 1.01, 0.01)
drug_ratioMatrix = drug_ratioMatrix.T.repeat(drugInput.shape[0],drugInput.shape[1],1)
drug_ratioMatrix = drug_ratioMatrix.double()
for i in range(numberOfModels):
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    ig = IntegratedGradients(model.drugLayer)
    Y_ALL[i,:,:] = Yhat.detach()
    Xin = model.drugLayer(X)
    
    # Get starting activity not steady state
    Yin = model.inputLayer(Xin)
    print2log('Begin score calculation with integrated gradients for model %s'%i)
    # Per output latent variable input importance translation captum
    # 1st dimesion input
    # 2nd dimesion output
    scores = torch.zeros (model.drugLayer.mask.T.shape)
    for target in range(scores.shape[1]):
        attr = ig.attribute(torch.eye(X.shape[1]).double(),target=target,n_steps=100)
        scores[:,target] = torch.diagonal(attr, 0)

    interactions = pd.DataFrame(scores.numpy())
    interactions.index = drugInput.columns
    interactions.columns = drugTargets.columns
    interactions.to_csv(save_path+'error_based_inference/'+cell+'_interactionScores_lamda'+str(i)+'.csv')
    scores = torch.abs(scores).detach()
    
    print2log('Begin calculating output using masking')
    for j in range(len(thresholds)):
        th = thresholds[j]
        mask = torch.mm((1.0 * (X != 0)).double(), (1.0*(scores>=th)).double())
        filtered_all[i,j] = torch.sum(1.0*(scores>=th))
        Xin_masked = torch.mul(Xin, mask)
        fullX_masked = model.inputLayer(Xin_masked)
        YhatFull_masked = model.network(fullX_masked)
        Yhat_masked = model.projectionLayer(YhatFull_masked)
        Y_ALL_masked[i, j, :,:] = Yhat_masked.detach()

    print2log('Finished model %s'%i)


torch.save(Y_ALL,ensembles_path+'Y_'+cell+'_ALL.pt')
torch.save(Y_ALL_masked,ensembles_path+'Y_'+cell+'_ALL_MASKED.pt')
torch.save(filtered_all,ensembles_path+'filtered_'+cell+'_ALL_MASKED.pt')