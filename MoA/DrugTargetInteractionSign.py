import os
import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
import torch.nn.functional as F
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit import Chem
import scipy
from sklearn.cluster import KMeans
from captum.attr import IntegratedGradients
from captum.attr import LayerConductance
from captum.attr import NeuronConductance
from matplotlib import pyplot as plt
import networkx as nx
import seaborn as sns
import argparse
import logging
from pathlib import Path
from operator import itemgetter
sns.set()

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

def resetGradients(model):
    for var in model.parameters():
        var.grad = None


def inferDrugTarget(interactions, global_thresholds, prior_mask, drugInput, drugTargets, model_no,
                    grad_thresholds=list(np.logspace(-3.5, 3.5, num=45)), thresh=0.25):
    global_thresholds = global_thresholds.detach().numpy()
    scores = torch.abs(torch.tensor(interactions.values))
    grad_thresh = global_thresholds[model_no]
    # Get new mask with drug-targets interactions
    mask = 1.0 * (scores >= grad_thresh)
    mask = mask.detach().numpy()

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

parser = argparse.ArgumentParser(prog='Infer MoA')
parser.add_argument('--ensembles_path', action='store', default="CVL1000_Paper/A549_ensembles/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_model_modeltype4_model")
parser.add_argument('--cell_line', action='store', default="A549")
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

### Choose TF and drug to investigate
# interestingSamples = pd.read_csv(ensembles_path+'interestingSamples.csv',index_col=1)
node = "P10826"
node_gene = "RARB"
drug = "Nc1ccccc1NC(=O)c1ccc(CNc2nccc(n2)-c2cccnc2)cc1"
drug_name = "mocetinostat"
sample = "CPC014_A549_6H:BRD-K16485616-001-03-0:10"
#%%
### Get drug-target interactions
# merged_interaction = pd.read_csv(ensembles_path+'InteractionScores/l1000_modeltype4_lamda6_'+cell+'_mergedInteractions_ROC.csv',index_col=0)
# merged_interaction = merged_interaction[merged_interaction['drug']==drug]
# merged_interaction = merged_interaction[merged_interaction['Inferred']=='Interaction']

### Load network
#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('data/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('data/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for a375 was 120
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

# Χ_all = torch.tensor(drugInput.values.copy(), dtype=torch.double)
# Y_all = torch.tensor(TFOutput.values, dtype=torch.double)
# Keep only drug/sample of interest
sample_ind = np.where(TFOutput.index==sample)[0]
drugInput = pd.DataFrame(drugInput.loc[sample,:]).T
TFOutput = pd.DataFrame(TFOutput.loc[sample,:]).T

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

### Get mapping of network nodes to TFs
dictionary = dict(zip(nodeNames, list(range(len(nodeNames)))))
node_index_inNet = dictionary[node]

# Get global grad score
# global_grad_scores = torch.load(ensembles_path+"all_drugs_global_thresholds.pt")
global_grad_scores = pd.read_csv(ensembles_path+"all_drugs_global_thresholds.csv",index_col=0)
global_grad_scores = torch.tensor(global_grad_scores.loc[drug,:].values)
# global_grad_scores = torch.load(ensembles_path+"InteractionScores/global_gradient_scores_"+drug+"_sample_ind"+str(sample_ind[0])+"_all_models.pt")

### Mapping of names and ids in graph
map = dict(zip(nodeNames, list(range(len(nodeNames)))))
df_all = pd.DataFrame({'source':[],'target':[],'weight':[],'name_source':[],'name_target':[],'interaction':[],'model_no':[]})
all_sources = pd.DataFrame({'sources':[],'model':[]})
df_all_source_types = pd.DataFrame({'type':[],'name':[]})
df_all_target_types = pd.DataFrame({'type':[],'name':[]})
score = []
node_activity = []
node_start_signal= []
drug_ratioMatrix = torch.arange(0, 1.01, 0.01)
drug_ratioMatrix = drug_ratioMatrix.T.repeat(drugInput.shape[0],drugInput.shape[1],1)
drug_ratioMatrix = drug_ratioMatrix.double()
for i in range(numberOfModels):#range(1):
    model = torch.load(inputPath+str(i)+".pt")
    resetGradients(model)
    model.eval()
    X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
    
    # # Get range of input activity for different concetrations of drugs
    # X_binned =  drug_ratioMatrix * X.unsqueeze(2)
    # X_binned = X_binned.squeeze()
    # X_binned = X_binned.T
    # Xin_binned =  model.drugLayer(X_binned)
    # Xin_range = torch.abs(torch.max(Xin_binned,0)[0] - torch.min(Xin_binned,0)[0]).unsqueeze(0)
    
    ### Get drug-target interactions
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/l1000_modeltype4_lamda6_' + cell + '_interactionScores_%s.csv' % i,index_col=0)
    # drugs = interactions.index.values
    # target_space = interactions.columns.values
    # interactions = torch.tensor(interactions.values) * Xin_range
    # interactions = pd.DataFrame(interactions.detach().numpy())
    # interactions.index =drugs
    # interactions = interactions.T
    # interactions.index = target_space
    # interactions = interactions.T
    merged_interactions = inferDrugTarget(interactions,global_grad_scores,
                                          model.drugLayer.mask.T.detach(), drugInput,drugTargets,
                                          i)
    drug_target_nodes = merged_interactions[merged_interactions['drug'] == drug]
    drug_target_nodes = drug_target_nodes[drug_target_nodes['Inferred'] == 'Interaction']
    drug_target_nodes = drug_target_nodes.variable.unique()
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/l1000_modeltype4_lamda6_' + cell + '_interactionScores_%s.csv' % i,index_col=0)
    if node in drug_target_nodes:
        score.append(interactions.loc[drug,node])
        target_ind = np.where(interactions.columns.values==node)[0]
        drug_ind = np.where(interactions.index.values==drug)[0]
        _ , YhatFull = model(X)
        Xin = model.drugLayer(X)
        node_activity.append(YhatFull[:,node_index_inNet].item())
        node_start_signal.append(Xin[:,target_ind].item())
    #print2log('Interaction score %s'%interactions.loc[drug,node])

#%%
### See the sign of the scores when the target has been inferred as interacting
plt.figure()
plt.hist(score)
# print2log('Mean score: %s'%np.mean(score))
# print2log('Median score: %s'%np.median(score))
print2log('Positive scores: %s'%(np.sum(np.array(score)>0)/len(score)))
plt.figure()
plt.hist(node_activity)
# print2log('Mean node_activity: %s'%np.mean(node_activity))
# print2log('Median node_activity: %s'%np.median(node_activity))
print2log('Greater than 0.3 node_activity: %s'%(np.sum(np.array(node_activity)>0.3)/len(node_activity)))
print2log('Greater than 0.5 node_activity: %s'%(np.sum(np.array(node_activity)>0.5)/len(node_activity)))
plt.figure()
# plt.hist(node_start_signal,40)
print2log('Node start signal positive: %s'%(np.sum(np.array(node_start_signal)>0)/len(node_start_signal)))