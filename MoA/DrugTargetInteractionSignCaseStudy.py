import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
from matplotlib import pyplot as plt
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

## This function resets the gradients of a model
def resetGradients(model):
    for var in model.parameters():
        var.grad = None

## This function infers the binary drug-target interactions based on pre-calculated gradient score thresholds for each drug for each model in the ensemble
## WITHOUT SAVING THESE RESULTS
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

## All the following aguments have defaults that correspond to the 
## case study of extracting a mechanism of action of the off-target effect of
## Lestaurtinib on FOXM1 through the CDK2 signaling node in A375 cancer cell lines
parser = argparse.ArgumentParser(prog='Infer sign of interaction')
parser.add_argument('--inputPattern', action='store', default='l1000_latest_model_modeltype4_a375_case_study')
parser.add_argument('--ensembles_path', action='store', default='../../results/case_study/')
parser.add_argument('--DrugsIn', action='store', default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_A375-conditions_drugs.tsv')
parser.add_argument('--TargetsIn', action='store', default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_A375-drugs_targets.tsv')
parser.add_argument('--TFsOut', action='store', default='../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv')
parser.add_argument('--ChemicalSims', action='store', default='../preprocessing/preprocessed_data/ChemicalSims/lvl3_similarities_A375.csv')
parser.add_argument('--PKN', action='store', default='../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
parser.add_argument('--PknAnnotation', action='store', default='../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv')
parser.add_argument('--res_dir', action='store', default='../../results/case_study/')
parser.add_argument('--interactionScorePattern', action='store', default='interactionScores')
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
## Define case study parameters
parser.add_argument('--moa_off_target', action='store',default='any')
parser.add_argument('--cell_line', action='store',default='A375')
parser.add_argument('--node', action='store',default='P24941')
parser.add_argument('--node_gene', action='store',default='CDK2')
parser.add_argument('--drug', action='store',default='C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13')
parser.add_argument('--drug_name', action='store',default='lestaurtinib')
parser.add_argument('--sample', action='store',default='CPC014_A375_6H:BRD-K23192422-001-01-1:10')
args = parser.parse_args()
numberOfModels = int(args.numberOfModels)
ConvertToEmpProb = args.ConvertToEmpProb
if type(ConvertToEmpProb) == str :
    ConvertToEmpProb = eval(ConvertToEmpProb)
ensembles_path = args.ensembles_path
inputPattern = args.inputPattern
DrugsIn = args.DrugsIn
TargetsIn = args.TargetsIn
TFsOut = args.TFsOut
ChemicalSims = args.ChemicalSims
PknAnnotation = args.PknAnnotation
PKN = args.PKN
res_dir = args.res_dir
interactionScorePattern = args.interactionScorePattern
inputPath = ensembles_path + 'models/' + inputPattern
cell = args.cell_line

### Choose node and drug to investigate
node = args.node
node_gene = args.node_gene
drug = args.drug
drug_name = args.drug_name
sample = args.sample
moa_off_target = args.moa_off_target
# CDK1 = P06493
# CDK2 = P24941
# CDK6 = Q00534

### Load network
#Load network
#### BE CAREFUL!!!!
#### In networkList the sources and targets are flipped !!!
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork(PKN)
annotation = pd.read_csv(PknAnnotation, sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for a375 was 120
spectralCapacity = np.exp(np.log(1e-2)/bionetParams['iterations'])
### Load the data
drugInput = pd.read_csv(DrugsIn, sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pd.read_csv(TargetsIn, sep='\t', low_memory=False, index_col=0)
TFOutput = pd.read_csv(TFsOut, sep='\t', low_memory=False, index_col=0)
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
drugSim = pd.read_csv(ChemicalSims,index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

# Keep only drug/sample of interest
sample_ind = np.where(TFOutput.index==sample)[0]
drugInput = pd.DataFrame(drugInput.loc[sample,:]).T
TFOutput = pd.DataFrame(TFOutput.loc[sample,:]).T

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

### Get mapping of network nodes to TFs
dictionary = dict(zip(nodeNames, list(range(len(nodeNames)))))
node_index_inNet = dictionary[node]

# Get pre-calculated global grad score thresholds for each drug for each model in the ensemble, that will be used for inferring drug-target interactions
global_grad_scores = pd.read_csv(ensembles_path+"all_drugs_global_thresholds.csv",index_col=0)
global_grad_scores = torch.tensor(global_grad_scores.loc[drug,:].values)

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
    
    ### Get drug-target interactions
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/'+interactionScorePattern+'_%s.csv' % i,index_col=0)
    merged_interactions = inferDrugTarget(interactions,global_grad_scores,
                                          model.drugLayer.mask.T.detach(), drugInput,drugTargets,
                                          i)
    drug_target_nodes = merged_interactions[merged_interactions['drug'] == drug]
    drug_target_nodes = drug_target_nodes[drug_target_nodes['Inferred'] == 'Interaction']
    drug_target_nodes = drug_target_nodes.variable.unique()
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/'+interactionScorePattern+'_%s.csv' % i,index_col=0)
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
plt.show()
# print2log('Mean score: %s'%np.mean(score))
# print2log('Median score: %s'%np.median(score))
print2log('Positive scores: %s'%(np.sum(np.array(score)>0)/len(score)))
plt.figure()
plt.hist(node_activity)
plt.show()
# print2log('Mean node_activity: %s'%np.mean(node_activity))
# print2log('Median node_activity: %s'%np.median(node_activity))
print2log('Greater than 0.3 node_activity: %s'%(np.sum(np.array(node_activity)>0.3)/len(node_activity)))
print2log('Greater than 0.5 node_activity: %s'%(np.sum(np.array(node_activity)>0.5)/len(node_activity)))
plt.figure()
plt.hist(node_start_signal,40)
plt.show()
print2log('Node start signal positive: %s'%(np.sum(np.array(node_start_signal)>0)/len(node_start_signal)))
