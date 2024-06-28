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
import argparse
import logging
import time
start_time = time.time()

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Drug-Targets for drugs of interest')
parser.add_argument('--ensembles_path', action='store', default="../results/A375_ensembles/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_model_modeltype4_model")
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
parser.add_argument('--drugInputFile', action='store',default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv')
parser.add_argument('--drugTargetsFile', action='store',default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv')
parser.add_argument('--TFOutFile', action='store',default='../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv')
parser.add_argument('--drugSimilarityFile', action='store',default='../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv')
parser.add_argument('--interactionsPath', action='store',default=None)
parser.add_argument('--Y_ALL_path', action='store',default=None)
parser.add_argument('--Y_ALL_masked_path', action='store',default=None)

args = parser.parse_args()
ensembles_path = args.ensembles_path
numberOfModels = int(args.numberOfModels)
ConvertToEmpProb = args.ConvertToEmpProb
if type(ConvertToEmpProb) == str :
    ConvertToEmpProb = eval(ConvertToEmpProb)
inputPattern = args.inputPattern
inputPath = ensembles_path + inputPattern
drugInputFile = args.drugInputFile
drugTargetsFile = args.drugTargetsFile
TFOutFile = args.TFOutFile
drugSimilarityFile = args.drugSimilarityFile
interactionsPath = args.interactionsPath
Y_ALL_masked_path = args.Y_ALL_masked_path
Y_ALL_path = args.Y_ALL_path

#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
spectralCapacity = np.exp(np.log(1e-2)/bionetParams['iterations'])
### Load the data
drugInput = pd.read_csv(drugInputFile, sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pd.read_csv(drugTargetsFile, sep='\t', low_memory=False, index_col=0)
TFOutput = pd.read_csv(TFOutFile, sep='\t', low_memory=False, index_col=0)
TFOutput=TFOutput.loc[drugInput.index,:]
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
drugSim = pd.read_csv(drugSimilarityFile,index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

#drugInput = drugInput[drugInput.loc[:,'CS(C)=O']==0]
dmso_ind = np.where(drugInput.columns.values=='CS(C)=O')[0][0]
#all_drugs = list(drugInput.columns.values)

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

model = torch.load(inputPath+str(0)+".pt")
all_scores = np.zeros((numberOfModels,model.drugLayer.mask.T.shape[0],model.drugLayer.mask.T.shape[1]))
mask = model.drugLayer.mask.T
### Define gradient thresholds
thresholds = list(np.logspace(-3.5, 3.5, num=50))
Y_ALL = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
Y_ALL_masked = torch.zeros(numberOfModels,len(thresholds),Y.shape[0],Y.shape[1])
filtered_all = torch.zeros(numberOfModels,len(thresholds))
for i in range(numberOfModels):
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    Y_ALL[i, :, :] = Yhat.detach()
    Xin = model.drugLayer(X)
    # Get starting activity not steady state
    Yin = model.inputLayer(Xin)

    print2log('Begin score calculation with linear algebra for model %s'%i)
    W = torch.mul(model.drugLayer.drugSim, model.drugLayer.Wsim)
    A = torch.mul(model.drugLayer.A, model.drugLayer.mask)
    kappa = model.drugLayer.bn.weight / torch.sqrt(model.drugLayer.bn.running_var + model.drugLayer.bn.eps)
    K = kappa
    scores = torch.mm(A, (kappa * W).T)
    scores = scores.T.detach()

    interactions = pd.DataFrame(scores.numpy())
    interactions.index = drugInput.columns
    interactions.columns = drugTargets.columns
    interactions.to_csv(interactionsPath +'interactionScores_'+str(i)+'.csv')
    #interactions = pd.read_csv(outPattern+'_interactionScores_'+str(i)+'.csv',index_col=0)

    # ## alternative
    # scores_2 = model.drugLayer(torch.eye(X.shape[1]).double())
    # interactions_2 = pd.DataFrame(scores_2.detach().numpy())
    # interactions_2.index = drugInput.columns
    # interactions_2.columns = drugTargets.columns
    # interactions_2.to_csv(interactionsPath + 'interactionScores_v2_' + str(i) + '.csv')

    scores = torch.abs(scores).detach()
    for j in range(len(thresholds)):
        th = thresholds[j]
        mask = torch.mm((1.0 * (X != 0)).double(), (1.0 * (scores >= th)).double())
        filtered_all[i, j] = torch.sum(1.0 * (scores >= th))
        Xin_masked = torch.mul(Xin, mask)
        fullX_masked = model.inputLayer(Xin_masked)
        YhatFull_masked = model.network(fullX_masked)
        Yhat_masked = model.projectionLayer(YhatFull_masked)
        Y_ALL_masked[i, j, :, :] = Yhat_masked.detach()

    torch.save(Y_ALL, Y_ALL_path)
    torch.save(Y_ALL_masked, Y_ALL_masked_path)
    all_scores[i, :, :] = interactions.values
    np.save('all_scores.npy',all_scores)
    one_mod_time = time.time()
    if (i==0):
        print2log('One model time : %s'%(one_mod_time - start_time))
torch.save(Y_ALL, Y_ALL_path)
torch.save(Y_ALL_masked, Y_ALL_masked_path)
# Calculate average score
scores = torch.tensor(np.mean(all_scores,0))
scores_std = torch.tensor(np.std(all_scores,0))
interactions = pd.DataFrame(scores.numpy())
interactions.index = drugInput.columns
interactions.columns = drugTargets.columns
interactions.to_csv(interactionsPath+'interactionScoresEnsembled.csv')
interactions_std = pd.DataFrame(scores_std.numpy())
interactions_std.index = drugInput.columns
interactions_std.columns = drugTargets.columns
interactions_std.to_csv(interactionsPath+'interactionScoresSTDEnsembled.csv')

end_time = time.time()
execution_time = end_time - start_time
print2log("Execution time: %s"%execution_time)