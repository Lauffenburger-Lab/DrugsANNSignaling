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
import seaborn as sns
sns.set()
import argparse
import logging

def pearson_r(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = torch.mean(x, dim=0)
    my = torch.mean(y, dim=0)
    xm, ym = x - mx, y - my
    r_num = torch.sum(xm * ym,dim=0)
    x_square_sum = torch.sum(xm * xm,dim=0)
    y_square_sum = torch.sum(ym * ym,dim=0)
    r_den = torch.sqrt(x_square_sum * y_square_sum)
    r = r_num / r_den
    return r

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Infer off-target effect')
parser.add_argument('--ensembles_path', action='store', default="../results/FinalEnsemble/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_modeltype4_model")
parser.add_argument('--cell_line', action='store', default="VCAP")
parser.add_argument('--numberOfModels', action='store', default=51)
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

### Load network
#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for a375 was 120
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

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)
#Xin = torch.mm(X,drugLayer)
Υ_ALL = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
Υ_ALL_masked_1 = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
Υ_ALL_masked_2 = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
for i in range(numberOfModels):
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    Υ_ALL[i,:,:] = Yhat.detach()
    mask = torch.mm(1.0*(X!=0).double(),model.drugLayer.mask.T)

    Xin = model.drugLayer(X)

    Xin_masked = torch.mul(Xin, mask)
    fullX_masked = model.inputLayer(Xin_masked)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_1 = model.projectionLayer(YhatFull_masked)
    Υ_ALL_masked_1[i, :, :] = Yhat_masked_1.detach()

    Xin_masked = torch.mul(Xin, 1.0 - mask)
    fullX_masked = model.inputLayer(Xin_masked)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_2 = model.projectionLayer(YhatFull_masked)
    Υ_ALL_masked_2[i, :, :] = Yhat_masked_2.detach()

# Calculate mean predictions
Yhat = torch.mean(Υ_ALL,0)
Yhat_masked = torch.mean(Υ_ALL_masked_1,0)
Yhat_masked_2 = torch.mean(Υ_ALL_masked_2,0)

# Per TF performance in VCAP
performance = pearson_r(Y.detach(), Yhat.detach()).detach().numpy()
performance = pd.DataFrame(performance)
performance.index = TFOutput.columns
performance.columns = ['r']
performance.to_csv(ensembles_path+cell+'TrainEnsemblePerformance.csv')

# Calculate Delta1
Delta1 = Yhat_masked - Yhat
Delta1 = pd.DataFrame(Delta1.detach().numpy())
Delta1.columns =  TFOutput.columns
Delta1.index =  TFOutput.index
Delta1.to_csv(ensembles_path+'DeltaTF1.csv')

# Calculate Delta2
Delta2 = Yhat_masked_2 - Yhat
Delta2 = pd.DataFrame(Delta2.detach().numpy())
Delta2.columns =  TFOutput.columns
Delta2.index =  TFOutput.index
Delta2.to_csv(ensembles_path+ 'DeltaTF2.csv')