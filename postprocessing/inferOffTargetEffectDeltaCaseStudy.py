import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
import argparse
import logging
import time
start_time = time.time()

# Function to calculate pearson correlation between predicted and true activities per TF in pytorch
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
parser.add_argument('--inputPathPattern', action='store',required=True,help='Path to the ensembles folder')
parser.add_argument('--DrugsIn', action='store',required=True,help='Model`s Input: Drug concetrations file')
parser.add_argument('--TargetsIn', action='store',required=True,help='File containing drug-target interactions used to train the model')
parser.add_argument('--TFsOut', action='store',required=True,help='Model`s Outuput: TF activity file')
parser.add_argument('--ChemicalSims', action='store',required=True,help='Pre-calculated drug similarity matrix used to train the model')
parser.add_argument('--PKN', action='store',required=True,help='Path to the trimmed PKN file')
parser.add_argument('--PknAnnotation', action='store',required=True,help='Path to the trimmed PKN annotation file')
parser.add_argument('--res_dir', action='store',required=True,help='Path to the results folder')
parser.add_argument('--numberOfModels', action='store',required=True,help='Number of trained models in the ensemble')
parser.add_argument('--ConvertToEmpProb', action='store',default=False,help='Should we apply sigmoid-like transformation the TF activity file? (default=False, because it has already applied)')
args = parser.parse_args()
numberOfModels = int(args.numberOfModels)
ConvertToEmpProb = args.ConvertToEmpProb
if type(ConvertToEmpProb) == str :
    ConvertToEmpProb = eval(ConvertToEmpProb)
inputPathPattern = args.inputPathPattern
DrugsIn = args.DrugsIn
TargetsIn = args.TargetsIn
TFsOut = args.TFsOut
ChemicalSims = args.ChemicalSims
PknAnnotation = args.PknAnnotation
PKN = args.PKN
res_dir = args.res_dir

### Load network
#Load network
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

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)
Y_ALL = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
Y_ALL_masked_1 = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
Y_ALL_masked_2 = torch.zeros(numberOfModels,Y.shape[0],Y.shape[1])
model_times = []
for i in range(numberOfModels):
    prev_time = time.time()
    model = torch.load(inputPathPattern+str(i)+".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    Y_ALL[i,:,:] = Yhat.detach()
    mask = torch.mm(1.0*(X!=0).double(),model.drugLayer.mask.T)

    Xin = model.drugLayer(X)

    Xin_masked = torch.mul(Xin, mask)
    fullX_masked = model.inputLayer(Xin_masked)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_1 = model.projectionLayer(YhatFull_masked)
    Y_ALL_masked_1[i, :, :] = Yhat_masked_1.detach()
    print2log('Finished model %s'%i)
    model_times.append(time.time()-prev_time)

print2log('Average time per model = %s seconds'%np.mean(model_times))
# Calculate mean predictions
Yhat = torch.mean(Y_ALL,0)
Yhat_masked = torch.mean(Y_ALL_masked_1,0)

# Per TF performance in the cell line of interest
performance = pearson_r(Y.detach(), Yhat.detach()).detach().numpy()
performance = pd.DataFrame(performance)
performance.index = TFOutput.columns
performance.columns = ['r']
performance.to_csv(res_dir+'TrainEnsemblePerformance.csv')

# Calculate Delta1
Delta1 = Yhat_masked - Yhat
Delta1 = pd.DataFrame(Delta1.detach().numpy())
Delta1.columns =  TFOutput.columns
Delta1.index =  TFOutput.index
Delta1.to_csv(res_dir+'DeltaTF1.csv')

runtime = time.time() - start_time
print2log('Total process time = %s seconds'%runtime)