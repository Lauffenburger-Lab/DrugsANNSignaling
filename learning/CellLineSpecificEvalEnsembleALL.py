import torch
import numpy
from scipy.stats import beta
import matplotlib.pyplot as plt
import bionetworkWithDrugs as bionetwork
import pandas
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Cell-line-specific training evaluation')
parser.add_argument('--ensembles_path', action='store', required=True,help='Path to the ensembles folder')
parser.add_argument('--inputPattern', action='store', required=True,help='Input file pattern for trained models')
parser.add_argument('--numberOfModels', action='store', default=50,help='Number of trained models in the ensemble')
parser.add_argument('--ConvertToEmpProb', action='store',default=False,help='Should we apply sigmoid-like transformation the TF activity file? (default=False, because it has already applied)')
parser.add_argument('--drugInputFile', action='store',required=True,help='Model`s Input: Drug concetrations file')
parser.add_argument('--drugTargetsFile', action='store',required=True,help='File containing drug-target interactions used to train the model')
parser.add_argument('--TFOutFile', action='store',required=True,help='Model`s Outuput: TF activity file')
parser.add_argument('--drugSimilarityFile', action='store',required=True,help='Pre-calculated drug similarity matrix used to train the model')
parser.add_argument('--PKN', action='store', required=True,help='PKN file')
parser.add_argument('--PknAnnotation', action='store', required=True,help='PKN annotation file')
parser.add_argument('--res_dir', action='store', required=True,help='Results directory to save results')
parser.add_argument('--CellPrefix', action='store', required=True,help='Prefix for saving the files (the cell line name e.g. A375)')
args = parser.parse_args()
ConvertToEmpProb = args.ConvertToEmpProb
if type(ConvertToEmpProb) == str :
    ConvertToEmpProb = eval(ConvertToEmpProb)
no_models = args.numberOfModels
ensembles_path = args.ensembles_path
inputPattern = args.inputPattern
inputPath = ensembles_path + inputPattern
drugInputFile = args.drugInputFile
drugTargetsFile = args.drugTargetsFile
TFOutFile = args.TFOutFile
drugSimilarityFile = args.drugSimilarityFile
PknAnnotation = args.PknAnnotation
PKN = args.PKN
res_dir = args.res_dir
CellPrefix = args.CellPrefix

def Rsquared(y_pred,y_true):
    SS_res = torch.sum((y_true - y_pred) ** 2,0)
    SS_tot = torch.sum((y_true - torch.mean(y_true,0).repeat(y_true.shape[0],1)) ** 2,0)
    return torch.mean(1 - SS_res / SS_tot)

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
    return r #torch.mean(r)

def lineOfIdentity():
    xBounds = plt.xlim()
    yBounds = plt.ylim()
    minLevel = numpy.min([xBounds[0], yBounds[0]])
    maxLevel = numpy.max([xBounds[1], yBounds[1]])
    plt.plot([minLevel, maxLevel], [minLevel, maxLevel], 'black', label='_nolegend_')


#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork(PKN)
annotation = pandas.read_csv(PknAnnotation, sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
spectralCapacity = numpy.exp(numpy.log(1e-2)/bionetParams['iterations'])

### Load all data
drugInput = pandas.read_csv(drugInputFile, sep='\t', low_memory=False,index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pandas.read_csv(drugTargetsFile, sep='\t', low_memory=False, index_col=0)
TFOutput = pandas.read_csv(TFOutFile, sep='\t', low_memory=False, index_col=0)
TFOutput = TFOutput.loc[drugInput.index, :]
train_pear = numpy.zeros((no_models,TFOutput.shape[1]))
Y_ALL_TRAIN = torch.zeros((no_models,TFOutput.shape[0],TFOutput.shape[1]))

#Subset input and output to intersecting modeldes
druginName = drugInput.columns.values
inName = drugTargets.columns.values
outconds = drugInput.index.values
outName = TFOutput.columns.values
inName = numpy.intersect1d(nodeNames, inName)
outName = numpy.intersect1d(nodeNames, outName)
inNameGene = [uniprot2gene[x] for x in inName]
outNameGene = [uniprot2gene[x] for x in outName]
TFOutput = TFOutput.loc[outconds,outName]
drugTargets = drugTargets.loc[:,inName]
if ConvertToEmpProb==True:
    print2log('Convereted to Emprirical probabilities')
    TFOutput = 1/(1+numpy.exp(-TFOutput))

drugTargets = drugTargets.loc[drugInput.columns.values,:]
drugSim = pandas.read_csv(drugSimilarityFile,index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)
    
for i in range(no_models):
    # Here we load the pre-trained drug layer, which was trained on cell-line data
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()

    Yhat, YhatFull = model(X)
    Y_ALL_TRAIN[i, :, :] = Yhat
    
    train_pear[i,:] = pearson_r(Y.detach(), Yhat.detach()).detach().numpy()
    print2log('Finished model %s'%i)

train_pear = pandas.DataFrame(train_pear)
train_pear.columns = TFOutput.columns
train_pear.index = ['model_ '+str(ind) for ind in range(no_models)]
train_pear.to_csv(res_dir+CellPrefix+'trainPerformance_perTF.csv')

torch.save(Y_ALL_TRAIN,res_dir+'Y_'+CellPrefix+'_ensemble_train.pt')

Y_MEAN_TRAIN = torch.mean(Y_ALL_TRAIN, 0)
Y_STD_TRAIN = torch.std(Y_ALL_TRAIN, 0)

pearTrain = pearson_r(Y_MEAN_TRAIN.detach(), Y.detach()).detach().numpy()
dist = beta(Y_MEAN_TRAIN.shape[0] / 2 - 1, Y_MEAN_TRAIN.shape[0] / 2 - 1, loc=-1, scale=2)
pvaluesTrain = [2 * dist.cdf(-abs(r)) for r in list(pearTrain)]
pearTrain = pandas.DataFrame(pearTrain)
pearTrain.index = TFOutput.columns
pearTrain.to_csv(res_dir+CellPrefix+'_trainEnsemblePerformance.csv')
pvaluesTrain = pandas.DataFrame(numpy.array(pvaluesTrain))
pvaluesTrain.index = TFOutput.columns
pvaluesTrain.to_csv(res_dir+CellPrefix+'_trainEnsemblePvalues.csv')

Y_MEAN_TRAIN = pandas.DataFrame(Y_MEAN_TRAIN.detach().numpy())
Y_MEAN_TRAIN.index = TFOutput.index
Y_MEAN_TRAIN.columns = TFOutput.columns
Y_MEAN_TRAIN.to_csv(res_dir+'Y_'+CellPrefix+'_mean_train.csv')
Y_STD_TRAIN = pandas.DataFrame(Y_STD_TRAIN.detach().numpy())
Y_STD_TRAIN.index = TFOutput.index
Y_STD_TRAIN.columns = TFOutput.columns
Y_STD_TRAIN.to_csv(res_dir+'Y_'+CellPrefix+'_std_train.csv')