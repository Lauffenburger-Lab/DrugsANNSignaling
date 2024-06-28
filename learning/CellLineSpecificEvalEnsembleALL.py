import torch
import numpy
from scipy.stats import pearsonr,beta
import matplotlib.pyplot as plt
import bionetwork as SingleBionet
import bionetworkWithDrugs as bionetwork
#import bionetwork
import plotting
import pandas
import saveSimulations
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit import Chem
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

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


inputAmplitude = 3
projectionAmplitude = 1.2
ConvertToEmpProb = False

#Setup optimizer
batchSize = 25
MoAFactor = 0.1
spectralFactor = 1e-3
maxIter = 5000
noiseLevel = 10
L2 = 1e-6

no_models = 50


#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pandas.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
spectralCapacity = numpy.exp(numpy.log(1e-2)/bionetParams['iterations'])

### Load all data
drugInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_A375-conditions_drugs.tsv', 
                            sep='\t', low_memory=False,index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets_A375.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput = pandas.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
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
drugSim = pandas.read_csv('../preprocessing/preprocessed_data/ChemicalSims/lvl3_similarities_A375.csv',index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)
    
for i in range(no_models):
    # Here we load the pre-trained drug layer, which was trained on VCAP data
    model = torch.load('../../results/case_study/models/l1000_latest_model_modeltype4_a375_case_study%s.pt'%i)
    model.eval()

    Yhat, YhatFull = model(X)
    Y_ALL_TRAIN[i, :, :] = Yhat
    
    train_pear[i,:] = pearson_r(Y.detach(), Yhat.detach()).detach().numpy()
    print2log('Finished model %s'%i)

train_pear = pandas.DataFrame(train_pear)
train_pear.columns = TFOutput.columns
train_pear.index = ['model_ '+str(ind) for ind in range(no_models)]
train_pear.to_csv('../../results/case_study/performance/trainPerformance_A375_modeltype4_perTF.csv')

torch.save(Y_ALL_TRAIN,'../../results/case_study/Y_A375_train.pt')

Y_MEAN_TRAIN = torch.mean(Y_ALL_TRAIN, 0)
Y_STD_TRAIN = torch.std(Y_ALL_TRAIN, 0)

pearTrain = pearson_r(Y_MEAN_TRAIN.detach(), Y.detach()).detach().numpy()
dist = beta(Y_MEAN_TRAIN.shape[0] / 2 - 1, Y_MEAN_TRAIN.shape[0] / 2 - 1, loc=-1, scale=2)
pvaluesTrain = [2 * dist.cdf(-abs(r)) for r in list(pearTrain)]
pearTrain = pandas.DataFrame(pearTrain)
pearTrain.index = TFOutput.columns
pearTrain.to_csv('../../results/case_study/performance/trainEnsemblePerformance_A375.csv')
pvaluesTrain = pandas.DataFrame(numpy.array(pvaluesTrain))
pvaluesTrain.index = TFOutput.columns
pvaluesTrain.to_csv('../../results/case_study/performance/trainEnsemblePvalues_A375.csv')

Y_MEAN_TRAIN = pandas.DataFrame(Y_MEAN_TRAIN.detach().numpy())
Y_MEAN_TRAIN.index = TFOutput.index
Y_MEAN_TRAIN.columns = TFOutput.columns
Y_MEAN_TRAIN.to_csv('../../results/case_study/Y_A375_mean_train.csv')
Y_STD_TRAIN = pandas.DataFrame(Y_STD_TRAIN.detach().numpy())
Y_STD_TRAIN.index = TFOutput.index
Y_STD_TRAIN.columns = TFOutput.columns
Y_STD_TRAIN.to_csv('../../results/case_study/Y_A375_std_train.csv')