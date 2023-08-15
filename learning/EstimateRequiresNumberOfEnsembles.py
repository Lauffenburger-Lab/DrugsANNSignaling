import torch
import numpy
import random
import math
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import bionetwork as SingleBionet
import bionetworkWithDrugs as bionetwork
#import bionetwork
import plotting
import pandas
import saveSimulations
from itertools import combinations
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

def get_combinations(elements, k):
    """
    Returns a list of all possible combinations of k elements from the given list of elements.
    """
    return list(combinations(elements, k))

cell_lines = ['A375','A549','HA1E','HCC515','HT29','PC3','MCF7','HEPG2']

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

no_models = 51
models_set = list(range(no_models))
# R2 = numpy.zeros((len(cell_lines),5))
# model1 = torch.load('../results/FinalEnsemble/VCAP_l1000_latest_modeltype4_no0.pt')
# model1.eval()
pearVal_all = pandas.DataFrame({'r':[],'ensembles':[],'combo_id':[],'cell':[]})
pearTrain_all = pandas.DataFrame({'r':[],'ensembles':[],'combo_id':[],'cell':[]})
for c,cell in enumerate(cell_lines):

    Y_ALL_TRAIN = torch.load('../results/FinalEnsemble/preds/Y_' + cell + '_all_train.pt')
    Y_ALL_VAL = torch.load('../results/FinalEnsemble/preds/Y_' + cell + '_all_val.pt')

    #Load network
    networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('data/l1000_lvl3_withsignor-Model.tsv')
    annotation = pandas.read_csv('data/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
    uniprot2gene = dict(zip(annotation['code'], annotation['name']))
    bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
    spectralCapacity = numpy.exp(numpy.log(1e-2)/bionetParams['iterations'])

    ### Load all data
    drugInput = pandas.read_csv('data/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False,
                                index_col=0)
    drugSmiles = drugInput.columns.values
    drugTargets = pandas.read_csv('data/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
    TFOutput = pandas.read_csv('data/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
    cellInput = pandas.read_csv('data/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
    TFOutput = TFOutput.loc[cellInput[cellInput[cell] == 1].index, :]
    drugInput = drugInput.loc[cellInput[cellInput[cell] == 1].index, :]
    # sigs_to_drop = pandas.read_csv('../smiles_and_sigs_to_drop_easy.csv', index_col=0)
    sigs_to_drop = pandas.read_csv('../smiles_and_sigs_to_drop.csv', index_col=0)
    TFOutput_train = TFOutput.loc[~TFOutput.index.isin(sigs_to_drop.sig_id), :]
    TFOutput_val = TFOutput.loc[TFOutput.index.isin(sigs_to_drop.sig_id), :]

    TFOutput_train = TFOutput.loc[~TFOutput.index.isin(sigs_to_drop.sig_id), :]
    drugInput_train = drugInput.loc[~drugInput.index.isin(sigs_to_drop.sig_id), :]
    TFOutput_val = TFOutput.loc[TFOutput.index.isin(sigs_to_drop.sig_id), :]
    drugInput_val = drugInput.loc[drugInput.index.isin(sigs_to_drop.sig_id), :]
    # Subset input and output to intersecting modeldes
    druginName = drugInput_train.columns.values
    inName = drugTargets.columns.values
    outconds_train = drugInput_train.index.values
    outconds_val = drugInput_val.index.values
    outName = TFOutput_train.columns.values
    inName = numpy.intersect1d(nodeNames, inName)
    outName = numpy.intersect1d(nodeNames, outName)
    inNameGene = [uniprot2gene[x] for x in inName]
    outNameGene = [uniprot2gene[x] for x in outName]
    TFOutput_train = TFOutput_train.loc[outconds_train, outName]
    TFOutput_val = TFOutput_val.loc[outconds_val, outName]
    drugTargets = drugTargets.loc[:, inName]
    if ConvertToEmpProb == True:
        print2log('Convereted to Emprirical probabilities')
        TFOutput = 1 / (1 + numpy.exp(-TFOutput))
    drugTargets = drugTargets.loc[drugInput_train.columns.values, :]
    drugSim = pandas.read_csv('../out_lvl3_similaritiess.csv', index_col=0)
    drugSim = drugSim.loc[drugInput_train.columns.values, drugInput_train.columns.values]
    drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)
    X_train = torch.tensor(drugInput_train.values.copy(), dtype=torch.double)
    Y_train = torch.tensor(TFOutput_train.values, dtype=torch.double)
    X_val = torch.tensor(drugInput_val.values.copy(), dtype=torch.double)
    Y_val = torch.tensor(TFOutput_val.values, dtype=torch.double)

    for no_of_ensembles in range(1,no_models+1):
        pear = []
        for i in range(400):
            sampled_models = random.sample(models_set,no_of_ensembles)
            models = numpy.array(sampled_models)
            Y_hats_train = Y_ALL_TRAIN[models,:,:]
            Y_hats_val = Y_ALL_VAL[models,:,:]

            Y_MEAN_TRAIN = torch.mean(Y_hats_train,0)
            Y_MEAN_VAL = torch.mean(Y_hats_val, 0)

            pearVal = pearson_r(Y_MEAN_VAL.detach(), Y_val.detach()).detach().numpy()
            pear.append(numpy.nanmean(pearVal))
            pearVal = pandas.DataFrame(pearVal)
            pearVal.index = TFOutput_val.columns
            pearVal.columns = ['r']
            pearVal['ensembles'] = no_of_ensembles
            # pearVal['combo_id'] = j
            pearVal['trial'] = i
            pearVal['cell'] = cell
            pearVal_all = pearVal_all.append(pearVal)
            # pearVal.to_csv('../results/FinalEnsemble/EstimateEnsembles/valEnsemblePerformance_ensembles%s_combo%s_'%(no_of_ensembles,j) + cell  +'.csv')
            pearTrain = pearson_r(Y_MEAN_TRAIN.detach(), Y_train.detach()).detach().numpy()
            pearTrain = pandas.DataFrame(pearTrain)
            pearTrain.index = TFOutput_train.columns
            pearTrain.columns = ['r']
            pearTrain['ensembles'] = no_of_ensembles
            # pearTrain['combo_id'] = j
            pearTrain['trial'] = i
            pearTrain['cell'] = cell
            pearTrain_all = pearTrain_all.append(pearTrain)
            # pearTrain.to_csv('../results/FinalEnsemble/EstimateEnsembles/trainEnsemblePerformance_ensembles%s_combo%s_'%(no_of_ensembles,j) + cell + '.csv')
        print2log('Pearson: %s for number of ensembles %s'%(numpy.nanmean(pear),no_of_ensembles))
    print2log('Finished cell-line %s'%cell)
pearTrain_all.to_csv('../results/FinalEnsemble/EstimateEnsembles/trainEnsemblePerformance_ensembles_all.csv')
pearVal_all.to_csv('../results/FinalEnsemble/EstimateEnsembles/valEnsemblePerformance_ensembles_all.csv')

print2log(pearVal_all)