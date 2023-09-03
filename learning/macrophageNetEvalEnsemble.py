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
for c,cell in enumerate(cell_lines):
    #Load network
    networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
    annotation = pandas.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
    uniprot2gene = dict(zip(annotation['code'], annotation['name']))
    bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
    spectralCapacity = numpy.exp(numpy.log(1e-2)/bionetParams['iterations'])

    ### Load all data
    drugInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False,
                                index_col=0)
    drugSmiles = drugInput.columns.values
    drugTargets = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
    TFOutput = pandas.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
    cellInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
    TFOutput = TFOutput.loc[cellInput[cellInput[cell] == 1].index, :]
    drugInput = drugInput.loc[cellInput[cellInput[cell] == 1].index, :]
    sigs_to_drop = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/smiles_and_sigs_to_drop.csv', index_col=0)
    TFOutput_train = TFOutput.loc[~TFOutput.index.isin(sigs_to_drop.sig_id), :]
    TFOutput_val = TFOutput.loc[TFOutput.index.isin(sigs_to_drop.sig_id), :]

    train_pear = numpy.zeros((no_models,TFOutput.shape[1]))
    val_pear = numpy.zeros((no_models, TFOutput.shape[1]))
    train_var = numpy.zeros((no_models, TFOutput.shape[1]))
    Y_ALL_VAL = torch.zeros((no_models,TFOutput_val.shape[0],TFOutput_val.shape[1]))
    Y_ALL_TRAIN = torch.zeros((no_models,TFOutput_train.shape[0],TFOutput_train.shape[1]))
    for i in range(no_models):
        # Here we load the pre-trained drug layer, which was trained on VCAP data
        model1 = torch.load('../results/FinalEnsemble/models/VCAP_l1000_latest_modeltype4_model%s.pt'%i)
        model1.eval()
        model = torch.load('../results/FinalEnsemble/models/'+cell+'_l1000_latest_modeltype4_model'+str(i)+'.pt')
        model.eval()

        TFOutput_train=TFOutput.loc[~TFOutput.index.isin(sigs_to_drop.sig_id),:]
        drugInput_train=drugInput.loc[~drugInput.index.isin(sigs_to_drop.sig_id),:]
        TFOutput_val = TFOutput.loc[TFOutput.index.isin(sigs_to_drop.sig_id), :]
        drugInput_val = drugInput.loc[drugInput.index.isin(sigs_to_drop.sig_id), :]

        #Subset input and output to intersecting modeldes
        druginName = drugInput_train.columns.values
        inName = drugTargets.columns.values
        outconds_train = drugInput_train.index.values
        outconds_val = drugInput_val.index.values
        outName = TFOutput_train.columns.values
        inName = numpy.intersect1d(nodeNames, inName)
        outName = numpy.intersect1d(nodeNames, outName)
        inNameGene = [uniprot2gene[x] for x in inName]
        outNameGene = [uniprot2gene[x] for x in outName]
        TFOutput_train = TFOutput_train.loc[outconds_train,outName]
        TFOutput_val = TFOutput_val.loc[outconds_val, outName]
        drugTargets = drugTargets.loc[:,inName]
        if ConvertToEmpProb==True:
            print2log('Convereted to Emprirical probabilities')
            TFOutput = 1/(1+numpy.exp(-TFOutput))

        drugTargets = drugTargets.loc[drugInput_train.columns.values,:]
        drugSim = pandas.read_csv('../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv',index_col=0)
        drugSim = drugSim.loc[drugInput_train.columns.values,drugInput_train.columns.values]
        drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

        X_train = torch.tensor(drugInput_train.values.copy(), dtype=torch.double)
        Y_train = torch.tensor(TFOutput_train.values, dtype=torch.double)
        X_val = torch.tensor(drugInput_val.values.copy(), dtype=torch.double)
        Y_val = torch.tensor(TFOutput_val.values, dtype=torch.double)

        model1.eval()
        Xin_train = model1.drugLayer(X_train)
        Xin_train = Xin_train.detach()
        Xin_val = model1.drugLayer(X_val)
        Xin_val = Xin_val.detach()

        model.eval()
        Yhat_train, YhatFull_train = model(Xin_train)
        Yhat_val, YhatFull_val = model(Xin_val)
        Y_ALL_VAL[i,:,:] = Yhat_val
        Y_ALL_TRAIN[i, :, :] = Yhat_train


        train_pear[i,:] = pearson_r(Y_train.detach(), Yhat_train.detach()).detach().numpy()
        val_pear[i,:] = pearson_r(Y_val.detach(), Yhat_val.detach()).detach().numpy()
        # train_pear[i, :] = pearson_r(torch.clamp(Y_train.detach(), 0, 1),
        #                              torch.clamp(Yhat_train.detach(), 0, 1)).detach().numpy()
        # val_pear[i, :] = pearson_r(torch.clamp(Y_val.detach(), 0, 1),
        #                            torch.clamp(Yhat_val.detach(), 0, 1)).detach().numpy()

        train_var[i,:] = torch.std(Y_train.detach(),dim=0 ,unbiased=True).detach().numpy()

        # plt.figure()
        # plt.scatter(Yhat.detach().numpy(), Y.detach().numpy(), alpha=0.2)
        # plotting.lineOfIdentity()
        # plotting.addCorrelation(Yhat, Y)
        # plt.xlabel('Fit')
        # plt.ylabel('Data')
        # plt.gca().axis('equal')
        # plt.gca().set_xticks([0, 0.5, 1])
        # plt.gca().set_yticks([0, 0.5, 1])
        # plt.savefig('../results/FinalEnsemble/test/performance2_'+cell+'_l1000_latest_modeltype4_shuffleY'+ str(i)+'.png',dpi=600)

        # plt.figure()
        # rank = plotting.compareAllTFs(Yhat_val.detach().cpu(), Y_val.detach().cpu(), outNameGene)
        # plt.savefig('../results/FinalEnsemble/test/performancePerTF_' + cell +'_l1000_latest_modeltype4_shuffleY'+ str(i)+'.png', dpi=600)
        # plt.figure()
        # rank = plotting.compareAllTFs(Yhat_val.detach().cpu().T, Y_val.detach().cpu().T, outconds_val)
        # plt.savefig('../results/FinalEnsemble/test/performancePerSample_' + cell +'_l1000_latest_modeltype4_shuffleY'+ str(i)+'.png', dpi=600)


        # R2[c,i] = Rsquared(Yhat.detach().cpu(),Y.detach().cpu()).item()
        # if i==2:
        #     print2log(Yhat_val[:,38])
        #     print2log(Y_val[:,38])
        #     print2log(pearson_r(Y_val.detach(), Yhat_val.detach()).detach().numpy()[38])
        if len(numpy.where(numpy.isnan(val_pear[i,:]))[0]):
            inds = numpy.where(numpy.isnan(val_pear[i,:]))[0]
            r = pearson_r(Y_val[:,inds].detach(), Yhat_val[:,inds].detach() + 1e-8 * torch.randn(Yhat_val.shape[0],Yhat_val[:,inds].shape[1]))
            val_pear[i, inds] = r.detach().numpy()

        # uniq_index = uniq_index +1

        # print2log(len(numpy.where(numpy.isnan(train_pear))[0]))
    # plt.figure()
    # plt.scatter(train_pear, val_pear, alpha=0.2)
    # plt.xlim([numpy.min(numpy.concatenate((train_pear,val_pear))), 1])
    # plt.ylim([numpy.min(numpy.concatenate((train_pear,val_pear))), 1])
    # lineOfIdentity()
    # r, p = pearsonr(train_pear.flatten(), val_pear.flatten())
    # plt.text(numpy.min(train_pear), 0.8, 'r {:.2f}, p.value = {:.4f}'.format(r, p))
    # plt.xlabel('train correlation')
    # plt.ylabel('validation correlation')
    # plt.gca().axis('equal')
    # plt.savefig('../results/FinalEnsemble/test/performance_val_vs_train_' + cell + '_l1000_latest_modeltype4.png',dpi=600)
    #
    # plt.figure()
    # plt.scatter(train_var, val_pear, alpha=0.2)
    # r, p = pearsonr(train_var.flatten(), val_pear.flatten())
    # plt.text(numpy.min(train_var)*1.2, 0.6, 'r {:.2f}, p.value = {:.4f}'.format(r,p))
    # # lineOfIdentity()
    # plt.xlabel('train s.t.d')
    # plt.ylabel('validation correlation')
    # plt.savefig('../results/FinalEnsemble/test/performanceVal_vs_varTrain_' + cell + '_l1000_latest_modeltype4.png',dpi=600)
    #
    # plt.figure()
    # plt.scatter(train_var, train_pear, alpha=0.2)
    # r, p = pearsonr(train_var.flatten(), train_pear.flatten())
    # plt.text(numpy.min(train_var)*1.2, 0.6, 'r {:.2f}, p.value = {:.4f}'.format(r,p))
    # # lineOfIdentity()
    # plt.xlabel('train s.t.d')
    # plt.ylabel('train correlation')
    # plt.savefig('../results/FinalEnsemble/test/performanceTrain_vs_varTrain_' + cell + '_l1000_latest_modeltype4.png',dpi=600)

    train_pear = pandas.DataFrame(train_pear)
    train_pear.columns = TFOutput_train.columns
    train_pear.index = ['model_ '+str(ind) for ind in range(no_models)]
    train_pear.to_csv('../results/FinalEnsemble/test/trainPerformance_' + cell + '_modeltype4_perTF.csv')

    train_var = pandas.DataFrame(train_var)
    train_var.columns = TFOutput_train.columns
    train_var.index = ['model_ '+str(ind) for ind in range(no_models)]
    train_var.to_csv('../results/FinalEnsemble/test/trainVariance_' + cell + '_modeltype4_perTF.csv')

    val_pear = pandas.DataFrame(val_pear)
    val_pear.columns = TFOutput_val.columns
    val_pear.index = ['model_ '+str(ind) for ind in range(no_models)]
    val_pear.to_csv('../results/FinalEnsemble/test/valPerformance_' + cell + '_modeltype4_perTF.csv')

    torch.save(Y_ALL_TRAIN,'../results/FinalEnsemble/preds/Y_'+cell+'_all_train.pt')
    torch.save(Y_ALL_VAL, '../results/FinalEnsemble/preds/Y_' + cell + '_all_val.pt')

    # Y_ALL_TRAIN = torch.load('../results/FinalEnsemble/preds/Y_' + cell + '_all_train.pt')
    # Y_ALL_VAL = torch.load('../results/FinalEnsemble/preds/Y_' + cell + '_all_val.pt')

    Y_MEAN_TRAIN = torch.mean(Y_ALL_TRAIN, 0)
    Y_STD_TRAIN = torch.std(Y_ALL_TRAIN, 0)

    Y_MEAN_VAL = torch.mean(Y_ALL_VAL, 0)
    Y_STD_VAL = torch.std(Y_ALL_VAL, 0)

    pearVal = pearson_r(Y_MEAN_VAL.detach(), Y_val.detach()).detach().numpy()
    dist = beta(Y_MEAN_VAL.shape[0] / 2 - 1, Y_MEAN_VAL.shape[0] / 2 - 1, loc=-1, scale=2)
    pvaluesVal = [2*dist.cdf(-abs(r)) for r in list(pearVal)]
    pearVal = pandas.DataFrame(pearVal)
    pearVal.index = TFOutput_val.columns
    pearVal.to_csv('../results/FinalEnsemble/test/valEnsemblePerformance_' + cell + '.csv')
    pvaluesVal = pandas.DataFrame(numpy.array(pvaluesVal))
    pvaluesVal.index = TFOutput_val.columns
    pvaluesVal.to_csv('../results/FinalEnsemble/test/valEnsemblePvalues_' + cell + '.csv')
    pearTrain = pearson_r(Y_MEAN_TRAIN.detach(), Y_train.detach()).detach().numpy()
    dist = beta(Y_MEAN_TRAIN.shape[0] / 2 - 1, Y_MEAN_TRAIN.shape[0] / 2 - 1, loc=-1, scale=2)
    pvaluesTrain = [2 * dist.cdf(-abs(r)) for r in list(pearTrain)]
    pearTrain = pandas.DataFrame(pearTrain)
    pearTrain.index = TFOutput_train.columns
    pearTrain.to_csv('../results/FinalEnsemble/test/trainEnsemblePerformance_' + cell + '.csv')
    pvaluesTrain = pandas.DataFrame(numpy.array(pvaluesTrain))
    pvaluesTrain.index = TFOutput_train.columns
    pvaluesTrain.to_csv('../results/FinalEnsemble/test/trainEnsemblePvalues_' + cell + '.csv')

    Y_MEAN_TRAIN = pandas.DataFrame(Y_MEAN_TRAIN.detach().numpy())
    Y_MEAN_TRAIN.index = TFOutput_train.index
    Y_MEAN_TRAIN.columns = TFOutput_train.columns
    Y_MEAN_TRAIN.to_csv('../results/FinalEnsemble/preds/Y_' + cell + '_mean_train.csv')
    Y_STD_TRAIN = pandas.DataFrame(Y_STD_TRAIN.detach().numpy())
    Y_STD_TRAIN.index = TFOutput_train.index
    Y_STD_TRAIN.columns = TFOutput_train.columns
    Y_STD_TRAIN.to_csv('../results/FinalEnsemble/preds/Y_' + cell + '_std_train.csv')

    Y_MEAN_VAL = pandas.DataFrame(Y_MEAN_VAL.detach().numpy())
    Y_MEAN_VAL.index = TFOutput_val.index
    Y_MEAN_VAL.columns = TFOutput_val.columns
    Y_MEAN_VAL.to_csv('../results/FinalEnsemble/preds/Y_' + cell + '_mean_val.csv')
    Y_STD_VAL = pandas.DataFrame(Y_STD_VAL.detach().numpy())
    Y_STD_VAL.index = TFOutput_val.index
    Y_STD_VAL.columns = TFOutput_val.columns
    Y_STD_VAL.to_csv('../results/FinalEnsemble/preds/Y_' + cell + '_std_val.csv')

    print2log('Finished cell-line %s' % cell)
