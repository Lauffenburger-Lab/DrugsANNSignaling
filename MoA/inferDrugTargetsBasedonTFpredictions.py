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

parser = argparse.ArgumentParser(prog='Infer off-target drug-target interactions')
parser.add_argument('--ensembles_path', action='store', default="CVL1000_Paper/A375_ensembles/")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_model_modeltype4_model")
parser.add_argument('--cell_line', action='store', default="A375")
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
TF = "Q08050"
TF_gene = "FOXM1"
drug = "C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13"
drug_name = "lestaurtinib"
sample = "CPC014_A375_6H:BRD-K23192422-001-01-1:10"
### Define gradient thresholds
thresholds = list(np.logspace(-3.5, 3.5, num=50))
# thresholds = list(np.logspace(-5, 5, num=50))
# print2log('Thresholds:',thresholds)

### Load network
#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('data/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('data/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for A549 was 120
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

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)
#Xin = torch.mm(X,drugLayer)

# Keep only drug/sample of interest
TF_ind = np.where(TFOutput.columns==TF)[0]
sample_ind = np.where(TFOutput.index==sample)[0]

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
    Y_ALL[i,:,:] = Yhat.detach()
    Xin = model.drugLayer(X)
    # print2log(Xin.shape)
    
    # Get starting activity not steady state
    Yin = model.inputLayer(Xin)
    # X_binned =  drug_ratioMatrix * X.unsqueeze(2)
    # Xin_range = torch.zeros(X_binned.shape[1],Xin.shape[1]).double()
    # drug_sum = torch.zeros(X_binned.shape[1]).double()
    # for k in range(X_binned.shape[0]):
    #     Xin_binned =  model.drugLayer(X_binned[k,:,:].T)
    #     drug_ind = torch.where(X[k,:]!=0)[0]
    #     Xin_range[drug_ind,:] = Xin_range[drug_ind,:] + torch.abs(torch.max(Xin_binned,0)[0] - torch.min(Xin_binned,0)[0]).unsqueeze(0)
    #     drug_sum[drug_ind] = drug_sum[drug_ind] + 1.0
    # Xin_range = Xin_range.T / drug_sum
    # Xin_range = Xin_range.T
    # load interaction scores
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/l1000_modeltype4_lamda6_' + cell + '_interactionScores_%s.csv' % i,index_col=0)
    scores = torch.tensor(interactions.values).double()
    # scores = scores * Xin_range
    scores = torch.abs(scores).detach()
    # mask = torch.mm((1.0 * (X != 0)).double(), (1.0 * (scores >= 0.3)).double())
    # print2log(mask.shape)
    # print2log(mask)

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


torch.save(Y_ALL,ensembles_path+'preds/Y_'+cell+'_ALL.pt')
torch.save(Y_ALL_masked,ensembles_path+'preds/Y_'+cell+'_ALL_MASKED.pt')
torch.save(filtered_all,ensembles_path+'preds/filtered_'+cell+'_ALL_MASKED.pt')
# Y_ALL = torch.load(ensembles_path+'preds/Y_'+cell+'_ALL.pt')
# Y_ALL_masked = torch.load(ensembles_path+'preds/Y_'+cell+'_ALL_MASKED.pt')
# filtered_all = torch.load(ensembles_path+'preds/filtered_'+cell+'_ALL_MASKED.pt')

l = []
plt.figure()
for i in range(numberOfModels):
    # Calculate error seperetaly
    # global_error = torch.mean(torch.abs(Y_ALL_masked[i,:,sample_ind,:]-Y[sample_ind,:]),1).detach().squeeze().numpy()
    TF_error = torch.abs(Y_ALL_masked[i,:,sample_ind,TF_ind]-Y[sample_ind,TF_ind]).detach().squeeze().numpy()
    unmasked_TF_error = torch.abs(Y_ALL[i,sample_ind, TF_ind] - Y[sample_ind, TF_ind]).detach().squeeze().numpy()
    thresh = 0.25
    # TF_perc_err = (TF_error - TF_error[0]) / (TF_error[-1] - TF_error[0])
    TF_perc_err = (TF_error-unmasked_TF_error)/(TF_error[-1]-unmasked_TF_error)
    if np.where(TF_perc_err <= thresh)[0].shape[0]>0:
        tf_thresh_selection = np.where(TF_perc_err <= thresh)[0][-1]
    else:
        tf_thresh_selection = 0
    l.append(thresholds[tf_thresh_selection])
    #plot results
    # plt.plot(np.array(thresholds),global_error,'o-',label ="Global error %s"%i)
    plt.plot(np.array(thresholds),TF_error,'o-',label =TF_gene+" error %s"%i)
    plt.plot([0, thresholds[tf_thresh_selection]], np.array([1, 1]) * TF_error[tf_thresh_selection], 'red',
             linestyle='--')
    plt.plot(np.array([1, 1]) * thresholds[tf_thresh_selection], [0, TF_error[tf_thresh_selection]], 'red',
             linestyle='--')
plt.xscale('log')
plt.xlabel("absolute gradient theshold")
plt.ylabel("error")
plt.legend()
plt.savefig(ensembles_path+"InteractionScores/error_vs_threshold_allmodels_"+drug+"_"+TF+".png", bbox_inches='tight', dpi=600)
print2log('Average TF threshold %s'%np.mean(l))

l = []
plt.figure()
for i in range(numberOfModels):
    # Calculate error seperetaly
    global_error = torch.mean(torch.abs(Y_ALL_masked[i,:,sample_ind,:]-Y[sample_ind,:]).squeeze(),1).detach().squeeze().numpy()
    unmasked_global_error = torch.mean(torch.abs(Y_ALL[i,sample_ind, :] - Y[sample_ind, :]), 1).detach().squeeze().numpy()
    # global_perc_err = (global_error - global_error[0]) / (global_error[-1] - global_error[0])
    global_perc_err = (global_error-unmasked_global_error)/(global_error[-1]-unmasked_global_error)
    if np.where(global_perc_err <= thresh)[0].shape[0]>0:
        global_thresh_selection = np.where(global_perc_err <= thresh)[0][-1]
    else:
        global_thresh_selection = 0
    l.append(thresholds[global_thresh_selection])
    # TF_error = torch.abs(Y_ALL_masked[i,:,sample_ind,TF_ind]-Y[sample_ind,TF_ind]).detach().squeeze().numpy()
    #plot results
    plt.plot(np.array(thresholds),global_error,'o-',label ="Global error %s"%i)
    plt.plot([0, thresholds[global_thresh_selection]], np.array([1, 1]) * global_error[global_thresh_selection], 'red',
              linestyle='--')
    plt.plot(np.array([1, 1]) * thresholds[global_thresh_selection], [0, global_error[global_thresh_selection]], 'red',
              linestyle='--')
    # plt.plot(np.array(thresholds),TF_error,'o-',label =TF_gene+" error %s"%i)
plt.xscale('log')
plt.xlabel("absolute gradient theshold")
plt.ylabel("error")
plt.legend()
plt.savefig(ensembles_path+"InteractionScores/global_error_vs_threshold_allmodels_"+drug+"_"+TF+".png", bbox_inches='tight', dpi=600)
print2log('Average global threshold %s'%np.mean(l))
# torch.save(torch.tensor(np.array(l)),ensembles_path+"InteractionScores/global_gradient_scores_"+drug+"_sample_ind"+str(sample_ind[0])+"_all_models_withsignalRange.pt")

# Calculate mean predictions
Yhat = torch.mean(Y_ALL,0)
Yhat_masked = torch.mean(Y_ALL_masked,0)


Yhat = Yhat[sample_ind,:]
Yhat_masked = Yhat_masked[:,sample_ind,:]
Yhat_masked = Yhat_masked.squeeze()


# Calculate error
global_error = torch.mean(torch.abs(Yhat_masked-Y[sample_ind,:]),1).detach().squeeze().numpy()
unmasked_global_error = torch.mean(torch.abs(Yhat-Y[sample_ind,:]),1).detach().squeeze().numpy()
TF_error = torch.abs(Yhat_masked[:,TF_ind]-Y[sample_ind,TF_ind]).detach().squeeze().numpy()
unmasked_TF_error = torch.abs(Yhat[0,TF_ind]-Y[sample_ind,TF_ind]).detach().squeeze().numpy()

thresh = 0.25
TF_perc_err = (TF_error-unmasked_TF_error)/(TF_error[-1]-unmasked_TF_error)
# TF_perc_err = (TF_error-TF_error[0])/(TF_error[0])
if np.where(TF_perc_err <= thresh)[0].shape[0]>0:
    tf_thresh_selection = np.where(TF_perc_err <= thresh)[0][-1]
else:
    tf_thresh_selection = 0
global_perc_err = (global_error-unmasked_global_error)/(global_error[-1]-unmasked_global_error)
# global_perc_err = (global_error-global_error[0])/(global_error[0])
if np.where(global_perc_err <= thresh)[0].shape[0]>0:
    global_thresh_selection = np.where(global_perc_err <= thresh)[0][-1]
else:
    global_thresh_selection = 0
print2log('Ensemble global gradient threshold %s'%thresholds[global_thresh_selection])
print2log('Ensemble TF gradient threshold %s'%thresholds[tf_thresh_selection])

#plot results
plt.figure()
plt.plot(np.array(thresholds),global_error,'o-',label ="Global error")
plt.plot(np.array(thresholds),TF_error,'o-',label =TF_gene+" error")
plt.plot([0, thresholds[global_thresh_selection]], np.array([1, 1])*global_error[global_thresh_selection], 'red', linestyle='--')
plt.plot(np.array([1, 1])*thresholds[global_thresh_selection], [0, global_error[global_thresh_selection]], 'red', linestyle='--')
plt.plot([0, thresholds[tf_thresh_selection]], np.array([1, 1])*TF_error[tf_thresh_selection], 'red', linestyle='--')
plt.plot(np.array([1, 1])*thresholds[tf_thresh_selection], [0, TF_error[tf_thresh_selection]], 'red', linestyle='--')
plt.xscale('log')
plt.xlabel("absolute gradient theshold")
plt.ylabel("error")
plt.title('Error for masking increasingly more singal for: '+ drug_name[0].upper()+drug_name[1:len(drug_name)])
plt.legend()
plt.savefig(ensembles_path+"InteractionScores/error_vs_threshold_"+drug+"_"+TF+"_withsignalRange.png", bbox_inches='tight', dpi=600)

# plt.figure()
# plt.errorbar(thresholds,torch.mean(filtered_all,0).detach().numpy(),torch.std(filtered_all,0).detach().numpy(),
#               linestyle='--',marker='.',capsize=4, elinewidth=1.5,ecolor='black')
# plt.plot([0, thresholds[global_thresh_selection]], np.array([1, 1])*torch.mean(filtered_all,0).detach().numpy()[global_thresh_selection], 
#           'red', linestyle='--')
# plt.plot(np.array([1, 1])*thresholds[global_thresh_selection], [0, torch.mean(filtered_all,0).detach().numpy()[global_thresh_selection]],
#           'red', linestyle='--')
# plt.xscale('log')
# plt.xlabel("absolute gradient theshold")
# plt.ylabel("number of interactions inferred")
# plt.legend()
# plt.savefig(ensembles_path+"InteractionScores/filtered_vs_threshold_"+drug+"_withsignalRange.png", bbox_inches='tight', dpi=600)

# plt.hist(Yhat.flatten().detach().numpy())
# # plt.hist(Yin.flatten().detach().numpy()[Yin.flatten().detach().numpy()!=0],10000)
# plt.xscale('symlog')