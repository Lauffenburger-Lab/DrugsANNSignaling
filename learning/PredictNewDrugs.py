import torch
import numpy
import matplotlib.pyplot as plt
#import bionetwork as onlybionet
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
from rdkit import DataStructs
import argparse
import logging

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
    
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Evaluation of unseen drugs')
parser.add_argument('--inputPattern', action='store', default='A375_l1000_latest_model_modeltype4_model')
parser.add_argument('--ensembles_path', action='store', default='../results/A375_ensembles/models/')
parser.add_argument('--DrugsIn', action='store', default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv')
parser.add_argument('--TargetsIn', action='store', default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv')
parser.add_argument('--TFsOut', action='store', default='../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv')
parser.add_argument('--ChemicalSims', action='store', default='../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv')
parser.add_argument('--PKN', action='store', default='../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
parser.add_argument('--PknAnnotation', action='store', default='../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv')
parser.add_argument('--res_dir', action='store', default='../results/A375_ensembles/models/')
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
parser.add_argument('--new_drugs', action='store',default=None)

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
inputPath = ensembles_path + 'models/' + inputPattern


NewDrugs = pd.read_csv(new_drugs,index_col=0)
NewTFOutput = pd.read_csv(new_drugs_output,index_col=0)
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

### Calculate ECFP4 Similarity between new drugs and old drugs and make SimInput matrix
SimilarityMat = torch.zeros(NewDrugs.shape[1],drugSim.shape[1])
for new_smi in range(NewDrugs.shape[1]):
    new_smi = NewDrugs.columns.values[i]
    new_mol = Chem.MolFromSmiles(new_smi)
    new_fp = AllChem.GetMorganFingerprintAsBitVect(new_mol, 2, 1024)
    for j in range(drugSim.shape[1]):
        smi = drugSim.columns.values[j]
        mol = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
        SimilarityMat[i,j] = DataStructs.FingerprintSimilarity(new_fp, fp)
print2log('Calculated similarity between new drugs and ones used in training')

X_new =
Y_new = torch.tensor(TFOutput.values, dtype=torch.double)

#### Random shuffle X
#X = X[torch.randperm(X.size()[0]),:]

sampleName = drugInput.index.values
mask_dt =torch.tensor(numpy.array(drugTargets)).transpose(0,1)
val_pear = numpy.zeros((no_models, TFOutput.shape[1]))
Y_ALL_VAL = torch.zeros((no_models,TFOutput_val.shape[0],TFOutput_val.shape[1]))
for i in range(no_models):
    # Here we load the pre-trained drug layer, which was trained on VCAP data
    model = torch.load(inputPattern + 'model%s.pt'%i)
    model.eval()

    X_targets = model.drugLayer(torch.eye(SimilarityMat.shape[1]).double())
    Xin = torch.mm(SimilarityMat, X_targets) # to account for the similarity of drugs with other drugs, it is also the scores for drug-target interactions
    X_in_scaled = torch.mm(X_new,Xin)
    fullX = model.inputLayer(X_in_scaled)
    YhatFull_val = model.network(fullX)
    Yhat_val = model.projectionLayer(YhatFull_val)
    Y_ALL_VAL[i, :, :] = Yhat_val

    val_pear[i, :] = pearson_r(Y_new.detach(), Yhat_val.detach()).detach().numpy()
    if len(numpy.where(numpy.isnan(val_pear[i, :]))[0]):
        inds = numpy.where(numpy.isnan(val_pear[i, :]))[0]
    r = pearson_r(Y_new[:, inds].detach(),
                  Yhat_val[:, inds].detach() + 1e-8 * torch.randn(Yhat_val.shape[0], Yhat_val[:, inds].shape[1]))
    val_pear[i, inds] = r.detach().numpy()

val_pear = pandas.DataFrame(val_pear)
val_pear.columns = TFOutput_val.columns
val_pear.index = ['model_ ' + str(ind) for ind in range(no_models)]
val_pear.to_csv(res_dir + 'valPerformance_newDrugs_perTF.csv')

torch.save(Y_ALL_VAL, res_dir + 'Y_all_newDrugs.pt')

Y_MEAN_TRAIN = torch.mean(Y_ALL_TRAIN, 0)
Y_STD_TRAIN = torch.std(Y_ALL_TRAIN, 0)

Y_MEAN_VAL = torch.mean(Y_ALL_VAL, 0)
Y_STD_VAL = torch.std(Y_ALL_VAL, 0)

pearVal = pearson_r(Y_MEAN_VAL.detach(), Y_new.detach()).detach().numpy()
dist = beta(Y_MEAN_VAL.shape[0] / 2 - 1, Y_MEAN_VAL.shape[0] / 2 - 1, loc=-1, scale=2)
pvaluesVal = [2 * dist.cdf(-abs(r)) for r in list(pearVal)]
pearVal = pandas.DataFrame(pearVal)
pearVal.index = TFOutput_val.columns
pearVal.to_csv(res_dir + 'valEnsemblePerformance_newDrugs.csv')
pvaluesVal = pandas.DataFrame(numpy.array(pvaluesVal))
pvaluesVal.index = TFOutput_val.columns
pvaluesVal.to_csv(res_dir + 'valEnsemblePvalues_newDrugs.csv')
Y_MEAN_VAL = pandas.DataFrame(Y_MEAN_VAL.detach().numpy())
Y_MEAN_VAL.index = TFOutput_val.index
Y_MEAN_VAL.columns = TFOutput_val.columns
Y_MEAN_VAL.to_csv(res_dir + 'Y_newDrugs_mean.csv')
Y_STD_VAL = pandas.DataFrame(Y_STD_VAL.detach().numpy())
Y_STD_VAL.index = TFOutput_val.index
Y_STD_VAL.columns = TFOutput_val.columns
Y_STD_VAL.to_csv(res_dir + 'Y_newDrugs_std.csv')