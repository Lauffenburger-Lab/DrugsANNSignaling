import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
from captum.attr import IntegratedGradients
import argparse
import logging
import time
start_time = time.time()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(device)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Drug-Targets for drugs of interest')
parser.add_argument('--ensembles_path', action='store', default="../results/A375_ensembles/models/")
parser.add_argument('--inputPattern', action='store', default="A375_l1000_latest_model_modeltype4_model")
parser.add_argument('--numberOfModels', action='store', default=50)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
parser.add_argument('--drugInputFile', action='store',default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv')
parser.add_argument('--drugTargetsFile', action='store',default='../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv')
parser.add_argument('--TFOutFile', action='store',default='../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv')
parser.add_argument('--drugSimilarityFile', action='store',default='../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv')
parser.add_argument('--interactionsPath', action='store',default=None)
parser.add_argument('--Y_ALL_path', action='store',default=None)
parser.add_argument('--Y_ALL_masked_path', action='store',default=None)
parser.add_argument('--ig_n_steps', action='store', default=10)

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
ig_n_steps = int(args.ig_n_steps)

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
models_times = []
# Y_ALL = torch.load(Y_ALL_path)
# Y_ALL_masked = torch.load(Y_ALL_masked_path)
# all_scores = np.load('all_scores.npy')
for i in range(numberOfModels):
    prev_time = time.time()
    model = torch.load(inputPath + str(i) + ".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    Y_ALL[i, :, :] = Yhat.detach()
    Xin = model.drugLayer(X)
    # Get starting activity not steady state
    Yin = model.inputLayer(Xin)
    drug_module = model.drugLayer.to(device)
    drug_module.mask = drug_module.mask.to(device)
    drug_module.drugSim = drug_module.drugSim.to(device)
    ig = IntegratedGradients(drug_module)

    print2log('Score calculation with integrated gradients for model %s...'%i)
    ig_start_time = time.time()
    # Per output latent variable input importance translation captum
    # 1st dimesion input
    # 2nd dimesion output
    scores = torch.zeros (model.drugLayer.mask.T.shape)
    kk=1
    for target in range(scores.shape[1]):
        attr, _ = ig.attribute(torch.eye(X.shape[1]).double().to(device),target=target,n_steps=ig_n_steps,return_convergence_delta=True)
        scores[:,target] = torch.diagonal(attr, 0).cpu()
        if kk%50==0 or kk==1:
            print2log('Finished target %s'%kk)
        kk=kk+1
    ig_end_time = time.time()
    ig_runtime = (ig_end_time-ig_start_time)/60
    print2log('Score calculation with integrated gradients for model %s...Done in %s minutes' % (i,ig_runtime))

    interactions = pd.DataFrame(scores.numpy())
    interactions.index = drugInput.columns
    interactions.columns = drugTargets.columns
    interactions.to_csv(interactionsPath + 'interactionScores_' + str(i) + '.csv')
    # interactions = pd.read_csv(outPattern+'_interactionScores_'+str(i)+'.csv',index_col=0)

    scores = torch.abs(scores).detach()
    for j in range(len(thresholds)):
        th = thresholds[j]
        mask = torch.mm((1.0 * (X != 0)).double(), (1.0 * (scores >= th)).double())
        Xin_masked = torch.mul(Xin, mask)
        fullX_masked = model.inputLayer(Xin_masked)
        YhatFull_masked = model.network(fullX_masked)
        Yhat_masked = model.projectionLayer(YhatFull_masked)
        Y_ALL_masked[i, j, :, :] = Yhat_masked.detach()

    torch.save(Y_ALL, Y_ALL_path)
    torch.save(Y_ALL_masked, Y_ALL_masked_path)
    all_scores[i, :, :] = interactions.values
    np.save('all_scores.npy', all_scores)
    models_times.append(time.time() - prev_time)
    
avg_time = np.mean(models_times)/60
print2log('Average time per model = %s minutes'%avg_time)
torch.save(Y_ALL, Y_ALL_path)
torch.save(Y_ALL_masked, Y_ALL_masked_path)
# Calculate average score
scores = torch.tensor(np.mean(all_scores, 0))
scores_std = torch.tensor(np.std(all_scores, 0))
interactions = pd.DataFrame(scores.numpy())
interactions.index = drugInput.columns
interactions.columns = drugTargets.columns
interactions.to_csv(interactionsPath + 'interactionScoresEnsembled.csv')
interactions_std = pd.DataFrame(scores_std.numpy())
interactions_std.index = drugInput.columns
interactions_std.columns = drugTargets.columns
interactions_std.to_csv(interactionsPath + 'interactionScoresSTDEnsembled.csv')

end_time = time.time()
execution_time = end_time - start_time
print2log("Execution time: %s" % execution_time)