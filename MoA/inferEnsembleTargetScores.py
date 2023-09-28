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

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

parser = argparse.ArgumentParser(prog='Cell-line simulation')
parser.add_argument('--ensembles_path', action='store', default="../results/FinalEnsemble/")
parser.add_argument('--outPattern', action='store', default="InteractionScores/l1000_modeltype4_lamda6")
parser.add_argument('--inputPattern', action='store', default="l1000_latest_modeltype4_model")
parser.add_argument('--cell_line', action='store', default="VCAP")
parser.add_argument('--numberOfModels', action='store', default=51)
parser.add_argument('--ConvertToEmpProb', action='store',default=False)
args = parser.parse_args()
ensembles_path = args.ensembles_path
outPattern = args.outPattern
cell = args.cell_line
outPattern = ensembles_path + outPattern + "_" + cell
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

model = torch.load(inputPath+str(0)+".pt")
all_scores = np.zeros((numberOfModels,model.drugLayer.mask.T.shape[0],model.drugLayer.mask.T.shape[1]))
mask = model.drugLayer.mask.T
for i in range(numberOfModels):
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    ig = IntegratedGradients(model.drugLayer)

    print2log('Begin score calculation with integrated gradients for model %s'%i)
    # Per output latent variable input importance translation captum
    # 1st dimesion input
    # 2nd dimesion output
    scores = torch.zeros (model.drugLayer.mask.T.shape)
    for target in range(scores.shape[1]):
        attr, delta = ig.attribute(torch.eye(X.shape[1]).double(),target=target,n_steps=100,return_convergence_delta=True)
        scores[:,target] = torch.diagonal(attr, 0)

    interactions = pd.DataFrame(scores.numpy())
    interactions.index = drugInput.columns
    interactions.columns = drugTargets.columns
    interactions.to_csv(outPattern+'_interactionScores_'+str(i)+'.csv')
    # interactions = pd.read_csv(outPattern+'_interactionScores_'+str(i)+'.csv',index_col=0)

    scores = torch.tensor(interactions.values)
    df_mask = pd.DataFrame(mask.numpy())
    df_mask.index = drugInput.columns
    df_mask.columns = drugTargets.columns
    df_mask.reset_index(inplace=True)
    df_mask = df_mask.rename(columns={'index': 'drug'})
    df_mask = pd.melt(df_mask, id_vars=['drug'])
    df_mask['value'] = df_mask['value'].transform(lambda x: 'pos' if abs(x) > 0 else 'neg')
    df_mask = df_mask.rename(columns={'value': 'Prior knowledge'})

    thresholds = list(np.linspace(start=0, stop=60, num=1000))
    pvalues = []
    odds_ratios = []
    sensitivity = []
    specificity = []
    Gmean = []
    for th in thresholds:
        masked_scores = scores * (1. * torch.abs(scores) >= th)
        df_masked_scores = pd.DataFrame(masked_scores.numpy())
        df_masked_scores.index = drugInput.columns
        df_masked_scores.columns = drugTargets.columns
        df_masked_scores.reset_index(inplace=True)
        df_masked_scores = df_masked_scores.rename(columns={'index': 'drug'})
        df_masked_scores = pd.melt(df_masked_scores, id_vars=['drug'])
        df_masked_scores['value'] = df_masked_scores['value'].transform(lambda x: 'pos' if abs(x) > 0 else 'neg')
        df_masked_scores = df_masked_scores.rename(columns={'value': 'Inferred'})

        merged_df = pd.merge(df_mask, df_masked_scores, how='left', left_on=['drug', 'variable'],
                             right_on=['drug', 'variable'])
        contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])

        odds, p = scipy.stats.fisher_exact(contigency_table)
        pvalues.append(p)
        odds_ratios.append(odds)

        # Build also sensitivity, specificity curve and calculate G-mean max
        tp = contigency_table.loc['pos', 'pos']
        fn = contigency_table.loc['pos', 'neg']
        sens = tp / (tp + fn)
        tn = contigency_table.loc['neg', 'neg']
        fp = contigency_table.loc['neg', 'pos']
        spec = tn / (tn + fp)

        sensitivity.append(sens)
        specificity.append(spec)
        Gmean.append(np.sqrt(sens * spec))

    # optimalGIndex = np.argmax(Gmean)
    optimalGIndex = np.where(np.array(Gmean) == np.array(Gmean).max())[0][-1]
    # Plot contigency for ROC threshold
    df_mask = pd.DataFrame(mask.numpy())
    df_mask.index = drugInput.columns
    df_mask.columns = drugTargets.columns
    df_mask.reset_index(inplace=True)
    df_mask = df_mask.rename(columns={'index': 'drug'})
    df_mask = pd.melt(df_mask, id_vars=['drug'])
    df_mask['value'] = df_mask['value'].transform(lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
    df_mask = df_mask.rename(columns={'value': 'Prior knowledge'})
    masked_scores = scores * (1. * torch.abs(scores) >= np.max(thresholds[optimalGIndex]))
    df_masked_scores = pd.DataFrame(masked_scores.numpy())
    df_masked_scores.index = drugInput.columns
    df_masked_scores.columns = drugTargets.columns
    df_masked_scores.reset_index(inplace=True)
    df_masked_scores = df_masked_scores.rename(columns={'index': 'drug'})
    df_masked_scores = pd.melt(df_masked_scores, id_vars=['drug'])
    df_masked_scores['value'] = df_masked_scores['value'].transform(
        lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
    df_masked_scores = df_masked_scores.rename(columns={'value': 'Inferred'})
    merged_df = pd.merge(df_mask, df_masked_scores, how='left', left_on=['drug', 'variable'],
                         right_on=['drug', 'variable'])
    merged_df.to_csv(outPattern + '_mergedInteractions_ROC_'+str(i)+'.csv')


    # all_scores[i,:,:] = scores.numpy()
    all_scores[i, :, :] = interactions.values

# Calculate average score
scores = torch.tensor(np.mean(all_scores,0))
scores_std = torch.tensor(np.std(all_scores,0))
interactions = pd.DataFrame(scores.numpy())
interactions.index = drugInput.columns
interactions.columns = drugTargets.columns
interactions.to_csv(outPattern+'_interactionScoresEnsembled.csv')
interactions_std = pd.DataFrame(scores_std.numpy())
interactions_std.index = drugInput.columns
interactions_std.columns = drugTargets.columns
interactions_std.to_csv(outPattern+'_interactionScoresSTDEnsembled.csv')

mask = model.drugLayer.mask.T

### Run for different theresholds to consider positive interaction and calculate a p-value for the observed table
df_mask = pd.DataFrame(mask.numpy())
df_mask.index = drugInput.columns
df_mask.columns = drugTargets.columns
df_mask.reset_index(inplace=True)
df_mask = df_mask.rename(columns = {'index':'drug'})
df_mask = pd.melt(df_mask,id_vars=['drug'])
df_mask['value'] = df_mask['value'].transform(lambda x: 'pos' if abs(x)>0 else 'neg')
df_mask = df_mask.rename(columns = {'value':'Prior knowledge'})

thresholds = list(np.linspace(start=0, stop=60, num=1000))
pvalues = []
odds_ratios = []
sensitivity = []
specificity = []
Gmean = []
for th in thresholds:
    masked_scores = scores * (1.*torch.abs(scores)>=th)
    df_masked_scores = pd.DataFrame(masked_scores.numpy())
    df_masked_scores.index = drugInput.columns
    df_masked_scores.columns = drugTargets.columns
    df_masked_scores.reset_index(inplace=True)
    df_masked_scores = df_masked_scores.rename(columns = {'index':'drug'})
    df_masked_scores = pd.melt(df_masked_scores,id_vars=['drug'])
    df_masked_scores['value'] = df_masked_scores['value'].transform(lambda x: 'pos' if abs(x)>0 else 'neg')
    df_masked_scores = df_masked_scores.rename(columns = {'value':'Inferred'})

    merged_df = pd.merge(df_mask, df_masked_scores,  how='left', left_on=['drug','variable'], right_on = ['drug','variable'])
    contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])

    odds , p = scipy.stats.fisher_exact(contigency_table)
    pvalues.append(p)
    odds_ratios.append(odds)

    # Build also sensitivity, specificity curve and calculate G-mean max
    tp = contigency_table.loc['pos','pos']
    fn = contigency_table.loc['pos','neg']
    sens = tp/(tp+fn)
    tn = contigency_table.loc['neg', 'neg']
    fp = contigency_table.loc['neg', 'pos']
    spec = tn / (tn + fp)

    sensitivity.append(sens)
    specificity.append(spec)
    Gmean.append(np.sqrt(sens*spec))

# optimalGIndex = np.argmax(Gmean)
optimalGIndex = np.where(np.array(Gmean) == np.array(Gmean).max())[0][-1]
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.tight_layout(pad=1.75)
fig.suptitle('ROC Curve')
ax1.plot(1-np.array(specificity),sensitivity)
# ax1.scatter(1-np.array(specificity),sensitivity,s = 6)
ax1.plot([1-np.array(specificity)[optimalGIndex], 1-np.array(specificity)[optimalGIndex]], np.array([0, 1])*sensitivity[optimalGIndex], 'red', linestyle='--')
ax1.plot([0, 1-np.array(specificity)[optimalGIndex]], np.array([1, 1])*sensitivity[optimalGIndex], 'red', linestyle='--')
ax1.set_ylabel('sensitivity',fontsize = 11)
ax1.set_xlabel('1 - specificity',fontsize = 11)
ax2.plot(thresholds,Gmean)
# ax2.scatter(thresholds,Gmean,s = 6)
ax2.plot([0, thresholds[optimalGIndex]], np.array([1, 1])*Gmean[optimalGIndex], 'red', linestyle='--')
ax2.plot([thresholds[optimalGIndex], thresholds[optimalGIndex]], np.array([0, 1])*Gmean[optimalGIndex], 'red', linestyle='--')
ax2.set_ylabel('G-mean',fontsize = 11)
ax2.set_xlabel('absolute score threshold',fontsize = 11)
plt.savefig(outPattern+"_ROC_threshold_selection.png", bbox_inches='tight', dpi=600)

# Plot contigency for ROC threshold
df_mask = pd.DataFrame(mask.numpy())
df_mask.index = drugInput.columns
df_mask.columns = drugTargets.columns
df_mask.reset_index(inplace=True)
df_mask = df_mask.rename(columns = {'index':'drug'})
df_mask = pd.melt(df_mask,id_vars=['drug'])
df_mask['value'] = df_mask['value'].transform(lambda x: 'Interaction' if abs(x)>0 else 'No interaction')
df_mask = df_mask.rename(columns = {'value':'Prior knowledge'})

masked_scores = scores * (1.*torch.abs(scores)>=np.max(thresholds[optimalGIndex]))
df_masked_scores = pd.DataFrame(masked_scores.numpy())
df_masked_scores.index = drugInput.columns
df_masked_scores.columns = drugTargets.columns
df_masked_scores.reset_index(inplace=True)
df_masked_scores = df_masked_scores.rename(columns = {'index':'drug'})
df_masked_scores = pd.melt(df_masked_scores,id_vars=['drug'])
df_masked_scores['value'] = df_masked_scores['value'].transform(lambda x: 'Interaction' if abs(x)>0 else 'No interaction')
df_masked_scores = df_masked_scores.rename(columns = {'value':'Inferred'})
merged_df = pd.merge(df_mask, df_masked_scores,  how='left', left_on=['drug','variable'], right_on = ['drug','variable'])
merged_df.to_csv(outPattern+'_mergedInteractions_ROC.csv')
contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])
contigency_table = contigency_table.loc[['No interaction','Interaction'],:]
plt.clf()
plt.figure(figsize=(8, 6))
res = sns.heatmap(contigency_table, annot=True, fmt='.0f',
                  cmap="YlGnBu", vmin=0.0, vmax=50000.0)
plt.title('Interactions contigency table',fontsize=13)
plt.savefig(outPattern+"_contigency_ROC.png", bbox_inches='tight', dpi=600)