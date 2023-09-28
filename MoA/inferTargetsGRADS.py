# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 05:52:42 2022

@author: nmeim
"""

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
parser.add_argument('--pretrained', action='store', default=None)
parser.add_argument('--outPattern', action='store', default=None)
parser.add_argument('--cell_line', action='store', default=None)
args = parser.parse_args()
inputModel = args.pretrained
outPattern = args.outPattern
cell = args.cell_line

### Load network
#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('data/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('data/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for a375 was 120
spectralCapacity = np.exp(np.log(1e-2)/bionetParams['iterations'])
### Load the data
ConvertToEmpProb = False
drugInput = pd.read_csv('data/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pd.read_csv('data/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput = pd.read_csv('data/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
cellInput = pd.read_csv('data/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
#cell = 'A375'
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

model = torch.load('CVL1000_Paper/'+inputModel)
model.eval()
Yhat, YhatFull = model(X)
ig = IntegratedGradients(model.drugLayer)

print2log('Begin score calculation with integrated gradients')
# Per output latent variable input importance translation captum
# 1st dimesion input
# 2nd dimesion output
scores = torch.zeros (model.drugLayer.mask.T.shape)
for target in range(scores.shape[1]):
    #model.zero_grad()
    attr, delta = ig.attribute(torch.eye(X.shape[1]).double(),target=target,n_steps=100,return_convergence_delta=True)
    scores[:,target] = torch.diagonal(attr, 0)
# W = torch.mul(model.drugLayer.drugSim, model.drugLayer.Wsim)
# A = torch.mul(model.drugLayer.A,model.drugLayer.mask)
# kappa =  model.drugLayer.bn.weight /torch.sqrt(model.drugLayer.bn.running_var +  model.drugLayer.bn.eps)
# K = kappa
# scores = torch.mm(A,(kappa * W).T)
# scores = scores.T.detach()
#print2log(scores)
#print2log(scores.shape)
plt.figure()
plt.hist(scores.flatten().numpy())
plt.title('Gradient scores distribution')
plt.ylabel('counts')
plt.xlabel('gradient score')
plt.savefig('../'+outPattern+'_interactionScores.png', bbox_inches='tight',dpi=600)

interactions = pd.DataFrame(scores.numpy())
interactions.index = drugInput.columns
interactions.columns = drugTargets.columns
interactions.to_csv('../'+outPattern+'_interactionScores.csv')

### Create heatmap of ranks and gradients scores
nCuts = 200
top = scores.shape[1]
topScores , inds = torch.topk(torch.abs(scores),top,dim=1)
topScores = topScores.detach()
inds = inds.detach()
topScores = topScores.detach()
topScores = pd.DataFrame(topScores.detach().numpy())
topScores = pd.melt(topScores)
topScores.columns = ['rank','score']
topScores['rank'] = topScores['rank'] + 1
cuts = pd.DataFrame({'score bins' : pd.cut(topScores['score'], nCuts)})
topScores = pd.concat((topScores,cuts),1)
topScores['value'] = 1
binnedHeat = topScores.iloc[:,[0,2,3]].groupby(['score bins','rank']).sum()
binnedHeat = binnedHeat.unstack(level = 0)
#binnedHeat
plt.figure(figsize=(12,12))
plt.clf()
sns.heatmap(binnedHeat['value']) # mask=binnedHeat['value']==0
plt.title('Distribution of scores of jth hit')
plt.tight_layout()
plt.savefig("../"+outPattern+"_binnedHeatmapAbs.png", bbox_inches='tight', dpi=600)

df = topScores.loc[:,['score bins','rank']]
df = df.groupby("score bins").nunique().cumsum()
df.reset_index(inplace=True)
df['change'] = 1
for i in range(1,len(df)):
    df.iloc[i,2] = (df.iloc[i,1] - df.iloc[i-1,1])/df.iloc[i-1,1]
df = df.iloc[1:,[0,2]]
df.columns = ["score bins","% change in cumulative sum of ranks"]
plt.figure(figsize=(9,12))
plt.clf()
fig = sns.pointplot(data=df, x="score bins", y="% change in cumulative sum of ranks")
plt.xticks(list(range(0,nCuts,10)))
plt.setp(fig.get_xticklabels(), rotation=90)
plt.tight_layout()
plt.savefig("../"+outPattern+"_changeInCumSumRankScoreBinnedAbs.png", bbox_inches='tight', dpi=600)
ind = np.where(df["% change in cumulative sum of ranks"]==0)[0][0]
#print2log(df.iloc[ind,0])

# =============================================================================
# ## Change in mean rank per bin
# df = topScores.loc[:,['score','rank']].drop_duplicates()
# plt.figure(figsize=(9,12))
# plt.clf()
# fig = sns.scatterplot(data=df, x="score", y="rank")
# #plt.xticks(list(range(0,nCuts,10)))
# #plt.setp(fig.get_xticklabels(), rotation=90)
# plt.tight_layout()
# plt.savefig("../"+outPattern+"_RankScoreAbs.png", bbox_inches='tight', dpi=600)
# df = df.groupby('score').mean()
# df.reset_index(inplace=True)
# df = df.dropna()
# df =df.sort_values('score') # only for when using only score
# df['change'] = 1
# for i in range(1,len(df)):
#     df.iloc[i,2] = (df.iloc[i,1] - df.iloc[i-1,1])/df.iloc[i-1,1]
# df = df.iloc[1:,[0,2]]
# df.columns = ["score","% change in mean rank"]
# plt.figure(figsize=(9,12))
# plt.clf()
# fig = sns.scatterplot(data=df, x="score", y="% change in mean rank")
# #plt.xticks(list(range(0,nCuts,10)))
# #plt.setp(fig.get_xticklabels(), rotation=90)
# plt.tight_layout()
# plt.savefig("../"+outPattern+"_changeInRankScoreAbs.png", bbox_inches='tight', dpi=600)
# =============================================================================

## Non-absolute score
topScores , inds = torch.topk(torch.abs(scores),top,dim=1)
topScores = topScores.detach()
inds = inds.detach()
topScores = torch.cat((scores[list(range(scores.shape[0])),inds[:,0]].view(scores.shape[0],1),
                       scores[list(range(scores.shape[0])),inds[:,1]].view(scores.shape[0],1)),1)
for i in range(2,scores.shape[1]):
    topScores = torch.cat((topScores,scores[list(range(scores.shape[0])),inds[:,i]].view(scores.shape[0],1)),1)
topScores = topScores.detach()
topScores = pd.DataFrame(topScores.detach().numpy())
topScores = pd.melt(topScores)
topScores.columns = ['rank','score']
topScores['rank'] = topScores['rank'] + 1
cuts = pd.DataFrame({'score bins' : pd.cut(topScores['score'], nCuts)})
topScores = pd.concat((topScores,cuts),1)
topScores['value'] = 1
binnedHeat = topScores.iloc[:,[0,2,3]].groupby(['score bins','rank']).sum()
binnedHeat = binnedHeat.unstack(level = 0)
#binnedHeat
plt.figure(figsize=(12,12))
plt.clf()
sns.heatmap(binnedHeat['value']) # mask=binnedHeat['value']==0
plt.title('Distribution of scores of jth hit')
plt.tight_layout()
plt.savefig("../"+outPattern+"_binnedHeatmap.png", bbox_inches='tight', dpi=600)

mask = model.drugLayer.mask.T

### Threshold using k-means approach
df_mask = pd.DataFrame(mask.numpy())
df_mask.index = drugInput.columns
df_mask.columns = drugTargets.columns
df_mask.reset_index(inplace=True)
df_mask = df_mask.rename(columns = {'index':'drug'})
df_mask = pd.melt(df_mask,id_vars=['drug'])
df_mask['value'] = df_mask['value'].transform(lambda x: 'pos' if abs(x)>0 else 'neg')
df_mask = df_mask.rename(columns = {'value':'Prior knowledge'})
df_mask['Prior knowledge'] = df_mask['Prior knowledge'].transform(lambda x: 'Interaction' if x=='pos' else 'No interaction')

thresholds = []
for drug in interactions.index:
    ### Adopted from https://github.com/SBRG/pymodulon/blob/c333004a421f1f2f7ab32f5f4b29b6b2a708bb41/src/pymodulon/core.py#L1093
    data = interactions.loc[drug,:]
    model = KMeans(n_clusters=3, random_state=1)
    model.fit(abs(data).values.reshape(-1, 1))
    df = pd.DataFrame(abs(data))
    df["cluster"] = model.labels_
    # Get top two clusters
    counts = df.cluster.value_counts().sort_values(ascending=True)
    idx1 = counts.index[1]
    idx2 = counts.index[2]
    # idx1 = counts.index[0]
    # idx2 = counts.index[2]
    clust1 = df[df.cluster == idx1]
    clust2 = df[df.cluster == idx2]
    # Get midpoint between lowest interaction group and maximum insignificant
    if clust1[drug].min() >= clust2[drug].max():
        threshold = np.mean([clust1[drug].min(), clust2[drug].max()])
    else:
        threshold = np.mean([clust1[drug].max(), clust2[drug].mean()])
    thresholds.append(threshold)
plt.figure()
plt.hist(thresholds,100)
plt.title('Absolute gradient scores thresholds derived from a k-mean approach for each drug')
plt.ylabel('counts')
plt.xlabel('threshold')
plt.savefig('../'+outPattern+'_kmeans_thresholds.png', bbox_inches='tight',dpi=600)
print2log('Kmeans threshold %s'%np.mean(thresholds))

masked_scores = scores * (1.*torch.abs(scores)>=np.mean(thresholds))
df_masked_scores = pd.DataFrame(masked_scores.numpy())
df_masked_scores.index = drugInput.columns
df_masked_scores.columns = drugTargets.columns
df_masked_scores.reset_index(inplace=True)
df_masked_scores = df_masked_scores.rename(columns = {'index':'drug'})
df_masked_scores = pd.melt(df_masked_scores,id_vars=['drug'])
df_masked_scores['value'] = df_masked_scores['value'].transform(lambda x: 'Interaction' if abs(x)>0 else 'No interaction')
df_masked_scores = df_masked_scores.rename(columns = {'value':'Inferred'})
merged_df = pd.merge(df_mask, df_masked_scores,  how='left', left_on=['drug','variable'], right_on = ['drug','variable'])
merged_df.to_csv('../'+outPattern+'_mergedInteractions_kmeans.csv')
contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])
contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])
contigency_table = contigency_table.loc[['No interaction','Interaction'],:]
#print2log(contigency_table)
plt.clf()
plt.figure(figsize=(8, 6))
res = sns.heatmap(contigency_table, annot=True, fmt='.0f',
                  cmap="YlGnBu", vmin=0.0, vmax=50000.0)
plt.title('Interactions contigency table',fontsize=13)
plt.savefig("../"+outPattern+"_contigency_kmeans.png", bbox_inches='tight', dpi=600)

### Run for different theresholds to consider positive interaction and calculate a p-value for the observed table
df_mask = pd.DataFrame(mask.numpy())
df_mask.index = drugInput.columns
df_mask.columns = drugTargets.columns
df_mask.reset_index(inplace=True)
df_mask = df_mask.rename(columns = {'index':'drug'})
df_mask = pd.melt(df_mask,id_vars=['drug'])
df_mask['value'] = df_mask['value'].transform(lambda x: 'pos' if abs(x)>0 else 'neg')
df_mask = df_mask.rename(columns = {'value':'Prior knowledge'})

thresholds = list(np.linspace(start=0, stop=1.5, num=32))
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
    merged_df.to_csv('../'+outPattern+'_mergedInteractions.csv')
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

#plt.rcParams['figure.figsize'] = [10,6]
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.tight_layout(pad=1.75)
fig.suptitle('Fisher`s exact test for contegency tables for different score thresholds')
fig.text(0.5, 0.0, 'absolute score threshold', ha='center',fontsize = 13)
plt.setp((ax1, ax2), xticks=[0, 0.5, 1, 1.5])
ax1.plot(thresholds,pvalues,linestyle='dotted')
ax1.scatter(thresholds,pvalues,s = 12)
ax1.set_yscale('log')
ax1.plot([0, 1.5], np.array([1, 1])*0.05, 'red', linestyle='--')
ax1.set_ylabel('log(p.value)',fontsize = 13)
ax2.plot(thresholds,odds_ratios,linestyle='dotted')
ax2.scatter(thresholds,odds_ratios,s = 12)
ax2.set_ylabel('odds ratio',fontsize = 13)
plt.savefig("../"+outPattern+"_threshold_selection.png", bbox_inches='tight', dpi=600)

optimalGIndex = np.argmax(Gmean)
fig, (ax1, ax2) = plt.subplots(1, 2)
fig.tight_layout(pad=1.75)
fig.suptitle('ROC Curve')
plt.setp((ax1, ax2), xticks=[0, 0.5, 1, 1.5])
ax1.plot(1-np.array(specificity),sensitivity)
ax1.scatter(1-np.array(specificity),sensitivity,s = 8)
ax1.plot([1-np.array(specificity)[optimalGIndex], 1-np.array(specificity)[optimalGIndex]], np.array([0, 1])*sensitivity[optimalGIndex], 'red', linestyle='--')
ax1.plot([0, 1-np.array(specificity)[optimalGIndex]], np.array([1, 1])*sensitivity[optimalGIndex], 'red', linestyle='--')
ax1.set_ylabel('sensitivity',fontsize = 13)
ax1.set_xlabel('1 - specificity',fontsize = 13)
ax2.plot(thresholds,Gmean)
ax2.scatter(thresholds,Gmean,s = 8)
ax2.plot([0, thresholds[optimalGIndex]], np.array([1, 1])*Gmean[optimalGIndex], 'red', linestyle='--')
ax2.plot([thresholds[optimalGIndex], thresholds[optimalGIndex]], np.array([0, 1])*Gmean[optimalGIndex], 'red', linestyle='--')
ax2.set_ylabel('G-mean',fontsize = 13)
ax2.set_xlabel('absolute score threshold',fontsize = 13)
plt.savefig("../"+outPattern+"_ROC_threshold_selection.png", bbox_inches='tight', dpi=600)



#scores_scaled = (scores - torch.mean(scores,0))/torch.std(scores,0)
#scores_scaled = (scores - scores.mean())/scores.std()
#masked_scores = scores * (1.*torch.abs(scores)>1.0) # soecific threshold of one (this how the trainbales Wsim and Amasked are produced now)
masked_scores = scores * (1.*torch.abs(scores)>=np.max(thresholds[np.argmin(pvalues)]))
#masked_scores = scores * (1.*torch.abs(scores_scaled)>2.0)
#new_mask = mask * masked_scores

df_mask['Prior knowledge'] = df_mask['Prior knowledge'].transform(lambda x: 'Interaction' if x=='pos' else 'No interaction')

df_masked_scores = pd.DataFrame(masked_scores.numpy())
df_masked_scores.index = drugInput.columns
df_masked_scores.columns = drugTargets.columns
df_masked_scores.reset_index(inplace=True)
df_masked_scores = df_masked_scores.rename(columns = {'index':'drug'})
df_masked_scores = pd.melt(df_masked_scores,id_vars=['drug'])
df_masked_scores['value'] = df_masked_scores['value'].transform(lambda x: 'Interaction' if abs(x)>0 else 'No interaction')
df_masked_scores = df_masked_scores.rename(columns = {'value':'Inferred'})

merged_df = pd.merge(df_mask, df_masked_scores,  how='left', left_on=['drug','variable'], right_on = ['drug','variable'])
merged_df.to_csv('../'+outPattern+'_mergedInteractions.csv')
contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])
contigency_table = contigency_table.loc[['No interaction','Interaction'],:]
#print2log(contigency_table)

plt.clf()
plt.figure(figsize=(8, 6))
res = sns.heatmap(contigency_table, annot=True, fmt='.0f',
                  cmap="YlGnBu", vmin=0.0, vmax=50000.0)
plt.title('Interactions contigency table',fontsize=13)
plt.savefig("../"+outPattern+"_contigency.png", bbox_inches='tight', dpi=600)

# Plot contigency for ROC threshold
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
merged_df.to_csv('../'+outPattern+'_mergedInteractions_ROC.csv')
contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])
contigency_table = contigency_table.loc[['No interaction','Interaction'],:]
plt.clf()
plt.figure(figsize=(8, 6))
res = sns.heatmap(contigency_table, annot=True, fmt='.0f',
                  cmap="YlGnBu", vmin=0.0, vmax=50000.0)
plt.title('Interactions contigency table',fontsize=13)
plt.savefig("../"+outPattern+"_contigency_ROC.png", bbox_inches='tight', dpi=600)

# ## Plot contigency for different fixed thresholds
# fixed_thresh = [0.2,0.3,0.5,0,75,1,1.1,1.25,1.5]
# for th in fixed_thresh:
#     masked_scores = scores * (1. * torch.abs(scores) >= th)
#     df_masked_scores = pd.DataFrame(masked_scores.numpy())
#     df_masked_scores.index = drugInput.columns
#     df_masked_scores.columns = drugTargets.columns
#     df_masked_scores.reset_index(inplace=True)
#     df_masked_scores = df_masked_scores.rename(columns={'index': 'drug'})
#     df_masked_scores = pd.melt(df_masked_scores, id_vars=['drug'])
#     df_masked_scores['value'] = df_masked_scores['value'].transform(lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
#     df_masked_scores = df_masked_scores.rename(columns={'value': 'Inferred'})
#     merged_df = pd.merge(df_mask, df_masked_scores, how='left', left_on=['drug', 'variable'],right_on=['drug', 'variable'])
#     merged_df.to_csv('../' + outPattern + '_mergedInteractions_th'+str(th)+'.csv')
#     contigency_table = pd.crosstab(index=merged_df['Prior knowledge'], columns=merged_df['Inferred'])
#     contigency_table = contigency_table.loc[['No interaction', 'Interaction'], :]
#     plt.clf()
#     plt.figure(figsize=(8, 6))
#     res = sns.heatmap(contigency_table, annot=True, fmt='.0f',
#                       cmap="YlGnBu", vmin=0.0, vmax=50000.0)
#     plt.title('Interactions contigency table', fontsize=13)
#     plt.savefig("../" + outPattern + "_contigency_th"+str(th)+".png", bbox_inches='tight', dpi=600)