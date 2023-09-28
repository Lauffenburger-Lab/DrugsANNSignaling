# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 08:46:53 2022

@author: nmeim
"""

import torch
import pandas as pd
import numpy as np
import bionetwork2 as bionetwork
import torch.nn.functional as F
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import AllChem
from rdkit import Chem
from matplotlib import pyplot as plt
import seaborn as sns
sns.set()

#inputModel = 'FinalRotationProjectResults/A375_l1000_signalmodel_drop02_nomask_ep5kt120bs15.pt'
inputModel = 'FinalRotationProjectResults/A375_l1000_signalmodel_WsimECFP4_withRegA_v1_ep5kt120bs15.pt'
inputMask = None
#inputMask = 'FinalRotationProjectResults/A375_l1000_signalmodel_drop02_withRegA_ep5kt120bs15.pt'

model = torch.load(inputModel)
model.eval()

# drugInput = pd.read_csv('data/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
# drugSmiles = drugInput.columns.values
# drugTargets = pd.read_csv('data/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
# cellInput = pd.read_csv('data/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
# cell = 'A375'
# drugInput=drugInput.loc[cellInput[cellInput[cell]==1].index,:]
# ecfp4 = []
# for drug in drugSmiles:
#     if ((drug!='Ctrl_PBS') and (drug!='Ctrl_Cell')):
#         mol = Chem.MolFromSmiles(drug)
#         ecfp4.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
#     else:
#         ecfp4.append(np.zeros(1024))

# Xecfp4 = torch.Tensor(np.array(ecfp4).copy()).double()
# drugLayer = torch.mm(Xecfp4,Xecfp4.transpose(0,1))
# drugLayer = pd.read_csv('../similarity_matrix_lvl3_data.csv',index_col=0)
# drugTargets = drugTargets.loc[drugInput.columns.values,:]
# drugLayer = drugLayer.loc[drugInput.columns.values,drugInput.columns.values]
# drugLayer = torch.tensor(drugLayer.values.copy(), dtype=torch.double)
# W = drugLayer.detach()

#W = torch.mm(torch.mm( model.drugLayer.ecfp4_fingerprints, model.drugLayer.Wsim),
#              model.drugLayer.ecfp4_fingerprints.transpose(0, 1))

W = torch.mul(model.drugLayer.drugSim, model.drugLayer.Wsim)

W= W.detach()

if inputMask is not None:
    model2 = torch.load(inputMask)
    mask = model2.drugLayer.mask
    A = model.drugLayer.A
else:
    A = torch.mul( model.drugLayer.A, model.drugLayer.mask)
    mask = model.drugLayer.mask
A = A.detach()
mask = mask.detach()

T = torch.mm(A,W.transpose(0,1)).transpose(0,1)
#T = model.drugLayer(torch.eye(model.drugLayer.drugSim.shape[0]).double()).transpose(0,1)
T = T.detach()
plt.figure()
plt.hist(mask.flatten().numpy())
plt.title('Original known targets')
plt.ylabel('counts')
plt.show()
plt.figure()
plt.hist(T.flatten().numpy(),bins=25)
plt.title('Whole Drug Layer')
plt.ylabel('counts')
plt.xlabel('weight value')
plt.show()
# plt.figure()
# plt.hist((T/torch.abs(T).max()).flatten().numpy(),bins=20)
# plt.show()
# plt.figure()
# plt.hist((T/torch.abs(T).max())[mask.transpose(0,1)==1].flatten().numpy(),bins=20)
# plt.show()

Z = (T - T.mean()) / torch.sqrt(T.var())
plt.figure()
plt.hist(Z.flatten().numpy(),bins=25)
plt.axvline(-1,color='red',linestyle='dashed')
plt.axvline(1,color='red',linestyle='dashed')
plt.title('Whole Drug Layer scaled')
plt.ylabel('counts')
plt.xlabel('z-scored weight')
plt.show()

plt.figure()
plt.hist(Z[mask.transpose(0,1)==1].flatten().numpy(),bins=25)
plt.axvline(-1,color='red',linestyle='dashed')
plt.axvline(1,color='red',linestyle='dashed')
plt.title('Drug Layer scaled only for known interactions')
plt.ylabel('counts')
plt.xlabel('z-scored weight')
plt.show()

plt.figure()
plt.hist(1*(torch.abs(Z[mask.transpose(0,1)==1])>1).flatten().numpy(),bins=20)
plt.title('Drug Layer binary only for known interactions')
plt.ylabel('counts')
plt.show()

plt.figure()
plt.hist(W.flatten().numpy(),bins=20)
plt.title('Trainable drug similarity dot product matrix')
plt.ylabel('counts')
plt.xlabel('weight value')
plt.show()

plt.figure()
plt.hist(model.drugLayer.Wsim.detach().flatten().numpy(),bins=20)
plt.title('Trainable drug similarity weigths')
plt.ylabel('counts')
plt.xlabel('weight value')
plt.show()

plt.figure()
plt.hist(model.drugLayer.A.detach().flatten().numpy(),bins=20)
plt.title('Trainable drug-target matrix before masking')
plt.ylabel('counts')
plt.xlabel('weight value')
plt.show()