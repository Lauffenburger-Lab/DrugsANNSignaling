import torch
import numpy
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

cell_lines = ['A375','A549','HA1E','HCC515','HT29','PC3','MCF7','HEPG2']
parser = argparse.ArgumentParser(prog='Cell-line simulation')
parser.add_argument('--leaveIn', action='store', default=None)
parser.add_argument('--pretrained', action='store', default=None)
parser.add_argument('--model_type', action='store', default=None)
parser.add_argument('--outPattern', action='store', default=None)
args = parser.parse_args()
currentId = int(args.leaveIn)
cell = cell_lines[currentId]
pretrained_path = args.pretrained
model_type = int(args.model_type)
outPattern = args.outPattern

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



#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pandas.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01)
spectralCapacity = numpy.exp(numpy.log(1e-2)/bionetParams['iterations'])


drugInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput = pandas.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
cellInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput=TFOutput.loc[cellInput[cellInput[cell]==1].index,:]
drugInput=drugInput.loc[cellInput[cellInput[cell]==1].index,:]
sigs_to_drop = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/smiles_and_sigs_to_drop.csv', index_col=0)
TFOutput=TFOutput.loc[~TFOutput.index.isin(sigs_to_drop.sig_id),:]
drugInput=drugInput.loc[~drugInput.index.isin(sigs_to_drop.sig_id),:]

#Subset input and output to intersecting nodes
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
drugSim = pandas.read_csv('../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv',index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

sampleName = drugInput.index.values
model1 = torch.load('../results/models/'+pretrained_path)
model1.eval()

Xin = model1.drugLayer(X)
Xin = Xin.detach()

sampleName = drugInput.index.values
model = SingleBionet.model(networkList, nodeNames,modeOfAction, inputAmplitude, projectionAmplitude, inName, outName, bionetParams, torch.double)
model.inputLayer.weights.requires_grad = False
model.network.preScaleWeights()
#%%
criterion = torch.nn.MSELoss(reduction='mean')
optimizer = torch.optim.Adam(model.parameters(), lr=1,weight_decay=0) ###Changed to using RMSprop
resetState = optimizer.state.copy()

mLoss = criterion(torch.mean(Y, dim=0)*torch.ones(Y.shape), Y)
print2log(mLoss)


stats = plotting.initProgressObject(maxIter)
N = Xin.shape[0]
curState = torch.rand((X.shape[0], model.network.bias.shape[0]), dtype=torch.double, requires_grad=False) # model.network.bias.shape[0]

e = 0
for e in range(e, maxIter):
    curLr = bionetwork.oneCycle(e, maxIter, maxHeight = 2e-3, minHeight = 1e-8, peak = 1000)
    optimizer.param_groups[0]['lr'] = curLr

    curLoss = []
    curEig = []
    trainloader = bionetwork.getSamples(N, batchSize)
    for dataIndex in trainloader:
        model.train()
        model.network.weights.data = model.network.weights.data + 1e-8 * torch.randn(model.network.weights.shape) #breaks potential symmetries
        optimizer.zero_grad()

        dataIn = Xin[dataIndex, :].view(len(dataIndex), Xin.shape[1])
        dataOut = Y[dataIndex, :].view(len(dataIndex), Y.shape[1])

        Yin = model.inputLayer(dataIn)
        Yin = Yin + noiseLevel * curLr * torch.randn(Yin.shape)
        YhatFull = model.network(Yin)
        Yhat = model.projectionLayer(YhatFull)

        curState[dataIndex, :] = YhatFull.detach()

        fitLoss = criterion(dataOut, Yhat)

        signConstraint = MoAFactor * torch.sum(torch.abs(model.network.weights[model.network.getViolations(model.network.weights)]))
        stateLoss = 1e-5 * bionetwork.uniformLoss(curState, dataIndex, YhatFull, maxConstraintFactor = 50)
        biasLoss = L2 * torch.sum(torch.square(model.network.bias))
        weightLoss = L2 * (torch.sum(torch.square(model.network.weights)) + torch.sum(1/(torch.square(model.network.weights) + 0.5)))
        projectionLoss = 1e-6 * torch.sum(torch.square(model.projectionLayer.weights - projectionAmplitude))
        spectralRadiusLoss, spectralRadius = bionetwork.spectralLoss(model, YhatFull, model.network.weights, expFactor = 10)

        loss = fitLoss + signConstraint + biasLoss + weightLoss + spectralFactor * spectralRadiusLoss + stateLoss + projectionLoss
        loss.backward()
        optimizer.step()

        curEig.append(spectralRadius.item())
        curLoss.append(fitLoss.item())

        stats = plotting.storeProgress(stats, e, curLoss, curEig, curLr, violations=torch.sum(model.network.getViolations(model.network.weights)).item())

    if e % 50 == 0:
        plotting.printStats(e, stats)

    if numpy.logical_and(e % 200 == 0, e>0):
        optimizer.state = resetState.copy()



stats = plotting.finishProgress(stats)
model.eval()
Yhat, YhatFull = model(Xin)
torch.save(model, '../results/models/'+cell+'_'+outPattern+'.pt')

#%%
plt.rcParams["figure.figsize"] = (9,6)
plt.figure()


T = numpy.array(range(stats['loss'].shape[0]))
plotting.shadePlot(T, plotting.movingaverage(stats['loss'], 5), plotting.movingaverage(stats['lossSTD'], 10))
plt.plot([0, len(T)], numpy.array([1, 1])*mLoss.item(), 'black', linestyle='--')
plt.xlim([0, len(T)])
plt.ylim(bottom=1e-3)
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.yscale('log')
plt.savefig('../results/train/fig1_'+cell+'_'+outPattern+'.png',dpi=600)
#
plt.figure()
plt.plot(T, stats['rate'])
plt.plot([0, maxIter], [0, 0], 'black')
plt.legend(numpy.array(['lr', 'loss adjusted lr']), frameon=False)
plt.ylabel('Learning rate')
plt.xlabel('Epoch')
plt.savefig('../results/train/fig2_'+cell+'_'+outPattern+'.png',dpi=600)
#
plt.figure()
plt.plot([0, maxIter], [1, 1], 'black')
plt.plot([0, len(T)], spectralCapacity * numpy.array([1, 1]), 'black', linestyle='--')
plotting.shadePlot(T, plotting.movingaverage(stats['eig'], 5), plotting.movingaverage(stats['eigSTD'], 5))
plt.ylabel('Spectral radius')
plt.xlabel('Epoch')
plt.savefig('../results/train/fig3_'+cell+'_'+outPattern+'.png',dpi=600)

plt.figure()
plt.scatter(Yhat.detach().numpy(), Y.detach().numpy(), alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(Yhat, Y)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
plt.savefig('../results/train/fig4_'+cell+'_'+outPattern+'.png',dpi=600)

plt.figure()
srModel = bionetwork.getAllSpectralRadius(model, YhatFull)

plt.hist(srModel)
plt.ylabel('SR model')
plt.savefig('../results/train/fig6_'+cell+'_'+outPattern+'.png',dpi=600)
plt.rcParams["figure.figsize"] = (12,10)
plt.figure()
rank = plotting.compareAllTFs(Yhat, Y, outNameGene)
plt.savefig('../results/train/fig7_'+cell+'_'+outPattern+'.png',dpi=600)
plt.figure()
rank = plotting.compareAllTFs(Yhat.T, Y.T, sampleName)
plt.savefig('../results/train/fig8_'+cell+'_'+outPattern+'.png',dpi=600)
plotting.displayData(Y, sampleName, outNameGene)
plt.savefig('../results/train/fig9_'+cell+'_'+outPattern+'.png',dpi=600)


### Evaluate validation results
drugInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
drugSmiles = drugInput.columns.values
drugTargets = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput = pandas.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
cellInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
TFOutput=TFOutput.loc[cellInput[cellInput[cell]==1].index,:]
drugInput=drugInput.loc[cellInput[cellInput[cell]==1].index,:]
sigs_to_drop = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/smiles_and_sigs_to_drop.csv', index_col=0) # PROSEXE ALLO GIA TA KAINOURIA
TFOutput=TFOutput.loc[TFOutput.index.isin(sigs_to_drop.sig_id),:]
drugInput=drugInput.loc[drugInput.index.isin(sigs_to_drop.sig_id),:]

#Subset input and output to intersecting nodes
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
drugSim = pandas.read_csv('../preprocessing/preprocessed_data/ChemicalSims/out_lvl3_similaritiess.csv',index_col=0)
drugSim = drugSim.loc[drugInput.columns.values,drugInput.columns.values]
drugSim = torch.tensor(drugSim.values.copy(), dtype=torch.double)

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)

model1.eval()
Xin = model1.drugLayer(X)
# Xin = model1.drugLayer(Xin)
Xin = Xin.detach()

model.eval()
Yhat, YhatFull = model(Xin)

plt.figure()
plt.scatter(Yhat.detach().numpy(), Y.detach().numpy(), alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(Yhat, Y)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
plt.savefig('../results/test/performance_'+cell+'_'+outPattern+'.png',dpi=600)
