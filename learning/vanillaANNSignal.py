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

device = torch.device('cuda')

class DenseANN(torch.nn.Module):
    def __init__(self, in_channel, hidden_layers, outChannel,dropRate=0.1, activation=torch.nn.LeakyReLU(0.01), bias=True,dtype=torch.double):

        super(DenseANN, self).__init__()

        self.bias = bias
        self.num_hidden_layers = len(hidden_layers)
        self.bn = torch.nn.ModuleList()
        self.linear_layers = torch.nn.ModuleList()
        self.linear_layers.append(torch.nn.Linear(in_channel, hidden_layers[0], bias=bias,dtype = dtype))
        self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[0], momentum=0.6,dtype = dtype))
        for i in range(1, len(hidden_layers)):
            self.linear_layers.append(torch.nn.Linear(hidden_layers[i - 1], hidden_layers[i], bias=bias,dtype = dtype))
            self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[i], momentum=0.6,dtype = dtype))

        self.linear_latent = torch.nn.Linear(hidden_layers[-1],
                                             outChannel,
                                             bias=False,
                                             dtype = dtype)
        self.activation = activation
        self.dropout = torch.nn.Dropout(dropRate)
        self.drop_in = torch.nn.Dropout(0.3)

        self.init_emb()

    def init_emb(self):
        for m in self.modules():
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.fill_(0.0)

    def forward(self, x):
        x = self.drop_in(x)
        for i in range(self.num_hidden_layers):
            x = self.linear_layers[i](x)
            x = self.bn[i](x)
            x = self.activation(x)
            x = self.dropout(x)
        z_latent = self.linear_latent(x)
        return z_latent

    def L2Regularization(self, L2):

        weightLoss = 0.
        biasLoss = 0.
        for i in range(self.num_hidden_layers):
            weightLoss = weightLoss + L2 * torch.sum((self.linear_layers[i].weight)**2)
            if self.bias==True:
                biasLoss = biasLoss + L2 * torch.sum((self.linear_layers[i].bias)**2)
        weightLoss = weightLoss + L2 * torch.sum((self.linear_latent.weight)**2)
        L2Loss = biasLoss + weightLoss
        return(L2Loss)


cell_lines = ['A375','A549','HA1E','HCC515','HT29','PC3','MCF7','HEPG2']
#cell_lines = ['HEPG2']
parser = argparse.ArgumentParser(prog='Cell-line simulation')
# parser.add_argument('--leaveIn', action='store', default=None)
parser.add_argument('--pretrained', action='store', default=None)
parser.add_argument('--model_type', action='store', default=None)
parser.add_argument('--outPattern', action='store', default=None)
args = parser.parse_args()
# currentId = int(args.leaveIn)
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
maxIter = 200
noiseLevel = 10
L2 = 1e-6

for cell in cell_lines:

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
    model1 = torch.load('../results/VanillaANNSignal/models/'+pretrained_path)
    model1.eval()
    # if model_type==5:
    #     Xin = torch.mm(X,drugSim)
    #     Xin = model1.drugLayer(Xin)
    # else:
    Xin = model1.drugLayer(X)
    Xin = Xin.detach()

    Xin = Xin.to(device)
    Y = Y.to(device)
    ###Y = Y[torch.randperm(Y.size()[0]), :]

    model = DenseANN(Xin.shape[1],[120],Y.shape[1],).to(device)

    sampleName = drugInput.index.values
    criterion = torch.nn.MSELoss(reduction='mean')
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01,weight_decay=0) ###Changed to using RMSprop
    resetState = optimizer.state.copy()

    mLoss = criterion(torch.mean(Y, dim=0)*torch.ones(Y.shape).to(device), Y)
    print2log(mLoss)


    N = Xin.shape[0]

    e = 0
    trainLossMEAN = []
    trainLossSTD = []
    for e in range(e, maxIter):

        curLoss = []
        trainloader = bionetwork.getSamples(N, batchSize)
        for dataIndex in trainloader:
            model.train()
            optimizer.zero_grad()

            dataIn = Xin[dataIndex, :].view(len(dataIndex), Xin.shape[1])
            dataOut = Y[dataIndex, :].view(len(dataIndex), Y.shape[1])

            Yhat = model(dataIn)

            fitLoss = criterion(dataOut, Yhat)
            RegtLoss = model.L2Regularization(L2)

            loss = fitLoss + RegtLoss
            loss.backward()
            optimizer.step()

            curLoss.append(fitLoss.item())

        outString = 'Cell-line %s, Epoch=%s / %s'%(cell,e,maxIter)
        outString += ', loss={:.4f}'.format(loss.item())
        outString += ', regularization={:.4f}'.format(RegtLoss.item())
        trainLossMEAN.append(numpy.mean(curLoss))
        trainLossSTD.append(numpy.std(curLoss))
        if e % 50 == 0:
            print2log(outString)

        if numpy.logical_and(e % 100 == 0, e>0):
            optimizer.state = resetState.copy()
    print2log(outString)



    model.eval()
    Yhat = model(Xin)
    torch.save(model, '../results/vanillaANNSignal//train/'+cell+'_'+outPattern+'.pt')

    #%%
    plt.rcParams["figure.figsize"] = (9,6)
    plt.figure()
    T = numpy.array(range(len(trainLossMEAN)))
    plt.plot(T, trainLossMEAN)
    plt.xlim([0, len(T)])
    curColor = plt.gca().lines[-1].get_color()
    plt.fill_between(T,
                     numpy.array(trainLossMEAN) - numpy.array(trainLossSTD),
                     numpy.array(trainLossMEAN) + numpy.array(trainLossSTD),
                     color=curColor, alpha=0.2)
    plt.plot([0, len(T)], numpy.array([1, 1]) * mLoss.item(), 'black', linestyle='--')
    plt.ylabel('Loss')
    plt.xlabel('Epoch')
    plt.yscale('log')
    plt.savefig('../results/vanillaANNSignal//train/fig1_'+cell+'_'+outPattern+'.png',dpi=600)

    plt.figure()
    plt.scatter(Yhat.detach().cpu().numpy(), Y.detach().cpu().numpy(), alpha=0.2)
    plotting.lineOfIdentity()
    plotting.addCorrelation(Yhat.detach().cpu(), Y.detach().cpu())
    plt.xlabel('Fit')
    plt.ylabel('Data')
    plt.gca().axis('equal')
    plt.gca().set_xticks([0, 0.5, 1])
    plt.gca().set_yticks([0, 0.5, 1])
    plt.savefig('../results/vanillaANNSignal//train/fig4_'+cell+'_'+outPattern+'.png',dpi=600)

    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu(), Y.detach().cpu(), outNameGene)
    plt.savefig('../results/vanillaANNSignal//train/fig7_'+cell+'_'+outPattern+'.png',dpi=600)
    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu().T, Y.detach().cpu().T, sampleName)
    plt.savefig('../results/vanillaANNSignal//train/fig8_'+cell+'_'+outPattern+'.png',dpi=600)
    plotting.displayData(Y.detach().cpu(), sampleName, outNameGene)
    plt.savefig('../results/vanillaANNSignal//train/fig9_'+cell+'_'+outPattern+'.png',dpi=600)


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
    Y = torch.tensor(TFOutput.values, dtype=torch.double).to(device)

    model1.eval()
    Xin = model1.drugLayer(X)
    # Xin = model1.drugLayer(Xin)
    Xin = Xin.detach().to(device)

    model.eval()
    Yhat = model(Xin)

    plt.figure()
    plt.scatter(Yhat.detach().cpu().numpy(), Y.detach().cpu().numpy(), alpha=0.2)
    plotting.lineOfIdentity()
    plotting.addCorrelation(Yhat.detach().cpu(), Y.detach().cpu())
    plt.xlabel('Fit')
    plt.ylabel('Data')
    plt.gca().axis('equal')
    plt.gca().set_xticks([0, 0.5, 1])
    plt.gca().set_yticks([0, 0.5, 1])
    plt.savefig('../results/vanillaANNSignal//validation/performance_'+cell+'_'+outPattern+'.png',dpi=600)

    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu(), Y.detach().cpu(), outNameGene)
    plt.savefig('../results/vanillaANNSignal//validation/performancePerTF_'+cell+'_'+outPattern+'.png',dpi=600)
    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu().T, Y.detach().cpu().T, outconds)
    plt.savefig('../results/vanillaANNSignal//validation/performancePerSample_'+cell+'_'+outPattern+'.png',dpi=600)