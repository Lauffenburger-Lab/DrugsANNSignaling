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
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor as KNN
from sklearn.multioutput import MultiOutputRegressor
#from MSVR import MSVR
import features
import graph_layers
from features import one_of_k_encoding, one_of_k_encoding_unk, atom_features, bond_features, num_atom_features, num_bond_features, padaxis, tensorise_smiles #, concat_mol_tensors
from graph_layers import NeuralGraphHidden
import argparse
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

device = torch.device('cuda')

class DenseANN(torch.nn.Module):
    def __init__(self, in_channel, hidden_layers, outChannel,dropRate=0.15, drop_in=0,activation=torch.nn.ReLU(), bias=True,dtype=torch.double):

        super(DenseANN, self).__init__()

        self.bias = bias
        self.num_hidden_layers = len(hidden_layers)
        self.bn = torch.nn.ModuleList()
        self.linear_layers = torch.nn.ModuleList()
        self.linear_layers.append(torch.nn.Linear(in_channel, hidden_layers[0], bias=bias,dtype = dtype))
        self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[0], momentum=0.2,dtype = dtype))
        for i in range(1, len(hidden_layers)):
            self.linear_layers.append(torch.nn.Linear(hidden_layers[i - 1], hidden_layers[i], bias=bias,dtype = dtype))
            self.bn.append(torch.nn.BatchNorm1d(num_features=hidden_layers[i], momentum=0.2,dtype = dtype))

        self.linear_out = torch.nn.Linear(hidden_layers[-1],
                                             outChannel,
                                             bias=False,
                                             dtype = dtype)
        self.activation = activation
        self.dropout = torch.nn.Dropout(dropRate)
        self.drop_in_rate = drop_in
        if drop_in>0:
            self.drop_in = torch.nn.Dropout(drop_in)

        self.init_emb()

    def init_emb(self):
        for m in self.modules():
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.fill_(0.0)

    def forward(self, x):
        if self.drop_in_rate>0:
            x = self.drop_in(x)
        for i in range(self.num_hidden_layers):
            x = self.linear_layers[i](x)
            x = self.bn[i](x)
            x = self.activation(x)
            x = self.dropout(x)
        y = self.linear_out(x)
        y = self.activation(y)
        return y

    def L2Regularization(self, L2):

        weightLoss = 0.
        biasLoss = 0.
        for i in range(self.num_hidden_layers):
            weightLoss = weightLoss + L2 * torch.sum((self.linear_layers[i].weight)**2)
            if self.bias==True:
                biasLoss = biasLoss + L2 * torch.sum((self.linear_layers[i].bias)**2)
        weightLoss = weightLoss + L2 * torch.sum((self.linear_out.weight)**2)
        L2Loss = biasLoss + weightLoss
        return(L2Loss)

class GCNN(torch.nn.Module):
    def __init__(self,params):
        super(GCNN, self).__init__()
        self.params = params
        self.gcnn_layers = torch.nn.ModuleList()
        self.bn_layers = torch.nn.ModuleList()
        self.activations = torch.nn.ModuleList()
        self.dropouts = torch.nn.ModuleList()
        self.gcnn_layers.append(
            NeuralGraphHidden(params["num_atom_features"], params["num_bond_features"], params["graph_conv_width"],
                              params["max_degree"], activ=None, bias=True))
        self.bn_layers.append(torch.nn.BatchNorm1d(num_features=params["max_atoms"], momentum=0.2))
        self.activations.append(torch.nn.ReLU())
        self.dropouts.append(torch.nn.Dropout(params["dropout_graph_encoder"]))
        for i in range(1,params["number_of_gcnn_layers"]):
            self.gcnn_layers.append(NeuralGraphHidden( params["graph_conv_width"],params["num_bond_features"],params["graph_conv_width"],params["max_degree"] , activ = None, bias = True))
            self.bn_layers.append(torch.nn.BatchNorm1d(num_features=params["max_atoms"],momentum=0.2))
            self.activations.append(torch.nn.ReLU())
            self.dropouts.append(torch.nn.Dropout(params["dropout_graph_encoder"]))

        # self.conv1d_layers = torch.nn.ModuleList()
        # self.bn_layers_conv1d = torch.nn.ModuleList()
        # self.activations_conv1d = torch.nn.ModuleList()
        # self.dropout_conv1d = torch.nn.ModuleList()
        # self.num_cons = params["number_of_hidden_convolutions"]
        # tmp_conv_size = params["conv1d_in"]
        # if params["number_of_hidden_convolutions"]>0:
        #     self.conv1d_layers.append(torch.nn.Conv1d(params["conv1d_in"], params["conv1d_in"]//2, params["kernel_size"],bias=False))
        #     self.bn_layers_conv1d.append(torch.nn.BatchNorm1d(num_features=params["conv1d_in"]//2, momentum=0.2))
        #     self.activations_conv1d.append(torch.nn.ReLU())
        #     self.dropout_conv1d.append(torch.nn.Dropout(params["dropout_graph_encoder"]))
        #     tmp_conv_size =  params["conv1d_in"]//2
        # for i in range(1,params["number_of_hidden_convolutions"]):
        #     self.conv1d_layers.append(torch.nn.Conv1d(tmp_conv_size, tmp_conv_size//2, params["kernel_size"],bias=False))
        #     self.bn_layers_conv1d.append(torch.nn.BatchNorm1d(num_features=tmp_conv_size//2, momentum=0.2))
        #     self.activations_conv1d.append(torch.nn.ReLU())
        #     self.dropout_conv1d.append(torch.nn.Dropout(params["dropout_graph_encoder"]))
        #     tmp_conv_size = tmp_conv_size//2
        #
        # self.con1d_final = torch.nn.Conv1d(tmp_conv_size, 1, params["kernel_size"],bias=False)
        # self.act = torch.nn.ReLU()
        # self.dropout_conv1d.append(torch.nn.Dropout(params["dropout_graph_encoder"]))
        # self.ind_final_drop = len(self.dropout_conv1d)

        self.fc = torch.nn.Linear(params["graph_conv_width"],params["output_dim"],bias=True)
        self.fc_act = torch.nn.Hardtanh(min_val= 0.0, max_val=1.0)

        self.init_emb()

    def init_emb(self):
        for m in self.modules():
            if isinstance(m, torch.nn.Linear):
                torch.nn.init.xavier_uniform_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.fill_(0.0)
            if isinstance(m,torch.nn.Conv1d):
                torch.nn.init.xavier_normal_(m.weight.data)

    def forward(self,atoms,bonds,edges):
        #print2log(atoms.shape)
        x = self.gcnn_layers[0](atoms, bonds, edges)
        x = self.bn_layers[0](x)
        x = self.activations[0](x)
        x = self.dropouts[0](x)
        #print2log(x.shape)
        for i in range(1,self.params["number_of_gcnn_layers"]):
            x = self.gcnn_layers[i](x,bonds,edges)
            x = self.bn_layers[i](x)
            x = self.activations[i](x)
            x = self.dropouts[i](x)

        # for i in range(self.params["number_of_hidden_convolutions"]):
        #     x = self.conv1d_layers[i](x)
        #     x = self.bn_layers_conv1d[i](x)
        #     x = self.activations_conv1d[i](x)
        #     x = self.dropout_conv1d[i](x)

        # x = self.con1d_final(x)
        # x = self.act(x)
        # x = x.squeeze()
        # x = self.dropout_conv1d[-1](x)
        x = torch.mean(x,1)

        #print2log(x.shape)
        x = self.fc(x)
        x = self.fc_act(x)

        return x

    def L2Regularization(self, L2):
        weightLoss = 0.
        biasLoss = 0.
        for m in self.modules():
            if isinstance(m, torch.nn.Linear):
                weightLoss = weightLoss + L2 * torch.sum((m.weight)**2)
                if m.bias is not None:
                    biasLoss = biasLoss + L2 * torch.sum((m.bias)**2)
            if isinstance(m, torch.nn.Conv1d):
                weightLoss = weightLoss + L2 * torch.sum((m.weight) ** 2)
        L2Loss = biasLoss + weightLoss
        return(L2Loss)

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
parser = argparse.ArgumentParser(prog='Cell-line vanilla models')
# parser.add_argument('--leaveIn', action='store', default=None)
parser.add_argument('--drug_representation', action='store', default='ECFP4')
parser.add_argument('--max_atoms', action='store', default=115)
parser.add_argument('--num_atom_features', action='store', default=62)
parser.add_argument('--max_degree', action='store', default=5)
parser.add_argument('--num_bond_features', action='store', default=6)
parser.add_argument('--number_of_gcnn_layers', action='store', default=3)
parser.add_argument('--graph_conv_width', action='store', default=128)
parser.add_argument('--number_of_hidden_convolutions', action='store', default=0)
parser.add_argument('--kernel_size', action='store', default=1)
parser.add_argument('--dropout_graph_encoder', action='store', default=0.2,type=float)
parser.add_argument('--lr_gcnn', action='store', default=0.001,type=float)
parser.add_argument('--lr', action='store', default=0.001,type=float)
parser.add_argument('--L2', action='store', default=1e-04,type=float) ### For graph convolutions make it 1e-03
parser.add_argument('--epochs', action='store', default=500) # 400 for gcnns  and 200 in new validation set
parser.add_argument('--batchSize', action='store', default=25)
parser.add_argument('--ANNDense_hidden', action='store', default=2)
parser.add_argument('--model_type', action='store', default=None)
args = parser.parse_args()
# currentId = int(args.leaveIn)
drug_representation = args.drug_representation
model_type = args.model_type
max_atoms = args.max_atoms
num_atom_features = args.num_atom_features
num_bond_features = args.num_bond_features
number_of_gcnn_layers = args.number_of_gcnn_layers
graph_conv_width = args.graph_conv_width
number_of_hidden_convolutions = args.number_of_hidden_convolutions
kernel_size = args.kernel_size
dropout_graph_encoder = args.dropout_graph_encoder
max_degree = args.max_degree
lr = args.lr
lr_gcnn =args.lr_gcnn
L2 = args.L2
maxIter = int(args.epochs)
batchSize = int(args.batchSize)
ANNDense_hidden = int(args.ANNDense_hidden)

assert model_type in ['GCNN','ANN','KNN','SVM'], 'model_type can be only one of the following: [GCNN,ANN,knn,svm]'
assert drug_representation in ['Graph','ECFP4'], 'drug_representation can be only one of the following: [Graph,ECFP4] with default being ECFP4'
if drug_representation=='Graph':
    assert model_type == 'GCNN', 'For graph representation of drugs only graph concolutions can be used (`GCNN`)'
if model_type=='GCNN':
    assert drug_representation == 'Graph', 'For GCNN models only graph representation can be used for drugs (`Graph`)'
        
ConvertToEmpProb = False

TFOutput = pandas.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False, index_col=0)
train_pear = numpy.zeros((len(cell_lines),TFOutput.shape[1]))
val_pear = numpy.zeros((len(cell_lines), TFOutput.shape[1]))
for c,cell in enumerate(cell_lines):
    # Load data
    annotation = pandas.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
    uniprot2gene = dict(zip(annotation['code'], annotation['name']))
    drugInput = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False, index_col=0)
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
    outconds = drugInput.index.values
    outName = TFOutput.columns.values
    outNameGene = [uniprot2gene[x] for x in outName]
    TFOutput = TFOutput.loc[outconds,outName]
    if ConvertToEmpProb==True:
        print2log('Convereted to Emprirical probabilities')
        TFOutput = 1/(1+numpy.exp(-TFOutput))

    ### Validation data
    drugInput_val = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-conditions_drugs.tsv', sep='\t', low_memory=False,
                                    index_col=0)
    drugTargets_val = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3_allcells-drugs_targets.tsv', sep='\t', low_memory=False,
                                      index_col=0)
    TFOutput_val = pandas.read_csv('../preprocessing/preprocessed_data/TF_activities/TrimmedFinal_l1000_allgenes_lvl3_tfs.tsv', sep='\t', low_memory=False,
                                   index_col=0)
    cellInput_val = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/L1000_lvl3-conditions_cells.tsv', sep='\t', low_memory=False, index_col=0)
    TFOutput_val = TFOutput_val.loc[cellInput[cellInput_val[cell] == 1].index, :]
    drugInput_val = drugInput_val.loc[cellInput[cellInput_val[cell] == 1].index, :]
    sigs_to_drop = pandas.read_csv('../preprocessing/preprocessed_data/TrainingValidationData/smiles_and_sigs_to_drop.csv', index_col=0)  # PROSEXE ALLO GIA TA KAINOURIA
    TFOutput_val = TFOutput_val.loc[TFOutput_val.index.isin(sigs_to_drop.sig_id), :]
    drugInput_val = drugInput_val.loc[drugInput_val.index.isin(sigs_to_drop.sig_id), :]
    # Subset input and output to intersecting nodes
    druginName_val = drugInput_val.columns.values
    outconds_val = drugInput_val.index.values
    outName_val = TFOutput_val.columns.values
    outNameGene_val = [uniprot2gene[x] for x in outName_val]
    TFOutput_val = TFOutput_val.loc[outconds_val, outName_val]
    if ConvertToEmpProb == True:
        print2log('Convereted to Emprirical probabilities')
        TFOutput_val = 1 / (1 + numpy.exp(-TFOutput_val))

    Y_val = torch.tensor(TFOutput_val.values, dtype=torch.double).to(device)
    sampleName_val = drugInput_val.index.values

    if drug_representation=='ECFP4':
        ecfp4_fingerprints = []
        for condition in outconds:
            drug = drugInput.columns[numpy.where(drugInput.loc[condition,:].values!=0)][0]
            if ((drug != 'Ctrl_PBS') and (drug != 'Ctrl_Cell')):
                mol = Chem.MolFromSmiles(drug)
                ecfp4_fingerprints.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
            else:
                ecfp4_fingerprints.append(numpy.zeros(1024))
        X = torch.tensor(numpy.array(ecfp4_fingerprints), dtype=torch.double)

        ecfp4_fingerprints = []
        for condition in outconds_val:
            drug = drugInput_val.columns[numpy.where(drugInput_val.loc[condition, :].values != 0)][0]
            if ((drug != 'Ctrl_PBS') and (drug != 'Ctrl_Cell')):
                mol = Chem.MolFromSmiles(drug)
                ecfp4_fingerprints.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
            else:
                ecfp4_fingerprints.append(numpy.zeros(1024))
        X_val = torch.tensor(numpy.array(ecfp4_fingerprints), dtype=torch.double)
    else:
        model_params = {
            "max_atoms" : int(max_atoms), "num_atom_features" : int(num_atom_features), "max_degree" : int(max_degree), "num_bond_features" : int(num_bond_features),
            "number_of_gcnn_layers":int(number_of_gcnn_layers),"number_of_hidden_convolutions":number_of_hidden_convolutions,
            "graph_conv_width" : int(graph_conv_width), "conv1d_in" : int(max_atoms), "kernel_size" : int(kernel_size), "dropout_graph_encoder" : dropout_graph_encoder,
            "output_dim":TFOutput.shape[1]
        }

        train_smiles = []
        for condition in outconds:
            drug = drugInput.columns[numpy.where(drugInput.loc[condition, :].values != 0)][0]
            train_smiles.append(drug)
        val_smiles = []
        for condition in outconds_val:
            drug = drugInput_val.columns[numpy.where(drugInput_val.loc[condition, :].values != 0)][0]
            val_smiles.append(drug)

        X_atoms, X_bonds, X_edges = tensorise_smiles(train_smiles, model_params["max_degree"], model_params["max_atoms"])
        X_atoms = torch.tensor(X_atoms, dtype=torch.double).to(device)
        X_bonds = torch.tensor(X_bonds, dtype=torch.double).to(device)
        X_edges = torch.tensor(X_edges).to(device)

        X_atoms_val, X_bonds_val, X_edges_val = tensorise_smiles(val_smiles, model_params["max_degree"], model_params["max_atoms"])
        X_atoms_val = torch.tensor(X_atoms_val, dtype=torch.double).to(device)
        X_bonds_val = torch.tensor(X_bonds_val, dtype=torch.double).to(device)
        X_edges_val = torch.tensor(X_edges_val).to(device)

    Y = torch.tensor(TFOutput.values, dtype=torch.double)
    #print2log(Y.shape)
    sampleName = drugInput.index.values
    if model_type in ['KNN','SVM']:
        if model_type == 'SVM':
            model = MultiOutputRegressor(SVR())
        else:
            model = KNN()
        model.fit(X.numpy(),Y.numpy())
        Yhat = model.predict(X.numpy())
        Yhat = torch.tensor(Yhat, dtype=torch.double).to(device)
        Y = Y.clone().to(device)
        #print2log(Y)
        train_pear[c, :] = pearson_r(Y.detach(), Yhat.detach()).detach().cpu().numpy()
    else:
        if model_type == 'GCNN':
            Y = Y.float().to(device)
            Y_val = Y_val.float()
            model = GCNN(model_params).to(device)
            #print2log(model)
            criterion = torch.nn.MSELoss(reduction='mean')
            optimizer = torch.optim.Adam(model.parameters(), lr=lr_gcnn, weight_decay=0)  ###Changed to using RMSprop
            resetState = optimizer.state.copy()
            mLoss = criterion(torch.mean(Y, dim=0) * torch.ones(Y.shape).to(device), Y)
            print2log(mLoss)
            N = X_atoms.shape[0]
            e = 0
            trainLossMEAN = []
            trainLossSTD = []
            validationLossMEAN = []
            validationLossSTD = []
            for e in range(e, maxIter):
                curLoss = []
                curLoss_val =[]
                trainloader = bionetwork.getSamples(N, batchSize)
                for dataIndex in trainloader:
                    model.train()
                    optimizer.zero_grad()
                    
                    dataAtomsIn = X_atoms[dataIndex, :].view(len(dataIndex), model_params["max_atoms"],model_params["num_atom_features"])
                    dataEdgesIn = X_edges[dataIndex, :].view(len(dataIndex), model_params["max_atoms"],model_params["max_degree"])
                    dataBondsIn = X_bonds[dataIndex, :].view(len(dataIndex), model_params["max_atoms"],model_params["max_degree"],model_params["num_bond_features"])
                    dataOut = Y[dataIndex, :].view(len(dataIndex), Y.shape[1])

                    Yhat = model(dataAtomsIn,dataBondsIn,dataEdgesIn)

                    fitLoss = criterion(dataOut, Yhat)
                    RegtLoss = model.L2Regularization(L2)

                    loss = fitLoss + RegtLoss
                    loss.backward()
                    optimizer.step()

                    curLoss.append(fitLoss.item())
                    model.eval()
                    Yhat_val = model(X_atoms_val,X_bonds_val,X_edges_val)
                    curLoss_val.append(criterion(Y_val, Yhat_val).item())

                outString = 'Cell-line %s, Epoch=%s / %s' % (cell, e+1, maxIter)
                outString += ', loss={:.4f}'.format(loss.item())
                outString += ', regularization={:.4f}'.format(RegtLoss.item())
                trainLossMEAN.append(numpy.mean(curLoss))
                trainLossSTD.append(numpy.std(curLoss))
                validationLossMEAN.append(numpy.mean(curLoss_val))
                validationLossSTD.append(numpy.std(curLoss_val))
                if e % 50 == 0:
                    print2log(outString)

                if numpy.logical_and(e % 50 == 0, e > 0):
                    optimizer.state = resetState.copy()
            print2log(outString)
            model.eval()
            Yhat = model(X_atoms,X_bonds,X_edges)
            torch.save(model, '../results/vanilla/'+model_type+'/models/' + cell + '.pt')
        else:
            X = X.to(device)
            Y = Y.to(device)
            hiddens = []
            layer = 1024
            for i in range(ANNDense_hidden):
                hiddens.append(int(layer//2))
                layer = int(layer//2)
            model = DenseANN(1024,hiddens,TFOutput.shape[1]).to(device)
            criterion = torch.nn.MSELoss(reduction='mean')
            optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=0)  ###Changed to using RMSprop
            resetState = optimizer.state.copy()
            mLoss = criterion(torch.mean(Y, dim=0) * torch.ones(Y.shape).to(device), Y)
            print2log(mLoss)
            N = X.shape[0]
            e = 0
            trainLossMEAN = []
            trainLossSTD = []
            validationLossMEAN = []
            validationLossSTD = []
            for e in range(e, maxIter):
                curLoss = []
                curLoss_val = []
                trainloader = bionetwork.getSamples(N, batchSize)
                for dataIndex in trainloader:
                    model.train()
                    optimizer.zero_grad()

                    dataIn = X[dataIndex, :].view(len(dataIndex), X.shape[1])
                    dataOut = Y[dataIndex, :].view(len(dataIndex), Y.shape[1])

                    Yhat = model(dataIn)

                    fitLoss = criterion(dataOut, Yhat)
                    RegtLoss = model.L2Regularization(L2)

                    loss = fitLoss + RegtLoss
                    loss.backward()
                    optimizer.step()

                    curLoss.append(fitLoss.item())
                    model.eval()
                    Yhat_val = model(X_val.to(device))
                    curLoss_val.append(criterion(Y_val,Yhat_val).item())

                outString = 'Cell-line %s, Epoch=%s / %s' % (cell, e+1, maxIter)
                outString += ', loss={:.4f}'.format(loss.item())
                outString += ', regularization={:.4f}'.format(RegtLoss.item())
                trainLossMEAN.append(numpy.mean(curLoss))
                trainLossSTD.append(numpy.std(curLoss))
                validationLossMEAN.append(numpy.mean(curLoss_val))
                validationLossSTD.append(numpy.std(curLoss_val))
                if e % 50 == 0:
                    print2log(outString)

                if numpy.logical_and(e % 50 == 0, e > 0):
                    optimizer.state = resetState.copy()
            print2log(outString)
            model.eval()
            Yhat = model(X)
            torch.save(model, '../results/vanilla/' + model_type + '/models/' + cell +'.pt')

        #print2log(Y.shape)
        #print2log(Yhat.shape)
        #print2log(X.shape)
        # %%
        plt.rcParams["figure.figsize"] = (9, 6)
        plt.figure()
        T = numpy.array(range(len(trainLossMEAN)))
        plt.plot(T, trainLossMEAN,label='train')
        plt.xlim([0, len(T)])
        curColor = plt.gca().lines[-1].get_color()
        plt.fill_between(T,
                        numpy.array(trainLossMEAN) - numpy.array(trainLossSTD),
                        numpy.array(trainLossMEAN) + numpy.array(trainLossSTD),
                        color=curColor, alpha=0.2)
        plt.plot(T, validationLossMEAN,label='validation')
        plt.xlim([0, len(T)])
        curColor = plt.gca().lines[-1].get_color()
        plt.fill_between(T,
                         numpy.array(validationLossMEAN) - numpy.array(validationLossSTD),
                         numpy.array(validationLossMEAN) + numpy.array(validationLossSTD),
                         color=curColor, alpha=0.2)
        plt.plot([0, len(T)], numpy.array([1, 1]) * mLoss.item(), 'black', linestyle='--',label='train mean MSE')
        plt.legend(loc="upper right")
        plt.ylabel('Loss')
        plt.xlabel('Epoch')
        plt.yscale('log')
        plt.savefig('../results/vanilla/' + model_type + '/train/fig1_' + cell  +'.png', dpi=600)

    train_pear[c, :] = pearson_r(Y.detach(), Yhat.detach()).detach().cpu().numpy()
    plt.rcParams["figure.figsize"] = (9, 6)
    plt.figure()
    plt.scatter(Yhat.detach().cpu(), Y.detach().cpu(), alpha=0.2)
    plotting.lineOfIdentity()
    plotting.addCorrelation(Yhat.detach().cpu(), Y.detach().cpu())
    plt.xlabel('Fit')
    plt.ylabel('Data')
    plt.gca().axis('equal')
    plt.gca().set_xticks([0, 0.5, 1])
    plt.gca().set_yticks([0, 0.5, 1])
    plt.savefig('../results/vanilla/'+model_type+'/train/fig4_' + cell  + '.png', dpi=600)
    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu(), Y.detach().cpu(), outNameGene)
    plt.savefig('../results/vanilla/'+model_type+'/train/fig7_' + cell+ '.png', dpi=600)
    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu().T, Y.detach().cpu().T, sampleName)
    plt.savefig('../results/vanilla/'+model_type+'/train/fig8_' + cell+ '.png', dpi=600)

    # Performance in validation
    if model_type in ['KNN','SVM']:
        Yhat = model.predict(X_val)
        Yhat = torch.tensor(Yhat, dtype=torch.double).to(device)
    else:
        if model_type == 'ANN':
            X_val = X_val.to(device)
            model.eval()
            Yhat = model(X_val)
        else:
            model.eval()
            Yhat = model(X_atoms_val,X_bonds_val,X_edges_val)



    val_pear[c, :] = pearson_r(Y_val.detach(), Yhat.detach()).detach().cpu().numpy()
    if len(numpy.where(numpy.isnan(val_pear[c, :]))[0]):
        inds = numpy.where(numpy.isnan(val_pear[c, :]))[0]
        r = pearson_r(Y_val[:, inds].detach(),Yhat[:, inds].detach() + 1e-8 * torch.randn(Yhat.shape[0], Yhat[:, inds].shape[1]).to(device))
        val_pear[c, inds] = r.detach().cpu().numpy()
    plt.figure()
    plt.scatter(Yhat.detach().cpu().numpy(), Y_val.detach().cpu().numpy(), alpha=0.2)
    plotting.lineOfIdentity()
    plotting.addCorrelation(Yhat.detach().cpu(), Y_val.detach().cpu())
    plt.xlabel('Fit')
    plt.ylabel('Data')
    plt.gca().axis('equal')
    plt.gca().set_xticks([0, 0.5, 1])
    plt.gca().set_yticks([0, 0.5, 1])
    plt.savefig('../results/vanilla/'+model_type+'/test/performance_'+cell+'.png',dpi=600)

    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu(), Y_val.detach().cpu(), outNameGene_val)
    plt.savefig('../results/vanilla/'+model_type+'/test/performancePerTF_'+cell+'.png',dpi=600)
    plt.figure()
    rank = plotting.compareAllTFs(Yhat.detach().cpu().T, Y_val.detach().cpu().T, outconds_val)
    plt.savefig('../results/vanilla/'+model_type+'/test/performancePerSample_'+cell+'.png',dpi=600)

    print2log('Finished cell-line %s' % cell)

train_pear = pandas.DataFrame(train_pear)
train_pear.columns = TFOutput.columns
train_pear.index = cell_lines
train_pear.to_csv('../results/vanilla/'+model_type+'/train_Performance_perTF.csv')
print2log(train_pear.mean(1).mean())

val_pear = pandas.DataFrame(val_pear)
val_pear.columns = TFOutput.columns
val_pear.index = cell_lines
val_pear.to_csv('../results/vanilla/'+model_type+'/val_easy_Performance_perTF.csv')
print2log(val_pear.mean(1).mean())