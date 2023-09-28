import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
# from matplotlib import pyplot as plt
# import seaborn as sns
import argparse
import logging
# sns.set()

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

def getZeroBiasIndices(node_dictionary,Nodes2Keep):
    nodes_inds = []
    for node_name,node_ind in node_dictionary.items():
        if node_ind not in Nodes2Keep:
            nodes_inds.append(node_ind)
    return(np.array(nodes_inds))

def getZeroWeightsIndices(A,Edges2Keep):
    edges_inds = []
    nonzero_indices = np.array(A.nonzero())
    for index, (source, target) in enumerate(zip(*nonzero_indices)):
        if (source,target) not in Edges2Keep:
            edges_inds.append(index)
    return(np.array(edges_inds))

parser = argparse.ArgumentParser(prog='Infer MoA')
parser.add_argument('--ensembles_path', action='store', default="../results/A375_ensembles/")
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
# interestingSamples = pd.read_csv(ensembles_path+'interestingSamples.csv',index_col=1)
TF = "Q08050"
TF_gene = "FOXM1"
drug = "C[C@]12O[C@H](C[C@]1(O)CO)n1c3ccccc3c3c4C(=O)NCc4c4c5ccccc5n2c4c13"
drug_name = "lestaurtinib"
sample = 'CPC014_A375_6H:BRD-K23192422-001-01-1:10'
nodes = ['P06493','P24941']
nodes_names = ['CDK1','CDK2']
targets = ['P08588','P00533','P04629','P36888']
targets_name = ['ADRB1','EGFR','NTRK1','FLT3']


### Load interesting in-silico validations
experiments = pd.read_csv('../results/ExperimentalValidation/inSilicoValResults.csv',index_col=0)

#%%

### Load network
#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for A375 was 120
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

TF_ind = np.where(TFOutput.columns.values == TF)[0]
node_ind = [np.where(drugTargets.columns == node)[0][0] for node in nodes]
target_ind = [np.where(drugTargets.columns == target)[0][0] for target in targets]
#%% Get results for of , on, on+off target effect
dmso_ind = np.where(drugInput.columns.values=='CS(C)=O')[0][0]
X = torch.tensor(drugInput.loc[sample,:].values.copy(), dtype=torch.double)
X = X.unsqueeze(dim=0)
all_effects = []
on_target = []
off_target = []
cdk1_KO = []
cdk2_KO = []
cdk6_KO = []
adrb1_KO = []
ntrk1_KO = []
flt3_KO = []
egfr_KO = []
dmso_simulation = []
egfr_cdk2 = []
flt3_cdk2 = []
# print2log(targets_name)
# print2log(nodes_names)
for i in range(numberOfModels):
    model = torch.load(inputPath+str(i)+".pt")
    model.eval()
    Yhat, YhatFull = model(X)
    mask = torch.mm(1.0*(X!=0).double(),model.drugLayer.mask.T)
    
    Xin = model.drugLayer(X)
    
    # on-target only
    Xin_masked = torch.mul(Xin, mask)
    fullX_masked = model.inputLayer(Xin_masked)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_1 = model.projectionLayer(YhatFull_masked)
    
    # off-target only
    Xin_masked_2 = torch.mul(Xin, 1.0 - mask)
    fullX_masked = model.inputLayer(Xin_masked_2)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_2 = model.projectionLayer(YhatFull_masked)
    
    all_effects.append(Yhat[:,TF_ind].item())
    on_target.append(Yhat_masked_1[:,TF_ind].item())
    off_target.append(Yhat_masked_2[:,TF_ind].item())
    # node KO only
    # CDK1 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,node_ind[0]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    cdk1_KO.append(Yhat_ko[:,TF_ind].item())
    # CDK2 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,node_ind[1]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    cdk2_KO.append(Yhat_ko[:,TF_ind].item())
    
    # target KO only
    # ADRB1 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[0]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    adrb1_KO.append(Yhat_ko[:,TF_ind].item())
    # EGFR KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[1]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    egfr_KO.append(Yhat_ko[:,TF_ind].item())
    # NTRK1 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[2]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    ntrk1_KO.append(Yhat_ko[:,TF_ind].item())
    # FLT3 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[3]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    flt3_KO.append(Yhat_ko[:,TF_ind].item())
    
    # combinations KOs
    # CDK1+FLT3 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[3]] = -10
    Xin_ko[:,node_ind[0]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    flt3_cdk2.append(Yhat_ko[:,TF_ind].item())
    # CDK1+EGFR KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[1]] = -10
    Xin_ko[:,node_ind[0]] = -10
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    egfr_cdk2.append(Yhat_ko[:,TF_ind].item())
    
    #dmso simulation
    X_dmso = torch.zeros(X.shape)
    X_dmso[:,dmso_ind]=1.04139269
    X_dmso = X_dmso.double()
    Yhat_dmso, YhatFull_dmso = model(X_dmso)
    dmso_simulation.append(Yhat_dmso[:,TF_ind].item())
    
    print2log('Finished model %s'%i)

insilico_results = pd.DataFrame({'DMSO':dmso_simulation,'lestaurtinib':all_effects,
                                 'on-target':on_target,'off-target':off_target,
                                 'ADRB1':adrb1_KO,'EGFR':egfr_KO,
                                 'NTRK1':ntrk1_KO,'FLT3':flt3_KO,
                                 'CDK1':cdk1_KO,'CDK2':cdk2_KO,
                                 'CDK1+EGFR':egfr_cdk2,'CDK1+FLT3':flt3_cdk2})
insilico_results.to_csv('../results/ExperimentalValidation/inSilicoKOs_minus10.csv')
#%%
# Perform the same analysis but using the reduced network
NetworkPandas = pd.read_csv(ensembles_path+'MoA/'+drug_name+'_any/'+cell+'_'+drug_name+'_'+TF_gene+'_moa_model_ensembleFiltered.csv',index_col=0)
trimmed_nodes = np.union1d(NetworkPandas.source.values,NetworkPandas.target.values)
trimmed_edges = []
for i in range(len(NetworkPandas)):
    trimmed_edges.append((NetworkPandas.source[i],NetworkPandas.target[i]))
dmso_ind = np.where(drugInput.columns.values=='CS(C)=O')[0][0]
X = torch.tensor(drugInput.loc[sample,:].values.copy(), dtype=torch.double)
X = X.unsqueeze(dim=0)
all_effects = []
on_target = []
off_target = []
cdk1_KO = []
cdk2_KO = []
cdk6_KO = []
adrb1_KO = []
ntrk1_KO = []
flt3_KO = []
egfr_KO = []
dmso_simulation = []
egfr_cdk2 = []
flt3_cdk2 = []
print2log(targets_name)
print2log(nodes_names)
for i in range(numberOfModels):
    model = torch.load(inputPath+str(i)+".pt")
    model.network.weights.requires_grad = False
    model.network.bias.requires_grad = False
    model.eval()
    dictionary = dict(zip(model.network.nodeNames, list(range(len(model.network.nodeNames)))))
    bias_inds = getZeroBiasIndices(dictionary,trimmed_nodes)
    weight_inds = getZeroWeightsIndices(model.network.A,trimmed_edges)
    model.network.bias[bias_inds,:] = 0.
    model.network.weights[weight_inds] = 0.
    model.network.A.data = model.network.weights
    
    Yhat, YhatFull = model(X)
    mask = torch.mm(1.0*(X!=0).double(),model.drugLayer.mask.T)
    
    Xin = model.drugLayer(X)
    
    # on-target only
    Xin_masked = torch.mul(Xin, mask)
    fullX_masked = model.inputLayer(Xin_masked)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_1 = model.projectionLayer(YhatFull_masked)
    
    # off-target only
    Xin_masked_2 = torch.mul(Xin, 1.0 - mask)
    fullX_masked = model.inputLayer(Xin_masked_2)
    YhatFull_masked = model.network(fullX_masked)
    Yhat_masked_2 = model.projectionLayer(YhatFull_masked)
    
    all_effects.append(Yhat[:,TF_ind].item())
    on_target.append(Yhat_masked_1[:,TF_ind].item())
    off_target.append(Yhat_masked_2[:,TF_ind].item())
    # node KO only
    # CDK1 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,node_ind[0]] = -100
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    cdk1_KO.append(Yhat_ko[:,TF_ind].item())
    # CDK2 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,node_ind[1]] = -100
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    cdk2_KO.append(Yhat_ko[:,TF_ind].item())
    
    # # target KO only
    # # ADRB1 KOs
    # Xin_ko = torch.zeros(1,Xin.shape[1])
    # Xin_ko[:,target_ind[0]] = -100
    # fullX_ko = model.inputLayer(Xin_ko)
    # YhatFull_masked = model.network(fullX_ko)
    # Yhat_ko = model.projectionLayer(YhatFull_masked)
    # adrb1_KO.append(Yhat_ko[:,TF_ind].item())
    # # EGFR KOs
    # Xin_ko = torch.zeros(1,Xin.shape[1])
    # Xin_ko[:,target_ind[1]] = -100
    # fullX_ko = model.inputLayer(Xin_ko)
    # YhatFull_masked = model.network(fullX_ko)
    # Yhat_ko = model.projectionLayer(YhatFull_masked)
    # egfr_KO.append(Yhat_ko[:,TF_ind].item())
    # # NTRK1 KOs
    # Xin_ko = torch.zeros(1,Xin.shape[1])
    # Xin_ko[:,target_ind[2]] = -100
    # fullX_ko = model.inputLayer(Xin_ko)
    # YhatFull_masked = model.network(fullX_ko)
    # Yhat_ko = model.projectionLayer(YhatFull_masked)
    # ntrk1_KO.append(Yhat_ko[:,TF_ind].item())
    # # FLT3 KOs
    # Xin_ko = torch.zeros(1,Xin.shape[1])
    # Xin_ko[:,target_ind[3]] = -100
    # fullX_ko = model.inputLayer(Xin_ko)
    # YhatFull_masked = model.network(fullX_ko)
    # Yhat_ko = model.projectionLayer(YhatFull_masked)
    # flt3_KO.append(Yhat_ko[:,TF_ind].item())
    
    # combinations KOs
    # CDK1+FLT3 KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[3]] = -100
    Xin_ko[:,node_ind[0]] = -100
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    flt3_cdk2.append(Yhat_ko[:,TF_ind].item())
    # CDK1+EGFR KOs
    Xin_ko = torch.zeros(1,Xin.shape[1])
    Xin_ko[:,target_ind[1]] = -100
    Xin_ko[:,node_ind[0]] = -100
    fullX_ko = model.inputLayer(Xin_ko)
    YhatFull_masked = model.network(fullX_ko)
    Yhat_ko = model.projectionLayer(YhatFull_masked)
    egfr_cdk2.append(Yhat_ko[:,TF_ind].item())
    
    #dmso simulation
    X_dmso = torch.zeros(X.shape)
    X_dmso[:,dmso_ind]=1.04139269
    X_dmso = X_dmso.double()
    Yhat_dmso, YhatFull_dmso = model(X_dmso)
    dmso_simulation.append(Yhat_dmso[:,TF_ind].item())
    
    print2log('Finished model %s'%i)
    
insilico_results = pd.DataFrame({'DMSO':dmso_simulation,'lestaurtinib':all_effects,
                                 'on-target':on_target,'off-target':off_target,
                                 'ADRB1':adrb1_KO,'EGFR':egfr_KO,
                                 'NTRK1':ntrk1_KO,'FLT3':flt3_KO,
                                 'CDK1':cdk1_KO,'CDK2':cdk2_KO,
                                 'CDK1+EGFR':egfr_cdk2,'CDK1+FLT3':flt3_cdk2})
insilico_results.to_csv('../results/ExperimentalValidation/reduced_inSilicoKOs_minus100.csv')

#%% Get KO level sensitivity analysis
levels = [0,1,3,5,15,10,25,40,50,60,75,100,150,200]
insilico_results_all = pd.DataFrame()
insilico_node_ko_all = pd.DataFrame()
nodesNameGene = [uniprot2gene[x] for x in nodeNames]
for lvl in levels:
    print2log('Start KO level = %s'%lvl)
    dmso_ind = np.where(drugInput.columns.values=='CS(C)=O')[0][0]
    X = torch.tensor(drugInput.loc[sample,:].values.copy(), dtype=torch.double)
    X = X.unsqueeze(dim=0)
    all_effects = []
    on_target = []
    off_target = []
    cdk1_KO = []
    cdk2_KO = []
    cdk6_KO = []
    adrb1_KO = []
    ntrk1_KO = []
    flt3_KO = []
    egfr_KO = []
    dmso_simulation = []
    egfr_cdk2 = []
    flt3_cdk2 = []
    cdk1_node = []
    cdk2_node = []
    adrb1_node = []
    ntrk1_node = []
    flt3_node = []
    egfr_node = []
    # print2log(targets_name)
    # print2log(nodes_names)
    for i in range(numberOfModels):
        model = torch.load(inputPath+str(i)+".pt")
        model.eval()
        Yhat, YhatFull = model(X)
        mask = torch.mm(1.0*(X!=0).double(),model.drugLayer.mask.T)
        
        Xin = model.drugLayer(X)
        
        # node KO only
        # CDK1 KOs
        node_ind_ko = np.where(np.array(nodesNameGene)=='CDK1')[0][0]
        Xin_ko = torch.zeros(1,Xin.shape[1])
        Xin_ko[:,node_ind[0]] = -lvl
        fullX_ko = model.inputLayer(Xin_ko)
        YhatFull_masked = model.network(fullX_ko)
        Yhat_ko = model.projectionLayer(YhatFull_masked)
        cdk1_KO.append(Yhat_ko[:,TF_ind].item())
        cdk1_node.append(YhatFull_masked[:,node_ind_ko].item())
        # CDK2 KOs
        node_ind_ko = np.where(np.array(nodesNameGene)=='CDK2')[0][0]
        Xin_ko = torch.zeros(1,Xin.shape[1])
        Xin_ko[:,node_ind[1]] = -lvl
        fullX_ko = model.inputLayer(Xin_ko)
        YhatFull_masked = model.network(fullX_ko)
        Yhat_ko = model.projectionLayer(YhatFull_masked)
        cdk2_KO.append(Yhat_ko[:,TF_ind].item())
        cdk2_node.append(YhatFull_masked[:,node_ind_ko].item())
        
        # target KO only
        # ADRB1 KOs
        node_ind_ko = np.where(np.array(nodesNameGene)=='ADRB1')[0][0]
        Xin_ko = torch.zeros(1,Xin.shape[1])
        Xin_ko[:,target_ind[0]] = -lvl
        fullX_ko = model.inputLayer(Xin_ko)
        YhatFull_masked = model.network(fullX_ko)
        Yhat_ko = model.projectionLayer(YhatFull_masked)
        adrb1_KO.append(Yhat_ko[:,TF_ind].item())
        adrb1_node.append(YhatFull_masked[:,node_ind_ko].item())
        # EGFR KOs
        node_ind_ko = np.where(np.array(nodesNameGene)=='EGFR')[0][0]
        Xin_ko = torch.zeros(1,Xin.shape[1])
        Xin_ko[:,target_ind[1]] = -lvl
        fullX_ko = model.inputLayer(Xin_ko)
        YhatFull_masked = model.network(fullX_ko)
        Yhat_ko = model.projectionLayer(YhatFull_masked)
        egfr_KO.append(Yhat_ko[:,TF_ind].item())
        egfr_node.append(YhatFull_masked[:,node_ind_ko].item())
        # NTRK1 KOs
        node_ind_ko = np.where(np.array(nodesNameGene)=='NTRK1')[0][0]
        Xin_ko = torch.zeros(1,Xin.shape[1])
        Xin_ko[:,target_ind[2]] = -lvl
        fullX_ko = model.inputLayer(Xin_ko)
        YhatFull_masked = model.network(fullX_ko)
        Yhat_ko = model.projectionLayer(YhatFull_masked)
        ntrk1_KO.append(Yhat_ko[:,TF_ind].item())
        ntrk1_node.append(YhatFull_masked[:,node_ind_ko].item())
        # FLT3 KOs
        node_ind_ko = np.where(np.array(nodesNameGene)=='FLT3')[0][0]
        Xin_ko = torch.zeros(1,Xin.shape[1])
        Xin_ko[:,target_ind[3]] = -lvl
        fullX_ko = model.inputLayer(Xin_ko)
        YhatFull_masked = model.network(fullX_ko)
        Yhat_ko = model.projectionLayer(YhatFull_masked)
        flt3_KO.append(Yhat_ko[:,TF_ind].item())
        flt3_node.append(YhatFull_masked[:,node_ind_ko].item())
        
        
        # if i % 5 == 0:
        #     print2log('Finished model %s'%i)

    insilico_results = pd.DataFrame({'ADRB1':adrb1_KO,'EGFR':egfr_KO,
                                     'NTRK1':ntrk1_KO,'FLT3':flt3_KO,
                                     'CDK1':cdk1_KO,'CDK2':cdk2_KO,
                                     'level':np.repeat(lvl,len(cdk1_KO))})
    insilico_node_ko = pd.DataFrame({'ADRB1':adrb1_node,'EGFR':egfr_node,
                                     'NTRK1':ntrk1_node,'FLT3':flt3_node,
                                     'CDK1':cdk1_node,'CDK2':cdk2_node,
                                     'level':np.repeat(lvl,len(cdk1_KO))})
    insilico_results_all = insilico_results_all.append(insilico_results)
    insilico_node_ko_all = insilico_node_ko_all.append(insilico_node_ko)
    
# Save results
insilico_results_all.to_csv('../results/ExperimentalValidation/inSilicoKOs_dose_response.csv')
insilico_node_ko_all.to_csv('../results/ExperimentalValidation/inSilicoKOs_node_response.csv')

