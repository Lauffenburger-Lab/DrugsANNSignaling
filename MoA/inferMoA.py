import torch
import pandas as pd
import numpy as np
import bionetworkWithDrugs as bionetwork
import scipy
from matplotlib import pyplot as plt
import networkx as nx
import seaborn as sns
import argparse
import logging
from pathlib import Path
from operator import itemgetter
sns.set()

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
print2log = logger.info

def remakeNetworkList(networkList,model_weights,weight_mask,nodes):
    import numpy as np
    sources = []
    targets = []
    weights = []
    for n in range(networkList.shape[1]):
        if networkList[0,n] in nodes:
            if networkList[1,n] in nodes:
                if weight_mask[n]>0:
                    sources.append(networkList[0,n])
                    targets.append(networkList[1, n])
                    weights.append(model_weights[n])
    networkList_new = np.array([sources,targets])
    weights = np.array(weights)
    return networkList_new,weights

def resetGradients(model):
    for var in model.parameters():
        var.grad = None


def inferDrugTarget(interactions, global_thresholds, prior_mask, drugInput, drugTargets, model_no, save_file_pattern,
                    grad_thresholds=list(np.logspace(-3.5, 3.5, num=45)), thresh=0.25):
    global_thresholds = global_thresholds.detach().numpy()
    scores = torch.abs(torch.tensor(interactions.values))
    grad_thresh = global_thresholds[model_no]
    # Get new mask with drug-targets interactions
    mask = 1.0 * (scores >= grad_thresh)
    mask = mask.detach().numpy()

    # Create merged dataframe
    df_mask = pd.DataFrame(prior_mask.numpy())
    df_mask.index = drugInput.columns
    df_mask.columns = drugTargets.columns
    df_mask.reset_index(inplace=True)
    df_mask = df_mask.rename(columns={'index': 'drug'})
    df_mask = pd.melt(df_mask, id_vars=['drug'])
    df_mask['value'] = df_mask['value'].transform(lambda x: 'Interaction' if abs(x) > 0 else 'No interaction')
    df_mask = df_mask.rename(columns={'value': 'Prior knowledge'})
    # New interactios
    df_masked_scores = pd.DataFrame(mask)
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
    merged_df.to_csv(save_file_pattern + '_mergedInteractions_%s.csv' % model_no)
    return (merged_df)

def pathFinder(source,target,network,mode='Shortest',max_depth = 5):
    paths = []
    hasPath = nx.has_path(network,source,target)
    if hasPath==True:
        if mode=='Shortest':
            paths = [p for p in nx.all_shortest_paths(network, source, target)]
        else:
            paths = [p for p in nx.all_simple_paths(network, source, target,max_depth)]
    return paths

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
sample = "CPC014_A375_6H:BRD-K23192422-001-01-1:10"
moa_off_target = 'any'
#'inhibit'

# deltaTF = interestingSamples.loc[sample,"delta"]

# Make drug directory
Path(ensembles_path+'MoA/'+drug_name).mkdir(parents=True, exist_ok=True)
Path(ensembles_path + 'InteractionScores/MergedInteractions/'+drug_name).mkdir(parents=True, exist_ok=True)

#%%
### Get drug-target interactions
# merged_interaction = pd.read_csv(ensembles_path+'InteractionScores/l1000_modeltype4_lamda6_'+cell+'_mergedInteractions_ROC.csv',index_col=0)
# merged_interaction = merged_interaction[merged_interaction['drug']==drug]
# merged_interaction = merged_interaction[merged_interaction['Inferred']=='Interaction']

### Load network
#Load network
networkList, nodeNames, modeOfAction = bionetwork.loadNetwork('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv')
annotation = pd.read_csv('../preprocessing/preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t')
uniprot2gene = dict(zip(annotation['code'], annotation['name']))
bionetParams = bionetwork.trainingParameters(iterations = 120, clipping=1, targetPrecision=1e-6, leak=0.01) # for A549 was 120
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

# Χ_all = torch.tensor(drugInput.values.copy(), dtype=torch.double)
# Y_all = torch.tensor(TFOutput.values, dtype=torch.double)
# Keep only drug/sample of interest
TF_ind = np.where(TFOutput.columns==TF)[0]
sample_ind = np.where(TFOutput.index==sample)[0]
drugInput = pd.DataFrame(drugInput.loc[sample,:]).T
TFOutput = pd.DataFrame(TFOutput.loc[sample,:]).T

X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
Y = torch.tensor(TFOutput.values, dtype=torch.double)


### Create matrix to range doses from 0 to 100% of the input signal after the drugLayer
ratioMatrix = torch.arange(0, 1.01, 0.01)
ratioMatrix = ratioMatrix.T.repeat(drugTargets.shape[1],1)
ratioMatrix = ratioMatrix.double()

drug_ratioMatrix = torch.arange(0, 1.01, 0.01)
drug_ratioMatrix = drug_ratioMatrix.T.repeat(drugInput.shape[1],1)
drug_ratioMatrix = drug_ratioMatrix.double()

### Tresholds for removing edges based on the weight gradients
thresholds = list(np.linspace(start=0, stop=60, num=3000))
thresholds = thresholds[1:]
# if moa_off_target == 'inhibit':
#     thresholds_w = list(np.linspace(start=-20, stop=0, num=2000))
# else:
#     thresholds_w = list(np.linspace(start=0, stop=20, num=2000))
thresholds_w = list(np.linspace(start=0, stop=20, num=2000))
thresholds_w = thresholds_w[1:]

### Get mapping of network nodes to TFs
dictionary = dict(zip(nodeNames, list(range(len(nodeNames)))))
nodeOrder = np.array([dictionary[x] for x in TFOutput.columns])
TF_index_inNet = nodeOrder[TF_ind]

# Get global grad score
# global_grad_scores = torch.load(ensembles_path+"all_drugs_global_thresholds.pt")
global_grad_scores = pd.read_csv(ensembles_path+"all_drugs_global_thresholds.csv",index_col=0)
global_grad_scores = torch.tensor(global_grad_scores.loc[drug,:].values)
# global_grad_scores = torch.load(ensembles_path+"InteractionScores/global_gradient_scores_"+drug+"_sample_ind"+str(sample_ind[0])+"_all_models.pt")

### Mapping of names and ids in graph
map = dict(zip(nodeNames, list(range(len(nodeNames)))))
df_all = pd.DataFrame({'source':[],'target':[],'weight':[],'name_source':[],'name_target':[],'interaction':[],'model_no':[]})
all_sources = pd.DataFrame({'sources':[],'model':[]})
df_all_source_types = pd.DataFrame({'type':[],'name':[]})
df_all_target_types = pd.DataFrame({'type':[],'name':[]})
mean_bias = np.zeros((len(nodeNames)))
for i in range(numberOfModels):#range(1):
    model = torch.load(inputPath+str(i)+".pt")
    resetGradients(model)
    model.eval()
    mean_bias = mean_bias + model.network.bias.detach().numpy().squeeze()
    X = torch.tensor(drugInput.values.copy(), dtype=torch.double)
    Y = torch.tensor(TFOutput.values, dtype=torch.double)
    Xin =  model.drugLayer(X)
    Xin =  ratioMatrix.T * Xin.squeeze()
    Xin = Xin.detach()

    # merged_interaction = pd.read_csv(ensembles_path+'InteractionScores/l1000_modeltype4_lamda6_'+cell+'_mergedInteractions_ROC.csv',index_col=0)
    # merged_interaction = merged_interaction[merged_interaction['drug']==drug]
    # merged_interaction = merged_interaction[merged_interaction['Inferred']=='Interaction']

    Yin = model.inputLayer(Xin)
    YhatFull = model.network(Yin)
    Yhat = model.projectionLayer(YhatFull)
    L_objective = torch.sum(Yhat[:,TF_ind])
    L_objective.backward()
    
    # # Get range of input activity for different concetrations of drugs
    # X_binned =  drug_ratioMatrix.T * X.squeeze()
    # Xin_binned =  model.drugLayer(X_binned)
    # Xin_range = torch.abs(torch.max(Xin_binned,0)[0] - torch.min(Xin_binned,0)[0])
    

    ### Get drug-target interactions
    interactions = pd.read_csv(ensembles_path + 'InteractionScores/l1000_modeltype4_lamda6_' + cell + '_interactionScores_%s.csv' % i,index_col=0)
    # drugs = interactions.index.values
    # target_space = interactions.columns.values
    # interactions = torch.tensor(interactions.values) * Xin_range
    # interactions = pd.DataFrame(interactions.detach().numpy())
    # interactions.index =drugs
    # interactions = interactions.T
    # interactions.index = target_space
    # interactions = interactions.T
    merged_interactions = inferDrugTarget(interactions,global_grad_scores,
                                          model.drugLayer.mask.T.detach(), drugInput,drugTargets,
                                          i,
                                          ensembles_path + 'InteractionScores/MergedInteractions/'+drug_name+'/'+'l1000_modeltype4_lamda6_' + cell+'_'+drug_name+"_sample_ind"+str(sample_ind[0]))
    drug_target_nodes = merged_interactions[merged_interactions['drug'] == drug]
    drug_target_nodes = drug_target_nodes[drug_target_nodes['Inferred'] == 'Interaction']
    drug_target_nodes = drug_target_nodes.variable.unique()
    # print2log('Unique drug-targets after thresholding %s'%len(drug_target_nodes))
    
    if len(drug_target_nodes)==0:
        continue

    ### calculate gradients
    db = model.network.bias.grad.T.squeeze()
    dw = model.network.weights.grad
    if moa_off_target == 'inhibit':
        ind = np.where(db[np.array([map[map_id] for map_id in drug_target_nodes])]<0)[0]
    elif moa_off_target == 'any':
        ind = np.arange(len(drug_target_nodes))
    else:
        ind = np.where(db[np.array([map[map_id] for map_id in drug_target_nodes])]>0)[0]
    drug_target_nodes = drug_target_nodes[ind]
    print2log('Unique drug-targets after thresholding with desired MoA on the TF %s'%len(drug_target_nodes))
    if len(drug_target_nodes)==0:
        continue
    ### calculate variance of YhatFull.
    ### I want to keep nodes with both high variance due to input signal change and high gradient
    Range = torch.abs(torch.max(YhatFull,0)[0] - torch.min(YhatFull,0)[0])
    #Range = torch.abs(YhatFull[-1,:] - YhatFull[0,:])
    nodeScore = torch.abs(db) * Range
    #var = YhatFull.var(0)
    #nodeScore = torch.abs(db) * var
    weightScore = torch.abs(dw) * torch.abs(model.network.weights)
    
    sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    # drug_target_nodes = drug_target_nodes[torch.abs(db[sources]).detach().numpy() >=0.01]
    # sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    target = itemgetter(TF)(map)
    
    total_edges = []
    keep_source = []
    keep_source_w = []
    
    for source_id,source in enumerate(sources):
        ### Remove nodes
        A = scipy.sparse.csr_matrix((model.network.weights.detach(), networkList))#,shape=(YhatFull.shape[1], YhatFull.shape[1]), dtype='float64')
        SignalingNet = nx.from_numpy_array(A.toarray(),parallel_edges=True,create_using= nx.DiGraph())
        j = 0
        flag = True
        while (flag==True and j<len(thresholds)):
            th = thresholds[j]
            RemoveNode = list(torch.where(nodeScore < th)[0].numpy())
            if source in RemoveNode:
                RemoveNode.remove(source)
            if target in RemoveNode:
                RemoveNode.remove(target)
            SignalingNet.remove_nodes_from(RemoveNode)
            has_path = nx.has_path(SignalingNet, source, target)
            if has_path==False:
                flag = False
            else:
                j += 1
        if (j-1)<0:
            keep_source.append(False)
            #continue
        else:
            keep_source.append(True)
            nodesToKeep = list(torch.where(nodeScore >= thresholds[j-1])[0].numpy())
            RemoveNode = list(torch.where(nodeScore < thresholds[j-1])[0].numpy())
            if source in RemoveNode:
                RemoveNode.remove(source)
                nodesToKeep.append(source)
            if target in RemoveNode:
                RemoveNode.remove(target)
                nodesToKeep.append(source)
            SignalingNet = nx.from_numpy_array(A.toarray(),parallel_edges=True,create_using= nx.DiGraph())
            SignalingNet.remove_nodes_from(RemoveNode)
            G_tmp = SignalingNet.copy()
            
            graphs = list(nx.connected_components(G_tmp.to_undirected()))
            graphsWithTF = []
            for graph in graphs:
                if TF_index_inNet in list(graph):
                    graphsWithTF.append(graph)
            largest_cc = max(graphsWithTF, key=len)
            H = G_tmp.subgraph(largest_cc)
            H = H.to_directed()
            nodesToKeepFinal = np.array(H.nodes)
            SignalingNet = nx.from_numpy_array(A.toarray(), parallel_edges=True, create_using=nx.DiGraph())
            SignalingNet = SignalingNet.edge_subgraph(list(H.edges)).copy()
            SignalingNet = nx.freeze(SignalingNet).copy()
            SignalingNet = nx.DiGraph(SignalingNet).copy()
            nodes_all = list(SignalingNet.nodes())
            remove_final = []
            for node in nodes_all:
                if node != target:
                    if node!=source:
                        if not nx.has_path(SignalingNet, node, target):
                            remove_final.append(node)
            SignalingNet.remove_nodes_from(remove_final)
            nodes_all = list(SignalingNet.nodes())
            # print('Finished removing nodes from %s to %s'%(drug_target_nodes[source_id],TF))
        ### Remove edges
        j = 0
        flag = True
        H = SignalingNet.copy()
        while (flag==True and j<len(thresholds_w)):
            G_tmp = SignalingNet.copy()
            th = thresholds_w[j]
            dw_mask = weightScore < th
            dw_mask = dw_mask.numpy()
            newNetworkList = networkList[:, dw_mask]
            l = []
            for n in range(newNetworkList.shape[1]):
                if (G_tmp.has_edge(newNetworkList[0, n], newNetworkList[1, n])):
                    l.append((newNetworkList[0, n], newNetworkList[1, n]))
            G_tmp.remove_edges_from(l)
            has_path = nx.has_path(G_tmp, source, target)
            if ~np.any(has_path):
                flag = False
            else:
                j += 1
        if (j-1)<0:
            keep_source_w.append(False)
            #total_edges = total_edges + list(H.edges())
            #continue
        else:
            keep_source_w.append(True)
            G_tmp = H.copy()
            dw_mask = weightScore < thresholds_w[j-1]
            dw_mask = dw_mask.numpy()
            newNetworkList = networkList[:, dw_mask]
            l = []
            for n in range(newNetworkList.shape[1]):
                if (G_tmp.has_edge(newNetworkList[0, n], newNetworkList[1, n])):
                    l.append((newNetworkList[0, n], newNetworkList[1, n]))
            G_tmp.remove_edges_from(l)
        
            SignalingNet = G_tmp.copy()
            
            graphs = list(nx.connected_components(G_tmp.to_undirected()))
            graphsWithTF = []
            for graph in graphs:
                if TF_index_inNet in list(graph):
                    graphsWithTF.append(graph)
            largest_cc = max(graphsWithTF, key=len)
            G_tmp = G_tmp.subgraph(largest_cc)
            G_tmp = G_tmp.to_directed()
            nodesToKeepFinal = np.array(G_tmp.nodes)
            SignalingNet = nx.from_numpy_array(A.toarray(), parallel_edges=True, create_using=nx.DiGraph())
            SignalingNet = SignalingNet.edge_subgraph(list(G_tmp.edges)).copy()
            SignalingNet = nx.freeze(SignalingNet).copy()
            SignalingNet = nx.DiGraph(SignalingNet).copy()
            nodes_all = list(SignalingNet.nodes())
            remove_final = []
            for node in nodes_all:
                if node != target:
                    if node!=source:
                        if not nx.has_path(SignalingNet, node, target):
                            remove_final.append(node)
            SignalingNet.remove_nodes_from(remove_final)
            nodes_all = list(SignalingNet.nodes())        
            ### Remove nodes not accessed by the drug
            undrugable_nodes = []
            for node in nodes_all:
                if node!=source:
                    if node!=target:
                        paths_bool = []
                        if source in nodes_all:
                            paths_bool=nx.has_path(SignalingNet,source,node)
                        if paths_bool==False:
                            undrugable_nodes.append(node)
            if len(undrugable_nodes)>0:
                SignalingNet.remove_nodes_from(undrugable_nodes)
                nodes_all = list(SignalingNet.nodes())
            total_edges = total_edges + list(SignalingNet.edges())
            # print('Finished removing edges from %s to %s'%(drug_target_nodes[source_id],TF))
            # print(len(SignalingNet.nodes()))
    SignalingNet = nx.from_numpy_array(A.toarray(), parallel_edges=True, create_using=nx.DiGraph())
    total_edges = list(set(total_edges))
    SignalingNet = SignalingNet.edge_subgraph(total_edges).copy()
    SignalingNet = nx.freeze(SignalingNet).copy()
    SignalingNet = nx.DiGraph(SignalingNet).copy()
    nodes_all = list(SignalingNet.nodes())
    
    
    keep_source = list(np.array(keep_source)*np.array(keep_source_w))
    drug_target_nodes = drug_target_nodes[keep_source]
    sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    # targets_to_remove = drug_target_nodes[torch.abs(db[sources]).detach().numpy() <0.01]
    # sources_to_remove = [itemgetter(s)(map) for s in list(targets_to_remove)]
    # SignalingNet.remove_nodes_from(sources_to_remove)
    # nodes_all = list(SignalingNet.nodes())
    
    ### Remove nodes not accessed by the drug
    undrugable_nodes = []
    # sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    #print2log('Source drug-targets:%s'%list(drug_target_nodes))
    #print2log('TF of interest: %s'%TF)
    for node in nodes_all:
        if node not in sources:
            if node!=target:
                paths_bool = []
                for source in sources:
                    if source in nodes_all:
                        paths_bool.append(nx.has_path(SignalingNet,source,node))
                if not np.any(paths_bool):
                    undrugable_nodes.append(node)
    if len(undrugable_nodes)>0:
        SignalingNet.remove_nodes_from(undrugable_nodes)
        nodes_all = list(SignalingNet.nodes())
    keep_source = [] 
    for s in sources:
        if s not in nodes_all:
            keep_source.append(False)
        else:
            keep_source.append(True)
    drug_target_nodes = drug_target_nodes[keep_source]
    sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    # for s in sources:
    #     print(nx.has_path(SignalingNet,s,target))
    ### Remove nodes not accessing the TF
    uninteresting = []
    # sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    for node in nodes_all:
        if node!=target:
            paths_bool = nx.has_path(SignalingNet,node,target)
            if paths_bool==False:
                uninteresting.append(node)
    if len(uninteresting)>0:
        SignalingNet.remove_nodes_from(uninteresting)
    nodes_all = list(SignalingNet.nodes())
    keep_source = [] 
    for s in sources:
        if s not in nodes_all:
            keep_source.append(False)
        else:
            keep_source.append(True)
    drug_target_nodes = drug_target_nodes[keep_source]
    sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
    
    ### Save network
    # print2log('Is there a path %s'%nx.has_path(SignalingNet, sources[0], target))
    NetworkPandas = nx.to_pandas_edgelist(SignalingNet)
    if len(NetworkPandas)<=0:
        continue
    df_map = pd.DataFrame(map, index=[0]).T
    df_map = df_map.reset_index()
    df_map.columns = ["uniprot_names", "id"]
    df_annot = annotation.loc[:,['code','name']]
    df_annot.columns = ['uniprot_names','name']
    df_map = df_map.merge(df_annot, on='uniprot_names', how='left')
    df_map = df_map.loc[:,['id','name']]
    NetworkPandas = NetworkPandas.merge(df_map.rename({'id': 'source'}, axis=1), on='source', how='left')
    NetworkPandas = NetworkPandas.rename({'name': 'name_source'}, axis=1)
    NetworkPandas = NetworkPandas.merge(df_map.rename({'id': 'target'}, axis=1), on='target', how='left')
    NetworkPandas = NetworkPandas.rename({'name': 'name_target'}, axis=1)
    sign = lambda a: '-' if (a<0) else '+'
    NetworkPandas['interaction'] = [sign(w) for w in NetworkPandas["weight"].values]
    # NetworkPandas['interaction'] = np.sign(NetworkPandas["weight"])
    NetworkPandas.to_csv(ensembles_path+'MoA/'+drug_name+'/'+cell+'_'+drug_name+'_'+TF_gene+'_moa_model_%s.csv'%i)
    NetworkPandas['model_no'] = i
    node_type = []
    nodes_all = list(SignalingNet.nodes())
    for node in nodes_all:
        if node in sources:
            node_type.append('source')
        elif (node==target):
            node_type.append('target')
        else:
            node_type.append('mid_node')
    df_node_types = pd.DataFrame({'id':nodes_all,'type':node_type})
    df_node_types = df_node_types.merge(df_map, on='id', how='left')
    df_node_types = df_node_types.drop(['id'],axis=1)
    df_node_types.to_csv(ensembles_path+'MoA/'+drug_name+'/'+cell+'_'+drug_name+'_'+TF_gene+'_node_types_%s.csv'%i)
    
    if (i==0):
        df_all = NetworkPandas.copy()
        df_all_source_types = df_node_types[df_node_types['type']=='source'].reset_index(drop=True)
        df_all_target_types = df_node_types[df_node_types['type']=='target'].reset_index(drop=True)
        all_sources = pd.DataFrame({"sources":sources})
        all_sources["model"] = i
    else:
        df_all = pd.concat([df_all,NetworkPandas],0).reset_index(drop=True)
        df_all_source_types = pd.concat([df_all_source_types,
                                         df_node_types[df_node_types['type']=='source'].reset_index(drop=True)],0).reset_index(drop=True)
        df_all_target_types = pd.concat([df_all_target_types,
                                         df_node_types[df_node_types['type']=='target'].reset_index(drop=True)],0).reset_index(drop=True)
        all_sources_tmp = pd.DataFrame({"sources":sources})
        all_sources_tmp["model"] = i
        all_sources = pd.concat([all_sources,all_sources_tmp],0).reset_index(drop=True)
    ### Save shorter network
    max_shortest_source_paths = nx.shortest_path(SignalingNet, sources[0], target)
    for source in sources:
        tmp = nx.shortest_path(SignalingNet, source, target)
        if (len(tmp)>len(max_shortest_source_paths)):
            max_shortest_source_paths = nx.shortest_path(SignalingNet, source, target)
    shortpath = max_shortest_source_paths.copy()
    #shortpath = nx.shortest_path(SignalingNet, sources[0], target)
    shortNet = SignalingNet.subgraph(shortpath)
    NetworkPandas = nx.to_pandas_edgelist(shortNet)
    NetworkPandas = NetworkPandas.merge(df_map.rename({'id': 'source'}, axis=1), on='source', how='left')
    NetworkPandas = NetworkPandas.rename({'name': 'name_source'}, axis=1)
    NetworkPandas = NetworkPandas.merge(df_map.rename({'id': 'target'}, axis=1), on='target', how='left')
    NetworkPandas = NetworkPandas.rename({'name': 'name_target'}, axis=1)
    NetworkPandas['interaction'] = [sign(w) for w in NetworkPandas["weight"].values]
    NetworkPandas.to_csv(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_moa_short_model_%s.csv' % i)

    print2log('Connected %s' %nx.is_connected(SignalingNet.to_undirected()))
    print2log('Unique nodes %s'%len(SignalingNet.nodes()))
    print2log('Unique edges %s' % len(SignalingNet.edges()))
    print2log('Finished MoA for model %s'%i)

df_all.to_csv(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_all_models.csv')
mean_bias = mean_bias/numberOfModels
#%%
#df_all = pd.read_csv(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_all_models.csv',index_col=0)
# Ensemble by model
ensembled = df_all.groupby(['name_source','name_target']).agg(model_counts = pd.NamedAgg(column ='model_no', aggfunc=lambda x: x.nunique()))
all_sources_counted = all_sources['sources'].value_counts()
all_sources_counted = all_sources_counted[all_sources_counted>=0.5*numberOfModels]
all_sources_counted = all_sources_counted.reset_index()
all_sources_counted.columns = ['sources','counts']
df_map = pd.DataFrame(map, index=[0]).T
df_map = df_map.reset_index()
df_map.columns = ['sources_name','sources']
all_sources_counted = all_sources_counted.merge(df_map,on=['sources'],how='left')

df_all = df_all.iloc[:,[0,1,2,3,4,6]]
avg_weights = df_all.groupby(['name_source','name_target']).agg(mean_weight = pd.NamedAgg(column ='weight', aggfunc=np.mean))
df_all = df_all.merge(avg_weights,on=['name_source','name_target'],how='left')
df_all = df_all.iloc[:,[0,1,3,4,5,6]]
l = list(df_all.columns)
l[-1] = 'weight'
df_all.columns = l
df_all = df_all.drop_duplicates()
ensembled = ensembled.merge(df_all,on=['name_source','name_target'],how='left')
ensembled['interaction'] = [sign(w) for w in ensembled["weight"].values]
# ensembled = ensembled[ensembled['model_counts']>=0.3*numberOfModels]
edge_thresh = 0.5
ensembled_tmp = ensembled[ensembled['model_counts']>=edge_thresh*numberOfModels]
EnsembledNet = nx.from_pandas_edgelist(ensembled_tmp,source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
has_path = []
for s in all_sources_counted.sources:
    if s in list(EnsembledNet.nodes()):
        paths_bool = nx.has_path(EnsembledNet,s ,target)
        has_path.append(paths_bool)
    else:
        has_path.append(False)
while np.any(has_path)!=True:
    edge_thresh = edge_thresh - 0.05
    ensembled_tmp = ensembled[ensembled['model_counts']>edge_thresh*numberOfModels]
    EnsembledNet = nx.from_pandas_edgelist(ensembled_tmp,source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
    has_path = []
    for s in all_sources_counted.sources:
        if s in list(EnsembledNet.nodes()):
            paths_bool = nx.has_path(EnsembledNet,s,target)
            has_path.append(paths_bool)
        else:
            has_path.append(False)
ensembled = ensembled[ensembled['model_counts']>edge_thresh*numberOfModels]

EnsembledNet = nx.from_pandas_edgelist(ensembled,source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
G_tmp = EnsembledNet.copy()
graphs = list(nx.connected_components(G_tmp.to_undirected()))
graphsWithTF = []
for graph in graphs:
    if TF_index_inNet in list(graph):
        graphsWithTF.append(graph)
largest_cc = max(graphsWithTF, key=len)
G_tmp = G_tmp.subgraph(largest_cc)
G_tmp = G_tmp.to_directed()
nodesToKeepFinal = np.array(G_tmp.nodes)
SignalingNet = EnsembledNet.copy()
SignalingNet = SignalingNet.edge_subgraph(list(G_tmp.edges)).copy()
SignalingNet = nx.freeze(SignalingNet).copy()
SignalingNet = nx.DiGraph(SignalingNet).copy()
nodes_all = list(SignalingNet.nodes())

all_sources_counted = all_sources_counted[all_sources_counted.sources.isin(nodes_all)]
### Remove nodes not accessed by the drug
undrugable_nodes = []
# sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
for node in nodes_all:
    if node not in  list(all_sources_counted.sources):
        if node!=target:
            paths_bool = []
            for source in list(all_sources_counted.sources):
                if source in nodes_all:
                    paths_bool.append(nx.has_path(SignalingNet,source,node))
            if not np.any(paths_bool):
                undrugable_nodes.append(node)
if len(undrugable_nodes)>0:
    SignalingNet.remove_nodes_from(undrugable_nodes)
nodes_all = list(SignalingNet.nodes())

### Remove nodes not accessing the TF
uninteresting = []
# sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
for node in nodes_all:
    if node!=target:
        paths_bool = nx.has_path(SignalingNet,node,target)
        if paths_bool==False:
            uninteresting.append(node)
if len(uninteresting)>0:
    SignalingNet.remove_nodes_from(uninteresting)
nodes_all = list(SignalingNet.nodes())

#### Remove edges with frequency lower than 0.5 but without breaking connection from sources to TFs
ensembled  = nx.to_pandas_edgelist(SignalingNet)
ensembled = ensembled.merge(df_all.iloc[:,[0,1,2,3,5]].drop_duplicates(),on=['source','target','weight'],how='left')
ensembled = ensembled.drop_duplicates()
df_all = df_all.groupby(['name_source','name_target']).agg(model_counts = pd.NamedAgg(column ='model_no', aggfunc=lambda x: x.nunique()))
df_all = df_all.merge(ensembled,on=['name_source','name_target'],how='left')
df_all = df_all[~np.isnan(df_all['weight'])]
df_all = df_all.iloc[:,[2,3,4]].drop_duplicates()
ensembled = ensembled.merge(df_all,on=['source','target'],how='left')
ensembled_tmp = ensembled[ensembled['model_counts']>0.5*numberOfModels]
EnsembledNet = nx.from_pandas_edgelist(ensembled_tmp,source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
has_path = []
for s in all_sources_counted.sources:
    if s in list(EnsembledNet.nodes()):
        paths_bool = nx.has_path(EnsembledNet,s,target)
        has_path.append(paths_bool)
    else:
        has_path.append(False)
if np.any(has_path):
    G_tmp = EnsembledNet.copy()
    graphs = list(nx.connected_components(G_tmp.to_undirected()))
    graphsWithTF = []
    for graph in graphs:
        if TF_index_inNet in list(graph):
            graphsWithTF.append(graph)
    largest_cc = max(graphsWithTF, key=len)
    G_tmp = G_tmp.subgraph(largest_cc)
    G_tmp = G_tmp.to_directed()
    nodesToKeepFinal = np.array(G_tmp.nodes)
    SignalingNet = EnsembledNet.copy()
    SignalingNet = SignalingNet.edge_subgraph(list(G_tmp.edges)).copy()
    SignalingNet = nx.freeze(SignalingNet).copy()
    SignalingNet = nx.DiGraph(SignalingNet).copy()
    nodes_all = list(SignalingNet.nodes())
else:
    ### Find edges required to connect one of the targets with the TF and keep the ones with the highest sum of frequencies
    all_paths = []
    all_paths_df = pd.DataFrame({'source':[],'target':[],'weight':[],'path':[]})
    EnsembledNet_larger = nx.from_pandas_edgelist(ensembled,source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
    p_ind = 0
    for s in all_sources_counted.sources:
        paths = pathFinder(s,target,EnsembledNet_larger,mode='All',max_depth=7)
        all_paths = all_paths + paths
        for i,p in enumerate(paths):
            subnet = EnsembledNet_larger.subgraph(p)
            subnet_df = nx.to_pandas_edgelist(subnet)
            subnet_df['path'] = p_ind
            all_paths_df = all_paths_df.append(subnet_df)
            p_ind +=1
    all_paths_df = all_paths_df.merge(ensembled,on=['source','target','weight'],how='left')
    avg_paths = all_paths_df.groupby(['path']).agg(avg_counts = pd.NamedAgg(column ='model_counts', aggfunc=np.mean))
    final_path = avg_paths.index[np.argmax(avg_paths.avg_counts.values)]
    final_path_df = all_paths_df[all_paths_df['path']==final_path].loc[:,ensembled_tmp.columns.values]
    ensembled_tmp = ensembled_tmp.append(final_path_df)
    ensembled_tmp = ensembled_tmp.drop_duplicates() 
    EnsembledNet = nx.from_pandas_edgelist(ensembled_tmp,source='source', target='target', edge_attr='weight', create_using=nx.DiGraph())
    G_tmp = EnsembledNet.copy()
    graphs = list(nx.connected_components(G_tmp.to_undirected()))
    graphsWithTF = []
    for graph in graphs:
        if TF_index_inNet in list(graph):
            graphsWithTF.append(graph)
    largest_cc = max(graphsWithTF, key=len)
    G_tmp = G_tmp.subgraph(largest_cc)
    G_tmp = G_tmp.to_directed()
    nodesToKeepFinal = np.array(G_tmp.nodes)
    SignalingNet = EnsembledNet.copy()
    SignalingNet = SignalingNet.edge_subgraph(list(G_tmp.edges)).copy()
    SignalingNet = nx.freeze(SignalingNet).copy()
    SignalingNet = nx.DiGraph(SignalingNet).copy()
    nodes_all = list(SignalingNet.nodes())
### Get final networked
all_sources_counted = all_sources_counted[all_sources_counted.sources.isin(nodes_all)]
### Remove nodes not accessed by the drug
undrugable_nodes = []
# sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
for node in nodes_all:
    if node not in  list(all_sources_counted.sources):
        if node!=target:
            paths_bool = []
            for source in list(all_sources_counted.sources):
                if source in nodes_all:
                    paths_bool.append(nx.has_path(SignalingNet,source,node))
            if not np.any(paths_bool):
                undrugable_nodes.append(node)
if len(undrugable_nodes)>0:
    SignalingNet.remove_nodes_from(undrugable_nodes)
nodes_all = list(SignalingNet.nodes())

### Remove nodes not accessing the TF
uninteresting = []
# sources = [itemgetter(s)(map) for s in list(drug_target_nodes)]
for node in nodes_all:
    if node!=target:
        paths_bool = nx.has_path(SignalingNet,node,target)
        if paths_bool==False:
            uninteresting.append(node)
if len(uninteresting)>0:
    SignalingNet.remove_nodes_from(uninteresting)
nodes_all = list(SignalingNet.nodes())
    
#### Save final network
NetworkPandas = nx.to_pandas_edgelist(SignalingNet)
df_map = pd.DataFrame(map, index=[0]).T
df_map = df_map.reset_index()
df_map.columns = ["uniprot_names", "id"]
df_annot = annotation.loc[:,['code','name']]
df_annot.columns = ['uniprot_names','name']
df_map = df_map.merge(df_annot, on='uniprot_names', how='left')
df_map = df_map.loc[:,['id','name']]
NetworkPandas = NetworkPandas.merge(df_map.rename({'id': 'source'}, axis=1), on='source', how='left')
NetworkPandas = NetworkPandas.rename({'name': 'name_source'}, axis=1)
NetworkPandas = NetworkPandas.merge(df_map.rename({'id': 'target'}, axis=1), on='target', how='left')
NetworkPandas = NetworkPandas.rename({'name': 'name_target'}, axis=1)
sign = lambda a: '-' if (a<0) else '+'
NetworkPandas['interaction'] = [sign(w) for w in NetworkPandas["weight"].values]
NetworkPandas.to_csv(ensembles_path+'MoA/'+drug_name+'/'+cell+'_'+drug_name+'_'+TF_gene+'_moa_model_ensembleFiltered.csv')

node_type = []
nodes_all = list(SignalingNet.nodes())
names = []
df_annot.columns = ['sources_name','name']
all_sources_counted = all_sources_counted.merge(df_annot,how='left')
for node in nodes_all:
    node_name = df_map[df_map['id']==node].name.values[0]
    # if (node_name in list(df_all_source_types.name)):
    if (node_name in list(all_sources_counted.name)):
        node_type.append('source')
    elif (node_name in list(df_all_target_types.name.unique())):
        node_type.append('target')
    else:
        node_type.append('mid_node')
    names.append(node_name)
df_node_types = pd.DataFrame({'name':names,'type':node_type})
df_node_types.to_csv(ensembles_path+'MoA/'+drug_name+'/'+cell+'_'+drug_name+'_'+TF_gene+'_node_types_ensembleFiltered.csv')

print2log('The whole trimmed network')
print2log(NetworkPandas)

### Save shortest path net
NetworkPandas_short = pd.DataFrame({'source':[],'target':[],'weight':[]})
for source in  all_sources_counted.sources:
    shortest_path = nx.shortest_path(SignalingNet, source, target)
    shortNet = SignalingNet.subgraph(shortest_path)
    subnet_df = nx.to_pandas_edgelist(shortNet)
    NetworkPandas_short = NetworkPandas_short.append(subnet_df)
NetworkPandas_short = NetworkPandas_short.merge(df_map.rename({'id': 'source'}, axis=1), on='source', how='left')
NetworkPandas_short = NetworkPandas_short.rename({'name': 'name_source'}, axis=1)
NetworkPandas_short = NetworkPandas_short.merge(df_map.rename({'id': 'target'}, axis=1), on='target', how='left')
NetworkPandas_short = NetworkPandas_short.rename({'name': 'name_target'}, axis=1)
sign = lambda a: '-' if (a<0) else '+'
NetworkPandas_short['interaction'] = [sign(w) for w in NetworkPandas_short["weight"].values]
NetworkPandas_short.to_csv(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_moa_short_model_ensembled.csv')
print2log('The shortest path trimmed network')
print2log(NetworkPandas_short)

### Save path of largest average frequency
NetworkPandas_frequent = pd.DataFrame({'source':[],'target':[],'path':[]})
all_paths = []
p_ind = 0
for s in all_sources_counted.sources:
    paths = pathFinder(s,target,SignalingNet,mode='All',max_depth=7)
    all_paths = all_paths + paths
    for i,p in enumerate(paths):
        subnet = SignalingNet.subgraph(p)
        subnet_df = nx.to_pandas_edgelist(subnet)
        subnet_df['path'] = p_ind
        NetworkPandas_frequent = NetworkPandas_frequent.append(subnet_df)
        p_ind +=1
NetworkPandas_frequent = NetworkPandas_frequent.merge(ensembled,on=['source','target','weight'],how='left')
avg_paths = NetworkPandas_frequent.groupby(['path']).agg(avg_counts = pd.NamedAgg(column ='model_counts', aggfunc=np.mean))
final_path = avg_paths.index[np.argmax(avg_paths.avg_counts.values)]
NetworkPandas_frequent = NetworkPandas_frequent[NetworkPandas_frequent['path']==final_path]
sign = lambda a: '-' if (a<0) else '+'
NetworkPandas_frequent['interaction'] = [sign(w) for w in NetworkPandas_frequent["weight"].values]
NetworkPandas_frequent.to_csv(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_moa_frequent_model_ensembled.csv')
print2log('The most frequent trimmed network')
print2log(NetworkPandas_frequent)

# Plot the most frequent path
subnet = nx.from_pandas_edgelist(NetworkPandas_frequent,source='name_source', target='name_target', edge_attr='weight', create_using=nx.DiGraph())
fig = plt.figure(num=None, figsize=(8, 8))
ax = fig.add_subplot(111)
pos = nx.spring_layout(subnet)
nx.draw_networkx_nodes(subnet, pos, node_size=500, ax=ax)
nx.draw_networkx_labels(subnet, pos,font_size=8, ax=ax)
for edge in subnet.edges(data=True):
	sign = edge[2]['weight']
	if sign < 0:
		sign = '-['
	else:
		sign = '-|>'
	nx.draw_networkx_edges(subnet, pos, edgelist=[(edge[0], edge[1])],
                        node_size=500,
                        arrowstyle=sign,
                        ax=ax)
# Turn off unnecessary x and y axis labels
ax.set_axis_off()
# Save figure with a tight bbox to ensure it isn't cut off
fig.savefig(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_moa_frequent_model_ensembled.png', 
            bbox_inches='tight', dpi=600)

# Plot the most shortest path
subnet = nx.from_pandas_edgelist(NetworkPandas_short,source='name_source', target='name_target', edge_attr='weight', create_using=nx.DiGraph())
fig = plt.figure(num=None, figsize=(8, 8))
ax = fig.add_subplot(111)
pos = nx.spring_layout(subnet)
nx.draw_networkx_nodes(subnet, pos, node_size=500, ax=ax)
nx.draw_networkx_labels(subnet, pos,font_size=8, ax=ax)
for edge in subnet.edges(data=True):
	sign = edge[2]['weight']
	if sign < 0:
		sign = '-['
	else:
		sign = '-|>'
	nx.draw_networkx_edges(subnet, pos, edgelist=[(edge[0], edge[1])],
                        node_size=500,
                        arrowstyle=sign,
                        ax=ax)
# Turn off unnecessary x and y axis labels
ax.set_axis_off()
# Save figure with a tight bbox to ensure it isn't cut off
fig.savefig(ensembles_path + 'MoA/'+drug_name+'/'+ cell +'_'+drug_name+'_'+TF_gene+ '_moa_short_model_ensembled.png', 
            bbox_inches='tight', dpi=600)

print2log(all_sources_counted)