import pandas as pd
import numpy
import networkx

def getConnectionToTF(model, allTFs, affectedNodes):
    g = networkx.from_pandas_edgelist(model, 'source', 'target', create_using=networkx.DiGraph())
    allNodes = numpy.array(list(g.nodes))
    includedTFs = numpy.intersect1d(allTFs, allNodes)

    connectedToTF = numpy.isin(affectedNodes, includedTFs)
    for i in range(len(affectedNodes)):
        if affectedNodes[i] in allNodes:
            for tf in includedTFs:
                if networkx.algorithms.shortest_paths.generic.has_path(g, affectedNodes[i], tf):
                    connectedToTF[i] = True
                    break

    return connectedToTF

def getConnectionToLigand(model, allLigands, affectedNodes):
    g = networkx.from_pandas_edgelist(model, 'target', 'source', create_using=networkx.DiGraph())
    allNodes = numpy.array(list(g.nodes))
    includedLigands = numpy.intersect1d(allLigands, allNodes)
    connectedToLigand = numpy.isin(affectedNodes, allLigands)
    for i in range(len(affectedNodes)):
        if affectedNodes[i] in allNodes:
            for ligand in includedLigands:
                if networkx.algorithms.shortest_paths.generic.has_path(g, affectedNodes[i], ligand):
                    connectedToLigand[i] = True
                    break
    return connectedToLigand

def trimDeadEnds(model, allTFs, allLigands):
    allNodes = numpy.union1d(model['source'].values, model['target'].values)

    connectedToTF = getConnectionToTF(model, allTFs, allNodes)
    connectedToLigand = getConnectionToLigand(model, allLigands, allNodes)
    connectedToBoth = numpy.logical_and(connectedToTF, connectedToLigand)
    disconectedNodes = allNodes[connectedToBoth == False]

    dissconectedSources = numpy.isin(model.source.values, disconectedNodes)
    disconnectedTargets = numpy.isin(model.target.values, disconectedNodes)
    disconnectedEdges = numpy.logical_or(dissconectedSources, disconnectedTargets)
    model = model.loc[disconnectedEdges==False,:]
    return model

def trimSelfConnecting(model, allTFs, allLigands):
    lastSize = numpy.inf
    curentSize = model.shape[0]
    while curentSize<lastSize:
        lastSize = curentSize
        sources, counts = numpy.unique(model['source'], return_counts=True)
        onlyOneInput = sources[counts == 1]
        targets, counts = numpy.unique(model['target'], return_counts=True)
        onlyOneOutput = targets[counts == 1]
        overlap = numpy.intersect1d(onlyOneInput, onlyOneOutput)
        #Exclude ligands and TFs
        overlap = numpy.setdiff1d(overlap, allLigands)
        overlap = numpy.setdiff1d(overlap, allTFs)
        selfLoop = numpy.full(len(overlap), False, dtype=bool)
        for i in range(len(selfLoop)):
            curSource = model.loc[numpy.isin(model['target'], overlap[i]), 'source'].values[0]
            curTarget = model.loc[numpy.isin(model['source'], overlap[i]), 'target'].values[0]
            selfLoop[i] = curSource == curTarget

        affectedProteins = overlap[selfLoop]

        affectedInteractions = numpy.logical_or(numpy.isin(model['source'], affectedProteins), numpy.isin(model['target'], affectedProteins))
        model = model.loc[affectedInteractions==False,:]
        curentSize = model.shape[0]

    return model

def subsetOnSource(df, coreSources):
    dfFilter = numpy.full(df.shape[0], False, dtype=bool)
    for i in range(len(dfFilter)):
        dfFilter[i] = len(numpy.intersect1d(df.iloc[i,:]['sources'].split(';'), coreSources))>0
    df = df.loc[dfFilter,:].copy()
    return df


coreSources = ['KEGG', 'InnateDB','SIGNOR']


#Load and trim PKN
PKN = pd.read_csv('preprocessed_data/PKN/pkn.tsv', sep='\t', low_memory=False)
PKNFull = PKN.copy()
PKN = subsetOnSource(PKN, coreSources)

DT = pd.read_csv('preprocessed_data/PKN/L1000_lvl3_DT.tsv', sep='\t', low_memory=False)
inferedReceptors = numpy.intersect1d(DT['target'].values, PKNFull['source'].values)

TFgene = pd.read_csv('preprocessed_data/TF_activities/Trimmed_l1000_allgenes_lvl3_tfs.tsv ', sep='\t', low_memory=False, index_col=0)
allTFs = numpy.intersect1d(TFgene.columns.values, PKN['target'].values)

PKN = trimDeadEnds(PKN, allTFs, inferedReceptors)
PKN = trimSelfConnecting(PKN, allTFs, inferedReceptors)
targetd_tfs_interactions = pd.read_csv('preprocessed_data/PKN/L1000_latest_Add_lvl3.tsv',sep='\t', low_memory=False, index_col=0)
targeted_tfs = pd.read_csv('preprocessed_data/TF_activities/tfs_targetd_alls_genes_lvl3.tsv',sep='\t', low_memory=False, index_col=0)
PKN = pd.concat((PKN,targetd_tfs_interactions))

allTFs = numpy.intersect1d(allTFs, numpy.union1d(PKN['target'],targeted_tfs['Entry']))


PKN.to_csv('preprocessed_data/PKN/l1000_lvl3_withsignor-Model.tsv', sep='\t', index=False)


#Build annotation file
uniprot = pd.read_csv('../data/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t', low_memory=False)
uniprot = uniprot.loc[:, ['Entry', 'Gene names  (primary )']].values
nodeNames = numpy.union1d(PKN['source'], PKN['target'])
annotation = uniprot[numpy.isin(uniprot[:,0], nodeNames),:]
missingNodes = nodeNames[numpy.isin(nodeNames, uniprot[:,0])==False]
missingNodeAnotation = numpy.array([missingNodes, missingNodes]).T
annotation = numpy.concatenate((annotation, missingNodeAnotation))
annotation = pd.DataFrame(annotation, columns=['code', 'name'])
annotation = annotation.drop_duplicates(subset='code', keep='first')
annotation['name'] = annotation['name'].str.replace('; ','/')
annotation['TF'] = numpy.isin(annotation['code'], allTFs)
annotation['ligand'] = False
annotation.to_csv('preprocessed_data/PKN/l1000_lvl3_withsignor-Annotation.tsv', sep='\t', index=False)

