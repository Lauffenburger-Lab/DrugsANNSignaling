import pandas as pd
import numpy
import argparse

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Extract receptor-ligands from the PKN')
parser.add_argument('--species', action='store', help='species_id',default=9606)
parser.add_argument('--add', action='store',help='interactions to manually add',default='preprocessed_data/RL/add.tsv')
parser.add_argument('--remove', action='store',help='interactions to manually remove',default='preprocessed_data/RL/remove.tsv')
parser.add_argument('--edit', action='store',help='interactions to manually edit',default='preprocessed_data/RL/edit.tsv')
parser.add_argument('--RLFull', help='all receptor-ligand interactions in .tsv format', default = 'preprocessed_data/PKN/RLFull.tsv')
parser.add_argument('--RL', help='receptors-ligands in .tsv format filtered', default='preprocessed_data/PKN/RL.tsv')
parser.add_argument('--WholePKN', help='all prior knowledge interactions tha', default='../data/omnipath_webservice_interactions__recent.tsv')
args = parser.parse_args()
species = int(args.species)
add = args.add
remove = args.remove
edit = args.edit
RLFull = args.RLFull
RL = args.RL
WholePKN = args.WholePKN

# Interaction sources to be kept
trustedSource = numpy.array(['KEGG',
              'InnateDB',
              'Ramilowski2015',
              'Baccin2019',
              'Reactome_LRdb',
              'UniProt_LRdb',
              'CellPhoneDB',
              'HuGeSiM'
              ])


# Load omnipath and keep only species of interest
omnipath = pd.read_csv(WholePKN, sep='\t', low_memory=False)
humanFilter = omnipath['ncbi_tax_id_target'] == species
omnipath = omnipath.loc[humanFilter, :]

#Only ligand-receptor interactions
LRFilter = omnipath['ligrecextra'].values == True
omnipath =  omnipath.loc[LRFilter, :]

#Only in omnipath
omnipathFilter = omnipath['omnipath'].values
omnipath =  omnipath.loc[omnipathFilter, :]

#Subset to only uniprot proteins
uniprot = pd.read_csv('../data/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab', sep='\t', low_memory=False)
uniprot = uniprot['Entry'].values
uniprotFilter = numpy.logical_and(numpy.isin(omnipath['source'].values, uniprot), numpy.isin(omnipath['target'].values, uniprot))
omnipath = omnipath.loc[uniprotFilter, :]

#Subset to relevant info
relevantInformation = ['source', 'target', 'consensus_direction',  'consensus_stimulation', 'consensus_inhibition', 'sources', 'references']
omnipath = omnipath[relevantInformation]
omnipath = omnipath.rename(columns={'consensus_direction': 'direction', 'consensus_stimulation': 'stimulation', 'consensus_inhibition': 'inhibition'})
omnipath[['references']] = omnipath[['references']].astype(str)

#Remove interactions without reference
referenceFilter = omnipath['references']=='nan'
omnipath = omnipath.loc[referenceFilter==False, :]

#Resolve paradoxes
paradox = numpy.logical_and(omnipath['stimulation'].values==1, omnipath['inhibition'].values==1)
omnipath.loc[paradox, ['stimulation', 'inhibition']] = 0

#Add interactions manually
currationAdd = pd.read_csv(add, sep='\t', low_memory=False)
for i in range(currationAdd.shape[0]):
    curSource = currationAdd.iloc[i,:]['source']
    curTarget = currationAdd.iloc[i,:]['target']
    inList = numpy.logical_and(numpy.isin(omnipath['source'], curSource), numpy.isin(omnipath['target'], curTarget))
    if sum(inList) == 0:  #add new
        omnipath = omnipath.append(currationAdd.iloc[i,:])
    else:         #add reference
        #Note, does not check for consistency of sign etc
        omnipath.loc[inList,'sources'] = omnipath.loc[inList,'sources'] + ';' + currationAdd.iloc[i,:]['sources']
        omnipath.loc[inList,'references'] = omnipath.loc[inList,'references'] + ';' + currationAdd.iloc[i,:]['references']

#Remove interactions manually
currationRemove = pd.read_csv(remove, sep='\t', low_memory=False)
for i in range(currationRemove.shape[0]):
    curSource = currationRemove.iloc[i,:]['source']
    curTarget = currationRemove.iloc[i,:]['target']
    inList = numpy.logical_and(numpy.isin(omnipath['source'], curSource), numpy.isin(omnipath['target'], curTarget))
    if sum(inList)>0:
        omnipath = omnipath.loc[inList==False,:]
    else:
        print('No match for remove', currationRemove.iloc[i,:])


#Edit interactions manually
currationEdit = pd.read_csv(edit, sep='\t', low_memory=False)
for i in range(currationEdit.shape[0]):
    curSource = currationEdit.iloc[i,:]['source']
    curTarget = currationEdit.iloc[i,:]['target']
    inList = numpy.logical_and(numpy.isin(omnipath['source'], curSource), numpy.isin(omnipath['target'], curTarget))
    #Note, should add a reference for the change
    if sum(inList)>0:
        if currationEdit.iloc[i,:]['action'] == 'set_stimulation':
            omnipath.loc[inList,'stimulation'] = currationEdit.iloc[i,:]['value']
        elif currationEdit.iloc[i,:]['action'] == 'set_inhibition':
            omnipath.loc[inList,'inhibition'] = currationEdit.iloc[i,:]['value']
        elif currationEdit.iloc[i,:]['action'] == 'reverse_direction':
            omnipath.loc[inList,'source'] = currationEdit.iloc[i,:]['target']
            omnipath.loc[inList,'target'] = currationEdit.iloc[i,:]['source']
            omnipath.loc[inList,'direction'] = 1 #for it to be meaningful to reverse the direction, the direction must be known
        elif currationEdit.iloc[i,:]['action'] == 'set_direction':
            omnipath.loc[inList,'direction'] = currationEdit.iloc[i,:]['value']
    else:
        print('No match for edit', currationEdit.iloc[i,:])


#Remove interactions with same source and target
sameSourceAndTargetFilter = omnipath['source'] == omnipath['target']
print('Removed interactions with same source and target', sum(sameSourceAndTargetFilter))
omnipath =  omnipath.loc[sameSourceAndTargetFilter==False, :]


#Duplicate reversible
reversibleFilter = omnipath['direction'] == 0
revOmni = omnipath.loc[reversibleFilter,:].copy()
revOmni = revOmni.rename(columns={'source': 'target', 'target': 'source'})
omnipath = pd.concat([omnipath.copy(), revOmni])
omnipath['direction'] = 1

# Save all receptor-ligand interactions
omnipath.to_csv(RLFull, sep='\t', index=False)

#Keep only trusted sources
pknFilter = numpy.full(omnipath.shape[0], False, dtype=bool)
for i in range(len(pknFilter)):
    pknFilter[i] = len(numpy.intersect1d(omnipath.iloc[i,:]['sources'].split(';'), trustedSource))>0
omnipath =  omnipath.loc[pknFilter, :]
# Save all receptor-ligand interactions that are in specific sources
omnipath.to_csv(RL, sep='\t', index=False)

