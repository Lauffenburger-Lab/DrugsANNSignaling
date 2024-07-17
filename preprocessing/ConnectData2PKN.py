import pandas as pd
import numpy as np
import argparse

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Connect data to PKN')
parser.add_argument('--pknPath', action='store', help='untrimmed PKN file in .tsv format',default='preprocessed_data/PKN/pkn.tsv')
parser.add_argument('--DTIpath', action='store',help='untrimmed drug-target interactions in long format saved in .tsv format',default='preprocessed_data/PKN/L1000_lvl3_DT.tsv')
parser.add_argument('--TFpath', action='store', help='untrimmed TF activity file in .tsv format', required=True)
parser.add_argument('--targetedTFs', help='path to save TFs that are directly targeted by a drug in .tsv format', required=True)
parser.add_argument('--forced2keep', help='path to save PKN parts to forcefuly keep because it contains targeted TFs', required=True)
parser.add_argument('--outTrimmedTFs', help='path to re-save TF activities, trimmed to include only TFs in the PKN', required=True)
parser.add_argument('--outTrimmedDTIs', help='path to re-save DTIs, trimmed to include only drugs with targets in the PKN', required=True)
args = parser.parse_args()
pknPath = args.pknPath
DTIpath = args.DTIpath
TFpath = args.TFpath
targetedTFs = args.targetedTFs
forced2keep = args.forced2keep
outTrimmedTFs = args.outTrimmedTFs
outTrimmedDTIs = args.outTrimmedDTIs

### Load files
pkn = pd.read_csv(pknPath, sep='\t', low_memory=False)
nodeNames = np.union1d(pkn['source'], pkn['target'])
DTI = pd.read_csv(DTIpath, sep='\t', low_memory=False)
TFactivities = pd.read_csv(TFpath, sep='\t', low_memory=False, index_col=0)

### Filter data
DTI = DTI[DTI['target'].isin(nodeNames)]
keptTFs = np.intersect1d(TFactivities.columns.values, nodeNames) 
TFactivities = TFactivities.loc[:,keptTFs]
tfs_targeted = np.intersect1d(TFactivities.columns.values,DTI.target.values)    
boolean_index = np.logical_or(pkn['target'].isin(tfs_targeted), pkn['source'].isin(tfs_targeted))
pkn_tfs_targeted = pkn[boolean_index]

### Save files
pkn_tfs_targeted.reset_index(drop=True).to_csv(forced2keep, sep='\t', index=True)
pd.DataFrame({'Entry':tfs_targeted}).reset_index(drop=True).to_csv(targetedTFs, sep='\t', index=True)
TFactivities.to_csv(outTrimmedTFs, sep='\t', index=True)
DTI.to_csv(outTrimmedDTIs, sep='\t', index=False)