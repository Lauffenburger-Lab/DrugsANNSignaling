import numpy as np
import pandas as pd
import pyreadr
import argparse

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Convert a TSV or CSV file from matrix to long format and vice versa')
parser.add_argument('--type', action='store',help='only TSV or CSV')
parser.add_argument('--full_file_path', action='store',help='The CSV or TSV input file name (with the extention .csv or .tsv) file and its full path')
parser.add_argument('--output_full_file_path', action='store',help='The full path to save the converted file')
parser.add_argument('--input_file_format', action='store',help='The original format of the file (one of matrix or long)')
args = parser.parse_args()
full_file_path= args.full_file_path
output_full_file_path= args.output_full_file_path
input_file_format= args.input_file_format

save_index = False
### Load the data
if type=='CSV':
    df = pd.read_csv(full_file_path,index_col = 0 )
else:
   df = pd.read_csv(full_file_path,sep='\t',index_col = 0)

if input_file_format=='matrix':
    print('Converting matrix to long format')
    # Reset index to move row names into a column
    df = df.reset_index()
    # Melt the dataframe
    df = pd.melt(df, id_vars=['index'], var_name='target', value_name='Value')
    # Rename 'index' to something meaningful, e.g., 'Row'
    df.rename(columns={'index': 'source'}, inplace=True)
    df = df[df['Value']==1]
    df = df.loc[:,['source','target']]
elif input_file_format=='long':
    df = df.reset_index()
    df['Value'] = 1
    df.columns = ['source','target','Value']
    # Convert back to the original wide format
    df = df.pivot(index='source', columns='target', values='Value')
    df = df.fillna(0)
    save_index = True

if type=='CSV':
    df.to_csv(output_full_file_path, index=save_index)
else:
    df.to_csv(output_full_file_path,sep='\t', index=save_index)