import numpy as np
import pandas as pd
import pyreadr
import argparse

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Convert CSV or TSV file data frame to RDS file')
parser.add_argument('--type', action='store',help='only TSV or CSV')
parser.add_argument('--full_file_path', action='store',help='The CSV or TSV input file name (with the extention .csv or .tsv) file and its full path to convert into RDS')
parser.add_argument('--output_full_file_path', action='store',help='The RDS (with the extention .rds) and its full path to convert from CSV/TSV')
args = parser.parse_args()
full_file_path= args.full_file_path
output_full_file_path= args.output_full_file_path

### Load the data
if type=='CSV':
    df = pd.read_csv(full_file_path,index_col = 0 )
else:
   df = pd.read_csv(full_file_path,sep='\t',index_col = 0)

df = df.reset_index()
# Rename the new column if needed (default will be 'index')
df.rename(columns={'index': 'rownames'}, inplace=True)

### Save the file
pyreadr.write_rds(output_full_file_path, df)