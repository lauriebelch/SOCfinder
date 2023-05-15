#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 09:34:22 2023

@author: laurie
"""

## required python packages
import pandas as pd
import numpy as np
import argparse
import os
import sys
import re
import glob

## required custom function
def belch7(x, y, z, a, b, c, d):
    return exec((x + y + z + a + b + c + d))

def belch5(x, y, z, a, b):
    return exec((x + y + z + a + b))

def belch3(x, y, z):
    return exec((x + y + z))

## Argument parsing

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='SOCfinder parser')

# Add command-line arguments with flags and help messages
parser._optionals.title = 'Required arguments'
parser.add_argument('-ki', '--KOFAMinput', type=str, metavar='.', required=True, help='Path to KOFAM input')
parser.add_argument('-b', '--BLASTinput', type=str, metavar='.', required=True, help='Path to BLAST input folder')
parser.add_argument('-ad', '--ANTISMASHRegionGBK', metavar='.', type=str, required=True, help='Path to directory of region_gbk files')
parser.add_argument('-ac', '--accession', type=str, metavar='.', required=True, help='Accession number of the genome')
parser.add_argument('-k', '--KO', type=str, required=True, metavar='.', help='Path to Social KO file')
parser.add_argument('-a', '--ANTISMASHtypes', type=str, required=True, metavar='.', help='Path to list of ANTISMASH types')
parser.add_argument('-so', '--SOCKfile', type=str, required=True, metavar='.', help='Name of final list of cooperative genes (csv)')
parser.add_argument('-ko', '--KOFAMoutput', type=str, required=True, metavar='.', help='Name of KOFAM output csv')
parser.add_argument('-ao', '--ANTISMASHoutput', type=str, required=True, metavar='.', help='Name of ANTISMASH output csv')
parser.add_argument('-bo', '--BLASToutput', type=str, required=True, metavar='.', help='Name of BLAST output csv')
#parser.add_argument('-s', '--score', type=float, help='Numerical score')
#parser.add_argument('-y', '--helloworld', action='store_true', help='Print "Hello World"')

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the command-line arguments
inputpath = args.KOFAMinput
outputpath = args.KOFAMoutput
kpath = args.KO
antismashtypes = args.ANTISMASHtypes
antismashoutputpath = args.ANTISMASHoutput
accession = args.accession
antismashdirectory = args.ANTISMASHRegionGBK
blastinputfolder = args.BLASTinput
blastoutputpath = args.BLASToutput
SOCKSS = args.SOCKfile
#score = args.score

# Code if we had an option like y that we can turn on or off
#if args.helloworld:
 #   print("Hello World")

# store directory that script is running in
script_path = os.path.dirname(os.path.abspath(sys.argv[0]))

######################################################################
################### KOFAM ############################################
######################################################################

# load social KO
sock = pd.read_csv(kpath)
# define input file
file = inputpath
# define output file
outfile = outputpath

# read in data and tidy columns
data = pd.read_csv(file, sep="\t", quoting=3, header=None, names=["#", "gene name", "KO", "thrshld", "score", "E-value", "KO definition"])

# extract significant hits
data = data.loc[data["#"] == "*", :]

# count social genes in kofamscan output
soc = data.loc[data["KO"].isin(sock["term"]), :]
print(len(soc["gene name"].unique())) # 103 bacillus social genes

# remove duplicates
socks = soc.iloc[:, 1].unique()

# save output
pd.DataFrame(socks).to_csv(outfile, index=False, header=True)

######################################################################
################### BLASTE ###########################################
######################################################################

# set working directory
os.chdir(blastinputfolder)

# define file paths
# main blast output
files = glob.glob("*_PSORT.txt")
if len(files) > 0:
    blaste_output_filename = files[0]
else:
    print("No _PSORT.txt file found")
# blaste to psortbC extracellular filename
files1 = glob.glob("*_PSORT_E.txt")
if len(files1) > 0:
    blaste_e_output_filename = files1[0]
else:
    print("No _PSORT_E.txt file found")
# blaste to psortbC non-extracellular filename
files2 = glob.glob("*_PSORT_NE.txt")
if len(files2) > 0:
    blaste_ne_output_filename = files2[0]
else:
    print("No _PSORT_NE.txt file found")


# name of output file
outfile = blastoutputpath


# load blast output for PSORTB COMPUTATIONAL EXTRACELLULAR
data = pd.read_table(blaste_output_filename, header=None, delim_whitespace=True)
# tidy column names
data = data.assign(match_length=data.iloc[:, 7] - data.iloc[:, 6] + 1)
data.columns = ["subject_seq_id", "query_acc", "query_length", "evalue",
                "bitscore", "subject_start", "subject_end", "subject_length", "match_length"]
# load blast output for PSORTB EXPERIMENTAL EXTRACELLULAR
data_e = pd.read_table(blaste_e_output_filename, header=None, delim_whitespace=True)
data_e = data_e.assign(match_length=data_e.iloc[:, 7] - data_e.iloc[:, 6] + 1)
data_e.columns = ["subject_seq_id", "query_acc", "query_length", "evalue",
                  "bitscore", "subject_start", "subject_end", "subject_length", "match_length"]
# load blast output for PSORTB EXPERIMENTAL NON-EXTRACELLULAR
data_ne = pd.read_table(blaste_ne_output_filename, header=None, delim_whitespace=True)
data_ne = data_ne.assign(match_length=data_ne.iloc[:, 7] - data_ne.iloc[:, 6] + 1)
data_ne.columns = ["subject_seq_id", "query_acc", "query_length", "evalue",
                   "bitscore", "subject_start", "subject_end", "subject_length", "match_length"]

SOCK = [] # empty list to store cooperative genes

## 1)  add to social if EXACT match in experimental extracellular
temp = np.where(data_e['evalue']==0)[0] # temp stores row ID of target
# if there is at least one target, we add it to the list
if len(temp)>0:
    SOCK = np.append(SOCK, np.unique(data_e.loc[temp, 'query_acc']))

#### 2) remove from consideration if EXACT match in experimental non-extracellular
removies = data_ne.query_acc[data_ne.evalue==0].unique()
data = data[~data['query_acc'].isin(removies)]
data_e = data_e[~data_e['query_acc'].isin(removies)]

### removes is an array, removies2 is a list. this is a problem


# 3) remove if significant match in experiment non-extracellular
removies2 = []
# database protein and query have to be same length +- 10%
data_ne = data_ne[~(data_ne.iloc[:, 7] > data_ne.iloc[:, 2]*1.10)]
data_ne = data_ne[~(data_ne.iloc[:, 7] < data_ne.iloc[:, 2]*0.90)]
removies2 += data_ne['query_acc'].unique().tolist()
data = data[~data['query_acc'].isin(removies2)]
data_e = data_e[~data_e['query_acc'].isin(removies2)]
SOCK = pd.DataFrame({'SOCK': SOCK})

a = SOCK.isin(removies)
a_list = a.values.tolist()
a_list.count([True])

if a_list.count([True]) > 0:
    indices = [i for i, x in enumerate(a_list) if x == [True]]
    SOCK = SOCK.drop(indices)
    
a = SOCK.isin(removies2)
a_list = a.values.tolist()
a_list.count([True])

if a_list.count([True]) > 0:
   indices = [i for i, x in enumerate(a_list) if x == [True]]
   SOCK = SOCK.drop(indices)
    
    ### up to here ###
    
#### 4) include if significant match in experimental extracellular
# filter based on e-value ### 10e-20
data_e = data_e[data_e['evalue'] < 10e-20]
data_e = data_e[data_e['subject_length'] < data_e['query_length']*1.20 ]
data_e = data_e[data_e['subject_length'] > data_e['query_length']*0.80 ]
    
temp = data_e['query_acc'].unique()
if len(temp) > 0:
    temp_df = pd.DataFrame(temp, columns=SOCK.columns)
    SOCK = pd.concat([SOCK, temp_df], ignore_index=True)

    
#### 5) include if EXACT match in computational psortb
temp = data[data['evalue'] == 0].index
if len(temp) > 0:
    data = data.loc[temp]
    new = data['query_acc'].unique()
    new = new.tolist()
    SOCK = SOCK.append(new)

def combine_columns(row):
    values = [str(val) for val in row if not pd.isna(val)]
    return ' '.join(values)

# apply the function to each row of the dataframe
SOCK['combined'] = SOCK.apply(combine_columns, axis=1)   
 
blastE = SOCK['combined'].unique()
print(len(blastE)) # how many extracellular proteins are there
# save output
blastE_list = blastE.tolist()
os.chdir(script_path)
pd.DataFrame(blastE_list).to_csv(outfile, index=False)

######################################################################
################### ANTISMASH ########################################
######################################################################

# load list of antismash types
types_dir = antismashtypes
antismash_types = pd.read_csv(types_dir,encoding='latin-1')

FILEOUT = antismashoutputpath
ACC = accession
base_dir = antismashdirectory

# move to region gbk directory
os.chdir(base_dir)
# make output matrix
BIGGERLIST = pd.DataFrame(columns=["gene_kind", "product", "locus_tag", "protein_id", "region"])
# make list of region gbk files
files = [f for f in os.listdir() if f.endswith('.gbk')]
#### extract antismash type first, as its hard to do later
PROD = [None] * len(files)

for j in range(len(files)):
    two_num = str(j + 1).zfill(2)
    filefile = ACC + ".region0" + two_num + ".gbk"
    
    # read a genbank file
    with open(filefile, "r") as f:
        data = f.read().splitlines()
    data = [line.strip() for line in data]
    data = [line for line in data if line != ""]
    # split into chunks based on protocluster
    start = [i for i, line in enumerate(data) if "protocluster" in line][0]
    end = start + 100
    vec1 = data[start:end]
    bb = [i for i, line in enumerate(vec1) if "product" in line]
    prod = re.sub('\\.*product=','',vec1[bb[0]])
    prod = re.sub(' ', '', prod)
    prod = re.sub('/', '', prod)
    ##  ## problem line
    if prod is not None:
        PROD[j] = prod
# extract other information
for j in range(len(files)):
    two_num = str(j + 1).zfill(2)
    filefile = ACC + ".region0" + two_num + ".gbk"
    
    # read a genbank file
    with open(filefile, "r") as f:
        data = f.read().splitlines()
    data = [line.strip() for line in data]
    data = [line for line in data if line != ""]
    start = [i for i, s in enumerate(data) if 'CDS' in s]
    check = [i for i, s in enumerate(data) if re.search("CDS_", s)]
    if len(check) > 0:
        start = [x for x in start if x not in check]
    end = [i-1 for i in start]
    end = end[1:]
    end.append(start[-1]+50)
    seq = range(0, len(start))
    var = list(seq)
    for k,l in enumerate(var, start=1):
        exec(f"vec{l} = k")
        exec(f"vec{l} = data[start[{l}]:end[{l}]]")
##

    BIGLIST = np.chararray([len(start), 5])
    BIGLIST = pd.DataFrame(BIGLIST, columns=["gene_kind", "product", "locus_tag", "protein_id", "region"])


    my_dict = {}
    max_lists = len(start)
    for i in range(max_lists):
        key = f'vec{i}'
        if key in locals():
            my_dict[i] = locals()[key]
        else:
            my_dict[i] = None

    for i in range(len(start)):
        VEC = my_dict[i]
        aa = [i for i, s in enumerate(VEC) if 'gene_kind' in s]
        bb = [i for i, s in enumerate(VEC) if 'product' in s]
        cc = [i for i, s in enumerate(VEC) if 'locus_tag' in s]
        dd = [i for i, s in enumerate(VEC) if 'protein_id' in s]
        vec = VEC
        # gene_kind
        gk = ""
        if len(aa) > 0:
            gk = vec[aa[0]].replace('.*gene_kind=','').replace(' ', '').replace('/', '')
            # product
        if len(bb) > 0:
            prod = re.sub(r'\.*product=', '', vec[bb[0]])
            prod = re.sub(r'\s+', '', prod)
            prod = re.sub('/', '', prod)
        else:
            prod = '' # or whatever default value you want
        # locus tag
        if len(bb) > 0:
            lt = re.sub('\\.*locus_tag=', '', vec[cc[0]], flags=re.IGNORECASE).strip().replace(' ', '').replace('/', '')
        #protein_id
        if len(dd) > 0:
            pid = re.sub('\\.*protein_id=', '', vec[dd[0]], flags=re.IGNORECASE).strip().replace(' ', '').replace('/', '')    
        if len(gk) > 0:
            BIGLIST.iloc[i,0] = gk
        if len(prod) > 0:
            BIGLIST.iloc[i,1] = PROD[j]
        if len(lt) > 0:
            BIGLIST.iloc[i,2] = lt
        if len(pid) > 0:
            BIGLIST.iloc[i,3] = pid
        BIGLIST.iloc[:,4] = j
    BIGGERLIST = pd.concat([BIGGERLIST, BIGLIST])

#### dataframe has random extra columns that need deleting
    
# Set working directory
os.chdir(script_path)

# Create DataFrame from BIGGERLIST
df = pd.DataFrame(BIGGERLIST)

# Write DataFrame to CSV file
df.to_csv(FILEOUT, index=False)   
    
# Read the CSV file
data = pd.read_csv(FILEOUT)
    
# Filter out genes of unknown type
data = data[data['gene_kind'].str.contains('gene_kind')]
        
# Filter out non-social secondary metabolite types
social_types = antismash_types.loc[antismash_types['Social'] == 1, 'Label'].tolist()
data['product'] = data['product'].str.replace('"', '')
data['locus_tag'] = data['locus_tag'].str.replace('"', '')
data['protein_id'] = data['protein_id'].str.replace('"', '')
data = data[data['product'].isin(social_types)]
data.to_csv(FILEOUT[:-4] + "_filtered.csv", index=False)
os.remove(FILEOUT)

######################################################################
################### COMBINE ##########################################
######################################################################

## code that loads the three output csv, and combines into one
os.chdir(script_path)
FILEOUT = "antismash.csv"
outputpath = "kofam.csv"
blastoutputpath = "blast.csv"
data_a = pd.read_csv(FILEOUT[:-4] + "_filtered.csv")
data_a = pd.DataFrame(data_a['protein_id'])
data_k = pd.read_csv(outputpath)
data_b = pd.read_csv(blastoutputpath)
data_k.rename(columns={'0': 'protein_id'}, inplace=True)
data_b.rename(columns={'0': 'protein_id'}, inplace=True)
# Concatenate the dataframes
combined_data = pd.concat([data_a, data_b, data_k])
# Extract the unique entries
unique_entries = combined_data.iloc[:, 0].drop_duplicates()
unique_entries.to_csv(SOCKSS, index=False)







