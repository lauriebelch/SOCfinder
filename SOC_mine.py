#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 10:40:00 2023

@author: laurie
"""

import subprocess
import argparse
import os
import sys
import gffutils
import pandas as pd

# store directory that script is running in
script_path = os.path.dirname(os.path.abspath(sys.argv[0]))

# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='SOCfinder parser')

# Add command-line arguments with flags and help messages
parser._optionals.title = 'Required arguments'
parser.add_argument('-g', '--GENOMEinput', type=str, metavar='.', required=True, help='Path to GENOME protein (.faa)')
parser.add_argument('-f', '--FASTAinput', type=str, metavar='.', required=True, help='Path to GENOME nucleotide (.fna)')
parser.add_argument('-gff', '--gffinput', metavar='.', type=str, required=True, help='Path to GENOME gff file')
parser.add_argument('-O', '--outputfolder', type=str, metavar='.', required=True, help='Name of output folder')

# Mutually exclusive group for Gram staining options with required=True
gram_group = parser.add_mutually_exclusive_group(required=True)
gram_group.add_argument('-p', '--gramPositive', action='store_true', help='Gram stain (positive)')
gram_group.add_argument('-n', '--gramNegative', action='store_true', help='Gram stain (negative)')
gram_group.add_argument('-both', '--gramBoth', action='store_true', help='Search both Gram-positive and Gram-negative')

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the command-line arguments
genome = os.path.abspath(args.GENOMEinput)
fasta = os.path.abspath(args.FASTAinput)
gff = os.path.abspath(args.gffinput)
OUTPUT = os.path.abspath(args.outputfolder)


######################################################################
################### Make bash script ##################################
######################################################################

## user will need to add to their path
## e.g. export PATH="/Users/laurie/Dropbox/POSTDOC/PROJECT 3_SOCIAL GENES/psort_simple/kofam_scan-1.3.0:$PATH"
## user will need to modify config.yml to include correct paths to prokaryotes.hal and ko.list

## general
input_genome = genome
output_name = OUTPUT
## blaste
output_folder = "blast_outputs"
database_folder = os.path.join(script_path, "blast_databases")
## anti
fna = fasta
gff = gff
base_name = os.path.splitext(os.path.basename(gff))[0]
gff1 = os.path.join(os.path.dirname(gff), f"{base_name}_1.gff")

print("FINDING SOCKS ....")

## define correct database depending on gram + or -
# Define the initial value of gramDB
gramDB = ""
# Use conditional checks to assign the appropriate database based on the Gram stain selection
if args.gramPositive:
    gramDB = "blastdbP"
elif args.gramNegative:
    gramDB = "blastdbN"
elif args.gramBoth:
    gramDB = "blastdbBoth"
else:
    # This else block is now redundant as argparse ensures one option is always selected
    print("WARNING - NO GRAM STAIN SELECTED - WARNING")


directory_name = output_name
if not os.path.exists(directory_name):
    os.makedirs(directory_name)
else:
    print("WARNING - The directory already exists! - WARNINGS")

blast_outputs_dir = os.path.join(OUTPUT, "blast_outputs")
kofam_outputs_dir = os.path.join(OUTPUT, "kofam.txt")
adir = os.path.join(OUTPUT, "anti_smash")
tempdir = os.path.join(OUTPUT, "tmp")


### fix gff
database_filename = "OG_GFF"

# Create a GFF database to parse the GFF file
db = gffutils.create_db(gff, database_filename,force=True, keep_order=True, merge_strategy='create_unique')

# Extract header lines and feature information
header_lines = []
data = []
for line in open(gff):
    if line.startswith("##FASTA"):
        break
    if line.startswith("#"):
        header_lines.append(line.strip())
    else:
        fields = line.strip().split("\t")
        seqid, source, featuretype, start, end, score, strand, frame, attributes = fields
        attributes = dict(item.split("=") for item in attributes.split(";"))
        feature_data = {
            "Sequence ID": seqid,
            "Source": source,
            "Feature Type": featuretype,
            "Start": int(start),
            "End": int(end),
            "Score": score,
            "Strand": strand,
            "Frame": frame,
            "Attributes": attributes
        }
        data.append(feature_data)

# Create a DataFrame from the list of feature data
df = pd.DataFrame(data)

# Function to convert attributes from new format to original format
def convert_attributes(attributes_dict):
    attributes_list = [f"{key}={value}" for key, value in attributes_dict.items()]
    return ";".join(attributes_list)

# Convert the 'Attributes' column from dictionaries to strings in the original format
df['Attributes'] = df['Attributes'].apply(convert_attributes)


# Generate the output file name
base_name = os.path.splitext(os.path.basename(gff))[0]
output_file = os.path.join(os.path.dirname(gff), f"{base_name}_1.gff")

# Dictionary to store 'region' ends for each Sequence ID
region_ends = {}

# List to store indices of rows to be removed
rows_to_remove = []

# Loop through each row in the DataFrame
for idx, row in df.iterrows():
    sequence_id = row["Sequence ID"]
    feature_type = row["Feature Type"]
    end = row["End"]

    # Check if the feature is a 'region'
    if feature_type == "region":
        # Store 'region' end for this Sequence ID
        region_ends[sequence_id] = end
    else:
        # Check if the feature's 'end' is greater than the 'region' end
        if sequence_id in region_ends and end > region_ends[sequence_id]:
            rows_to_remove.append(idx)

####### new
# List to store indices of rows to be removed based on the new condition
rows_to_remove2 = []

for idx, row in df.iterrows():
    start = row["Start"]
    end = row["End"]
    feature_type = row["Feature Type"]

     # Calculate the length of the feature
    feature_length = end - start + 1

    # Check if the length of the feature type is greater than 15
    if len(feature_type) > 15:
        rows_to_remove.append(idx)

# Check for overlapping genes with the same end position and keep the longest
cds_entries = df[df["Feature Type"] == "CDS"].copy()
cds_entries.sort_values(by=["Sequence ID", "End", "Start"], inplace=True)
to_remove_overlap = []

for seqid, group in cds_entries.groupby("Sequence ID"):
    end_positions = group.groupby("End")
    for end, sub_group in end_positions:
        if len(sub_group) > 1:
            longest_idx = sub_group.loc[(sub_group["End"] - sub_group["Start"]).idxmax()].name
            to_remove_overlap.extend(sub_group.index.difference([longest_idx]))

rows_to_remove.extend(to_remove_overlap)

# Combine the rows to be removed from both conditions
rows_to_remove_all = rows_to_remove + rows_to_remove2

# Save the removed rows as a text file named "gff_removed.txt"
removed_records = df.iloc[rows_to_remove_all]

output_removed_file = os.path.join(os.path.dirname(gff), "gff_removed.txt")
removed_records.to_csv(output_removed_file, sep="\t", index=False)


# Remove the rows from the DataFrame
df = df.drop(rows_to_remove_all)

# Print the message about the number of removed records
num_removed_records = len(rows_to_remove_all)
if num_removed_records == 0:
    print("No GFF records removed")
else:
    print(f"{num_removed_records} GFF records removed")


# Save the DataFrame as the new GFF file, appending the feature lines
with open(output_file, "w") as out_file:
    for header_line in header_lines:
        out_file.write(f"{header_line}\n")
    df.to_csv(out_file, sep="\t", header=False, index=False, mode="a")

os.remove(database_filename)
### end fix GFF

kofam_profiles = os.path.join(script_path, 'KOFAM', 'profiles', 'prokaryote.hal')
ko_list = os.path.join(script_path, 'KOFAM', 'ko_list')

### make a bash script
# Define the bash script contents with variables
bash_script = """#!/bin/bash
### if you need to add exec_annotation to path, do it here ##
mkdir {blast_outputs_dir}
exec_annotation {input} -o {outputK} -f detail-tsv --threshold-scale 0.75 --tmp-dir={DIR} --cpu=32 \
-p {kofam_profiles} -k {ko_list} &
### blast to gram
blastp -db {db}/{gram} -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {blast_outputs_dir}/file_PSORT.txt -num_threads 32 &
##### blast to proven extracellular
blastp -db {db}/blastdbCExtra -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {blast_outputs_dir}/file_PSORT_E.txt -num_threads 32 &
###### blast to proven non-extracellular
blastp -db {db}/blastdbCNonExtra -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {blast_outputs_dir}/file_PSORT_NE.txt -num_threads 32 &
antismash {fna} --genefinding-gff3 {gff} --output-dir {adir}
### save the accession number
file={fna}
first_line=$(head -n 1 "$file")
acc=$(echo "$first_line" | cut -d' ' -f1)
acc1=$(echo "$acc" | cut -c2- )
echo $acc1 > {adir}/accession.txt
"""

# Format the bash script with the input and output filenames
formatted_script = bash_script.format(input=input_genome, outputK=kofam_outputs_dir, blast_outputs_dir=blast_outputs_dir, db=database_folder, gram=gramDB, fna=fna, gff=gff1, adir=adir, DIR=tempdir,kofam_profiles=kofam_profiles,
    ko_list=ko_list
)
# Create a temporary file to store the bash script
print("Running SOCfinder modules......")
print("Why not grab a coffee?")
script_filename = os.path.join(directory_name, 'SOCfinder.sh')
with open(script_filename, 'w') as script_file:
    script_file.write(formatted_script)

# Execute the bash script
#subprocess.run(['bash', script_filename], check=True)

process = subprocess.Popen(['bash', script_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = process.communicate()

# Wait for the process to finish
process.wait()


# Remove the temporary script file
subprocess.run(['rm', script_filename])

# remove the temp file
rdir = os.path.join(tempdir, "tabular/tabular.txt")
subprocess.run(['rm', rdir])
subprocess.run(['rm', '-rf', tempdir])

print("SOCfinder complete. Check the files in", directory_name)
