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
parser.add_argument('-p', '--gramPositive', action='store_true', help='Gram stain (positive)')
parser.add_argument('-n', '--gramNegative', action='store_true', help='Gram stain (negative)')


# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the command-line arguments
genome = args.GENOMEinput
fasta = args.FASTAinput
gff = args.gffinput
OUTPUT = args.outputfolder

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
database_folder = "blast_databases"
## anti
fna = fasta
gff = gff

## define correct database depending on gram + or -
# Define the initial value of gramDB
gramDB = ""
# Use a conditional loop to assign the value of gramDB
if args.gramPositive:
    gramDB = "blastdbP"
elif args.gramNegative:
    gramDB = "blastdbN"
else:
    print("WARNING - NO GRAM STAIN SELECTED - WARNING")


directory_name = output_name
if not os.path.exists(directory_name):
    os.makedirs(directory_name)
else:
    print("WARNING - The directory already exists! - WARNINGS")

blast_outputs_dir = os.path.join(directory_name, "blast_outputs")
kofam_outputs_dir = os.path.join(directory_name, "kofam.txt")
adir = os.path.join(directory_name, "anti_smash")
tempdir = os.path.join(directory_name, "tmp")


### make a bash script
# Define the bash script contents with variables
bash_script = """#!/bin/bash
### if you need to add exec_annotation to path, do it here ##
mkdir {blast_outputs_dir}
exec_annotation {input} -o {outputK} -f detail-tsv --threshold-scale 0.75 --tmp-dir={DIR} --cpu=32 &

### blast to gram 
blastp -db {db}/{gram} -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {blast_outputs_dir}/file_PSORT.txt -num_threads 16 &
##### blast to proven extracellular
blastp -db {db}/blastdbCExtra -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {blast_outputs_dir}/file_PSORT_E.txt -num_threads 16 &
###### blast to proven non-extracellular
blastp -db {db}/blastdbCNonExtra -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {blast_outputs_dir}/file_PSORT_NE.txt -num_threads 16 &
antismash {fna} --genefinding-gff3 {gff} --output-dir {adir}
### save the accession number
file={fna}
first_line=$(head -n 1 "$file")
acc=$(echo "$first_line" | cut -d' ' -f1)
acc1=$(echo "$acc" | cut -c2- )
echo $acc1 > {adir}/accession.txt
"""

# Format the bash script with the input and output filenames
formatted_script = bash_script.format(input=input_genome, outputK=kofam_outputs_dir, blast_outputs_dir=blast_outputs_dir, db=database_folder, gram=gramDB, fna=fna, gff=gff, adir=adir, DIR=tempdir)
# Create a temporary file to store the bash script
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


