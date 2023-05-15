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
parser.add_argument('-gff', '--gffinput', metavar='.', type=str, required=True, help='Path to directory of region_gbk files')
parser.add_argument('-O', '--KOFAMoutput', type=str, metavar='.', required=True, help='KOFAM output file e.g. kofam.txt')
parser.add_argument('-p', '--gramPositive', action='store_true', help='Gram stain (positive)')
parser.add_argument('-n', '--gramNegative', action='store_true', help='Gram stain (negative)')


# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the command-line arguments
genome = args.GENOMEinput
fasta = args.FASTAinput
gff = args.gffinput
OUTPUT = args.KOFAMoutput


### include an option for building blast databases
## a flag called m that determines if the databases must be made
#fasta_folder = "/Users/laurie/Dropbox/POSTDOC/PROJECT_3_SOCIAL_GENES/psort_simple"
#database_folder = "blast_databases"

#bash_script = """#!/bin/bash
#mkdir blast_databases
#makeblastdb -in {ff}/psort_extracellular_gramP.fasta -parse_seqids -dbtype prot -out {db}/blastdbP &
#makeblastdb -in {ff}/psort_extracellular_gramN.fasta -parse_seqids -dbtype prot -out {db}/blastdbN &
#makeblastdb -in {ff}/EXTRA.fasta -parse_seqids -dbtype prot -out {db}/blastdbCExtra &
#makeblastdb -in {ff}/NON_EXTRA.fasta -parse_seqids -dbtype prot -out {db}/blastdbCNonExtra
#"""
# Format the bash script with the input and output filenames
#formatted_script1 = bash_script.format(db=database_folder, ff=fasta_folder)
# Create a temporary file to store the bash script
#script_filename1 = 'makedb.sh'
#with open(script_filename1, 'w') as script_file:
 #   script_file.write(formatted_script1)
# Execute the bash script
# only run this line if you need to build the databases
###subprocess.run(['bash', script_filename1], check=True)


######################################################################
################### Make bash script ##################################
######################################################################

## user will need to add to their path
## e.g. export PATH="/Users/laurie/Dropbox/POSTDOC/PROJECT 3_SOCIAL GENES/psort_simple/kofam_scan-1.3.0:$PATH"
## user will need to modify config.yml to include correct paths to prokaryotes.hal and ko.list

## general
input_genome = genome
## kofam
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

### make a bash script
# Define the bash script contents with variables
bash_script = """#!/bin/bash
mkdir {output}
export PATH="/Users/laurie/Dropbox/POSTDOC/PROJECT_3_SOCIAL_GENES/psort_simple/kofam_scan-1.3.0:$PATH"
exec_annotation {input} -o {outputK} -f detail-tsv --threshold-scale 0.75 &
### blast to gram 
blastp -db {db}/{gram} -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {output}/file_PSORT.txt &
##### blast to proven extracellular
blastp -db {db}/blastdbCExtra -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {output}/file_PSORT_E.txt &
###### blast to proven non-extracellular
blastp -db {db}/blastdbCNonExtra -query {input} -evalue 10e-8 -outfmt \
"6 sseqid qacc qlen evalue bitscore sstart send slen" -out {output}/file_PSORT_NE.txt
antismash {fna} --genefinding-gff3 {gff}
"""

# Format the bash script with the input and output filenames
formatted_script = bash_script.format(input=input_genome, outputK=output_name, output=output_folder, db=database_folder, gram=gramDB, fna=fna, gff=gff)
# Create a temporary file to store the bash script
script_filename = 'SOCfinder.sh'
with open(script_filename, 'w') as script_file:
    script_file.write(formatted_script)

# Execute the bash script
subprocess.run(['bash', script_filename], check=True)

# Remove the temporary script file
#subprocess.run(['rm', script_filename])






