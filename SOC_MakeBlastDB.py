#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:02:16 2023

@author: laurie
"""

import subprocess

fasta_folder = "blast_files"
database_folder = "blast_databases"

bash_script = """#!/bin/bash
mkdir blast_databases
makeblastdb -in {ff}/psort_extracellular_gramP.fasta -parse_seqids -dbtype prot -out {db}/blastdbP
makeblastdb -in {ff}/psort_extracellular_gramN.fasta -parse_seqids -dbtype prot -out {db}/blastdbN
makeblastdb -in {ff}/EXTRA.fasta -parse_seqids -dbtype prot -out {db}/blastdbCExtra
makeblastdb -in {ff}/NON_EXTRA.fasta -parse_seqids -dbtype prot -out {db}/blastdbCNonExtra
"""
# Format the bash script with the input and output filenames
formatted_script1 = bash_script.format(db=database_folder, ff=fasta_folder)
# Create a temporary file to store the bash script
script_filename1 = 'makedb.sh'
with open(script_filename1, 'w') as script_file:
    script_file.write(formatted_script1)
# Execute the bash script
subprocess.run(['bash', script_filename1], check=True)
