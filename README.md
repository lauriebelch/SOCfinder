**SOCFinder**
![SOCfinder](Soc_finder_v4.png)

SOCfinder is a bioinformatics tool for finding cooperative genes in bacterial genomes. SOCfinder combines information from several methods, considering if a gene is likely to: (1) code for an extracellular protein; (2) have a cooperative functional annotation; or (3) be part of the biosynthesis of a cooperative secondary metabolite. SOCfinder uses information on the quality and significance of database matches and annotations.

## Installation

The easiest way to install is to clone this github page, or download the zip file https://github.com/lauriebelch/SOCfinder/archive/refs/heads/main.zip

You can then use the environment.yml file to create a conda environment with most of the required packages and tools

## Download SOCfinder scripts

```bash
git clone https://github.com/lauriebelch/SOCfinder.git
cd SOCfinder
conda env create -f environment.yml
```

## Download KOFAMscan files

You will need to download some required files for KOFAMscan. It is reccommended that you do this within the SOCfinder folder

```bash
mkdir KOFAM
cd ./KOFAM
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
wget ftp://ftp.genome.jp/pub/tools/kofam_scan/kofam_scan-1.3.0.tar.gz
wget https://www.genome.jp/ftp/tools/kofam_scan/README.md
tar -xf profiles.tar.gz
tar -xf kofam_scan-1.3.0.tar.gz
gunzip ko_list.gz
```
You will also need to edit the config file

```bash
pwd 
```
For me its /drives/Laurie/KOFAM
```bash
cd ./kofam_scan-1.3.0/
nano config-template.yml
```
on line 4, change it to
```bash
profile: /drives/Laurie/KOFAM/profiles
```
on line 7, change it to 
```bash
ko_list: /drives/4tb/Laurie/KOFAM/ko_list
```
on line 18, change it to
```bash
cpu: 32
```
ctrl-O and save as ‘config.yml’

You will also need to add the exec_annotation file to your path

```bash
nano ~/.bash_profile
export PATH="/path/to/kofam_scan-1.3.0:$PATH"
```
ctrl-O to save

## Install antismash

You will need to download some required files for ANTISMASH. It is reccommended that you do this within the SOCfinder folder

```bash
mkdir ANTISMASH
cd ./ANTISMASH 
wget https://dl.secondarymetabolites.org/releases/6.1.1/antismash-6.1.1.tar.gz
tar -zxf antismash-6.1.1.tar.gz
pip install ./antismash-6.1.1
python antismash-6.1.1/antismash/download_databases.py
```

## make BLAST databases

You will need to build the databases that the BLAST search uses
```bash
./SOC_MakeBlastDB.py
```

## Tutorial

Here are some genomes blah blah

-- Part 1: Mine the Genome
```python
./SOC_parse.py -ac CP013821.1 -i P_salmonis/ -k inputs/SOCIAL_KO.csv -a inputs/antismash_types.csv
./SOC_parse.py -ac CP002703.1 -i B_aphidicola/ -k inputs/SOCIAL_KO.csv -a inputs/antismash_types.csv
```

-- Part 2: Extract the Social Genes
```python
./SOC_mine.py -g test/B_aphidicola.faa -f test/B_aphidicola.fna -gff test/B_aphidicola.gff -O B_aphidicola -n
./SOC_mine.py -g test/P_salmonis.faa -f test/P_salmonis.fna -gff test/P_salmonis.gff -O P_salmonis -n 
```

## User Options

**SOC_mine.py**

- `-g GENOMEinput`
  - Path to GENOME protein (.faa)
- `-f FASTAinput`
  - Path to GENOME nucleotide (.fna)
- `-gff GENOMEinput`
  - Path to GENOME gff file (.gff)
- `-O outputfolder`
  - Name of output folder
- `-p -n GramPositive | GramNegative`
  - Gram stain (positive | negative)

**SOC_parse.py**
- `-i inputfolder`
  - Path to input folder from SOC_mine
- `-ac accession`
  - Accession number
- `-k ko`
  - Path to list of social KO terms
- `-a ANTISMASHtypes`
  - Path to list of antismash types

## How to download genomes

The SOCfinder reccommended way to download the genome files you need is to use the [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=datasets-command-line-20221012) command line tool. This is so that gene ID is the same in protein fasta, nucleotide fasta, and gff. Otherwise users will have to check the gene ID .

```python
datasets download genome accession GCA_003798305.1 --include gff3,genome,protein
```

## Manuscript

-- link to preprint

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

Free to use, but please cite the SOCfinder paper
