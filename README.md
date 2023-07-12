**SOCFinder**
![SOCfinder](Soc_finder_v4.png)

SOCfinder is a bioinformatics tool for finding cooperative genes in bacterial genomes. SOCfinder combines information from several methods, considering if a gene is likely to: (1) code for an extracellular protein; (2) have a cooperative functional annotation; or (3) be part of the biosynthesis of a cooperative secondary metabolite. SOCfinder uses information on the quality and significance of database matches and annotations.

## Table of contents
- [Installation](#Installation)
- [How to Use](#Tutorial)
- [Input Options](#Options)

## Installation

You will need miniconda, which can be installed by following the instructions [here](https://docs.conda.io/en/latest/miniconda.html). You can check that conda has installed correctly by running `conda list` (you may need to restart your terminal first).

For an introduction to conda, see [here](https://www.machinelearningplus.com/deployment/conda-create-environment-and-everything-you-need-to-know-to-manage-conda-virtual-environment/).

If you are on a mac that has an M1 or M2 chip, you might have to adjust your conda architecture. Instructions can be found [here](#Mac).

## Download SOCfinder

You can download SOCfinder from github using the code below.

```bash
git clone https://github.com/lauriebelch/SOCfinder.git
cd SOCfinder
conda env create -f environment_noversion.yml
# activate conda environment
conda activate SOCfinder
```

You will then need to download some files for KOFAMscan and ANTISMASH. The easiest way to do this is to use the helper script.

```bash
chmod +x ./helper_script
source helper_script
```

When this script has finished running, it will tell you how to add the required programs to your path. For a simple explanation of the path, see [here](https://janelbrandon.medium.com/understanding-the-path-variable-6eae0936e976).

## Make BLAST databases

You will need to build the databases that the BLAST search uses. You only need to do this once, and can use the script provided.

```bash
cd blast_files
unzip Archive.zip
cd ..
chmod +x ./SOC_MakeBlastDB.py
./SOC_MakeBlastDB.py
```

## Tutorial

SOCfinder comes with two genomes that you can test the code with. For each genome, you need a protein fasta, a nucleotide fasta, and a gff. This is because the tools that SOCfinder uses require different inputs.

The folder `test` contains the files for a strain of *Buchnera aphidicola*.
The folder `test2` contains the files for a strain of *Piscirickettsia salmonis*.

**Part 1: Mine the Genome.**
In this section, the three modules of SOCfinder are run. The output files are stored in a folder
```python
./SOC_mine.py -g test/B_aphidicola.faa -f test/B_aphidicola.fna -gff test/B_aphidicola.gff -O B_aphidicola -n
./SOC_mine.py -g test2/P_salmonis.faa -f test2/P_salmonis.fna -gff test2/P_salmonis.gff -O P_salmonis -n 
```

**Part 2: Extract the Social Genes.**

In this section, the outputs of each modules are converted into lists of social genes.
```python
./SOC_parse.py -i B_aphidicola/ -k inputs/SOCIAL_KO.csv -a inputs/antismash_types.csv
./SOC_parse.py -i P_salmonis/ -k inputs/SOCIAL_KO.csv -a inputs/antismash_types.csv
```
The final list of cooperative genes is stored as `SOCKS.csv`. Outputs for each module are stored as `K_SOCK.csv` for the functional annotation social genes, `B_SOCK.csv` for the extracellular genes, and `A_SOCK_filtered.csv` for the antismash social genes. There is also a summary file `summary.csv` that gives you the counts of cooperative genes for each module.

*B. aphidicola* has nine social genes, and *P. salmonis* has 64.

## Options

Command-line options for SOCfinder

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
- `-k ko`
  - Path to list of social KO terms
- `-a ANTISMASHtypes`
  - Path to list of antismash types

## How to download genomes

The SOCfinder reccommended way to download the genome files you need is to use the [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/?utm_source=ncbi_insights&utm_medium=referral&utm_campaign=datasets-command-line-20221012) command line tool. This is so that gene ID is the same in the protein fasta, nucleotide fasta, and gff.

```python
datasets download genome accession GCA_003798305.1 --include gff3,genome,protein --filename GCA_003798305.1.zip
```

## Mac
Apple recently made the switch from Intel processors to their own Apple Silicon processors. This can cause package compatibility issues if your computer has one of the new M1 or M2 chips. Currently, the best solution is to create conda environments that still use the old architecture. You can do this by running the following command before creating the SOCfinder conda environment.

```bash
conda config --add subdirs osx-64
```

Further discussion of this issue can be found [here](https://towardsdatascience.com/how-to-manage-conda-environments-on-an-apple-silicon-m1-mac-1e29cb3bad12).

## Manuscript

The manuscript "SOCfinder: a genomic tool for identifying cooperative genes in bacteria" has been recently submitted to *Microbial Genomics*

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

SOCfinder code is open-source and free to use and distribute, but please cite the SOCfinder paper.
antiSMASH is an open source tool available under the GNU Affero General Public License version 3.0 or greater.
KOFAMscan is released under the MIT License.

