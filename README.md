
![SOCfinder](Soc_finder_v4.png)
SOCfinder is a tool for finding cooperative genes in bacterial genomes.

## Installation

There are several prerequisities:
--list them --
- python version3.10

## Instructions
- adding KOFAM to path
- altering KOFAM config file

## Download SOCfinder scripts
```bash
https://github.com/lauriebelch/SOCfinder/archive/refs/heads/main.zip
```

## make BLAST databases
```bash
./SOC_MakeBlastDB.py
```

## Tutorial

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

## Options

**SOC_mine.py**

- `-g GENOMEinput`
  - Path to GENOME protein (.faa)
- `-f FASTAinput`
  - Path to GENOME nucleotide (.fna)
- `-gff GENOMEinput`
  - Path to GENOME gff file (.gff)
- `-O outputfolder`
  - Name of output folder
- `-p | -n GramPositive | GramNegative`
  - Gram stain (positive | negative)

**SOC_parse.py
- `-i inputfolder`
  - Path to input folder from SOC_mine
- `-ac accession`
  - Accession number
- `-k ko`
  - Path to list of social KO terms
- `-a ANTISMASHtypes`
  - Path to list of antismash types

## Download genomes

-- Reccommended way is to use the download datasets package. This is so that gene ID is the same in protein fasta, nucleotide fasta, and gff. Otherwise users will have to check the gene ID .

## Manuscript

-- link to preprint

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

Free to use
