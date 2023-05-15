
![SOCfinder](Soc_finder_v4.png)
SOCfinder is a tool for finding cooperative genes in bacterial genomes.

## Installation

There are several prerequisities:
--list them --

## Download SOCfinder scripts
```bash
some code
```

## make BLAST databases
-- provide a python script to make the blast databases

## Tutorial

-- explain the files it needs, and where they are

```python
./SOC_mine.py -g test2/P_salmonis.faa -f test2/P_salmonis.fna -gff test2/P_salmonis.gff3 -O P_salmonis_kofam.txt -n

./SOC_parse.py -ki P_salmonis_kofam.txt -b blast_outputs/ -ad P_salmonis/ -ac NZ_CP013821.1 -k inputs/SOCIAL_KO.csv -a inputs/antismash_types.csv -so socks.csv -ko kofam.csv -ao antismash.csv -bo blast.csv
```

## Manuscript

-- link to preprint

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License

Free to use
