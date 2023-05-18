#!/bin/bash
mkdir blast_outputs
exec_annotation test2/P_salmonis.faa -o TEST -f detail-tsv --threshold-scale 0.75 &
### blast to gram 
blastp -db blast_databases/blastdbN -query test2/P_salmonis.faa -evalue 10e-8 -outfmt "6 sseqid qacc qlen evalue bitscore sstart send slen" -out blast_outputs/file_PSORT.txt &
##### blast to proven extracellular
blastp -db blast_databases/blastdbCExtra -query test2/P_salmonis.faa -evalue 10e-8 -outfmt "6 sseqid qacc qlen evalue bitscore sstart send slen" -out blast_outputs/file_PSORT_E.txt &
###### blast to proven non-extracellular
blastp -db blast_databases/blastdbCNonExtra -query test2/P_salmonis.faa -evalue 10e-8 -outfmt "6 sseqid qacc qlen evalue bitscore sstart send slen" -out blast_outputs/file_PSORT_NE.txt
antismash test2/P_salmonis.fna --genefinding-gff3 test2/P_salmonis.gff
