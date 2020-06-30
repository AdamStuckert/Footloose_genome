Downloaded data.
Checked md5sums.
Navigated to the folders with data. Converted good reads to fasta files.

```bash
samtools fasta m54312U_200614_135548.subreads.bam > $HOME/footloose_genome/raw_data/S_parvus_smrtcell_3.fasta

samtools fasta m64019_200622_194744.subreads.bam > $HOME/footloose_genome/raw_data/S_parvus_smrtcell_2.fasta
```

Quick diagnostics on read data:

```bash
cd $HOME/footloose_genome/raw_data

# use my python script
readlengths.py S_parvus_smrtcell_2.fasta S_parvus_smrtcell_2.fasta.txt
readlengths.py S_parvus_smrtcell_3.fasta S_parvus_smrtcell_3.fasta.txt
```

Data. 

SMRTcell | Reads | Bases | Average read length | Read N50
-- | -- | -- | -- | --
1 | 10,075,593| 104,320,373,130 | 10,353.8 | 13,626
2 | 6,467,364 | 61,060,153,324 | 9,441.3 | 17,321
3 | 6,595,111 | 62,098,555,689 | 9,415.8 | 17,826
Total | 23,138,068‬ | 227,479,082,143 | |
