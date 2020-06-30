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
1 | 10075593| 104320373130 | 10353.8 | 13626
2 | 6467364 | 61060153324 | 9441.3 | 17321
3 | 6595111 | 62098555689 | 9415.8 | 17826
