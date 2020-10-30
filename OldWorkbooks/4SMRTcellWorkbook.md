Downloaded data from the 4th SMRTcell and checked md5sums. Converted good reads (subreads) to fasta files.

Quick diagnostics on read data:

```bash
cd $HOME/footloose_genome/raw_data

# use my python script
readlengths.py S_parvus_smrtcell_4.fasta S_parvus_smrtcell_4.fasta.txt
```

### Overview of all data 

SMRTcell | Reads | Bases | Average read length | Read N50
-- | -- | -- | -- | --
1 | 10,075,593| 104,320,373,130 | 10,353.8 | 13,626
2 | 6,467,364 | 61,060,153,324 | 9,441.3 | 17,321
3 | 6,595,111 | 62,098,555,689 | 9,415.8 | 17,826
4 | 8,972,917 | 84,993,649,263 | 9,472.2 |  14,344
Total | ‬32,110,985 | 312,472,731,406 | 9,731.0 | 15,262

Ran an assembly using wtdbg2 on the combined datasets. 


### Genome assembly comparisons

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.1.0 | 2,143,282,740 | 28,208 | 10,174 | 0.00 | C:0.6%[S:0.6%,D:0.0%],F:1.0%,M:98.4%,n:3950
S_parvus.2.0 | 3,782,941,384 | 156,120 | 49,172 | 0.00 | C:28.7%[S:28.6%,D:0.1%],F:14.2%,M:57.1%,n:3950
S_parvus.2.0_lowcoverage  | 3,756,961,657 | 197,812 | 42,434 | 0.00 | C:39.8%[S:39.3%,D:0.5%],F:15.2%,M:45.0%,n:3950
S_parvus.3.0 | 3,947,455,185 | 401,758 | 27,869 | 0.0 | C:56.5%[S:55.7%,D:0.8%],F:14.6%,M:28.9%,n:3950

### Calculated depth:

Total number of bases from all 4 SMRT cells: ????

Average coverage of assembled genome: ?????

## Concluding thoughts:

