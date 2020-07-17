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

### Overview of data 

SMRTcell | Reads | Bases | Average read length | Read N50
-- | -- | -- | -- | --
1 | 10,075,593| 104,320,373,130 | 10,353.8 | 13,626
2 | 6,467,364 | 61,060,153,324 | 9,441.3 | 17,321
3 | 6,595,111 | 62,098,555,689 | 9,415.8 | 17,826
Total | 23,138,068‬ | 227,479,082,143 | |

Ran an assembly using wtdbg2 on the combined datasets. 


### Comparisons of first iteration of the genome to the second iteration:

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.1.0 | 2,143,282,740 | 28,208 | 105174 | 0.00 | C:0.6%[S:0.6%,D:0.0%],F:1.0%,M:98.4%,n:3950
S_parvus.2.0 | 3,782,941,384 | 156,120 | 49172 | 0.00 | C:28.7%[S:28.6%,D:0.1%],F:14.2%,M:57.1%,n:3950

### Calculated depth:

Total number of bases from all 3 SMRT cells: 117,978,338,550

Average coverage of assembled genome: 31.18

## Concluding thoughts:

Overall this is much better! I am a bit concerned that 50% of our data (in number of bases) is not represented in the depth files. This is indicative of missing a fairly substantial portion of the genome I think.

## Calculating genome size

Apparently those weren't my concluding thoughts. I'm trying to get a handle on the actual genome size so we can do the next steps. But traditional metrics of computationally doing this are impossible with noisy long reads. So I'm going to correct them with the canu assembler and reasses.

```bash
canu \
-p parvus -d canu \
genomeSize=6g \
correctedErrorRate=0.055 \
gnuplot=undef \
purgeOverlaps=aggressive \
corMaxEvidenceErate=0.15 \
corMhapFilterThreshold=0.0000000002 \
corMhapOptions="--threshold 0.80 --num-hashes 512 --num-min-matches 3 --ordered-sketch-size 1000 --ordered-kmer-size 14 --min-olap-length 2000 --repeat-idf-scale 50" \
mhapMemory=60g \
mhapBlockSize=500 \
ovlMerDistinct=0.975 \
gridEngineArrayOption="-a ARRAY_JOBS%30" \
-pacbio raw_data/S_parvus_smrtcell_1.fasta \
-pacbio raw_data/S_parvus_smrtcell_2.fasta \
-pacbio raw_data/S_parvus_smrtcell_3.fasta
```

Most of the parameters are from the FAQ in order to decrease overall use of the hard drive, because my lab group is already annoyed that I have blown up our disk space quota a few times...

Still no dice, blew up disk usage. On to using the PacBio `ccs` program to calculate corrected reads.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J ccs
#SBATCH --output ccs.%A_%a.log
#SBATCH --array=1-21%7   


echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

# environments
module purge
source activate ccs

movie="m54312U_200614_135548"

# run ccs in iterative chunks
ccs ${movie}.subreads.bam ${movie}.ccs.${SLURM_ARRAY_TASK_ID}.bam --chunk ${SLURM_ARRAY_TASK_ID}/21 -j 10

# merge files
pbmerge -o ${movie}.ccs.bam ${movie}.ccs.*.bam
pbindex ${movie}.ccs.bam
```
