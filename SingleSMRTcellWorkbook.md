### Genome assembly

We have gotten PacBio reads from Mt. Sinai. The data are from a single SMRTcell of the PacBio Sequel II (instead of the 3 we expected) and are placed into a directory called `raw_data/`. Data come in a `.bam.` format, and only the `subreads.bam` are good for genome assembly. The first task is to convert these `.bam` files into `.fasta` files for assembly. I have a small script setup to run a conversion from `.bam` to `.fasta.` and then assemble the genome with `wtdbg2`.


```bash
#!/bin/bash
#SBATCH --partition=macmanes
#SBATCH -J footloosegenome
#SBATCH --output wtdbg2.log
#SBATCH --cpus-per-task=40
#SBATCH --mem 700Gb
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

# convert subreads.bam to fasta file
samtools fasta raw_data/1_A01/m64019_200422_015111.subreads.bam > raw_data/S_parvus_smrtcell_1.fasta


module purge

mkdir S_parvus_wtdbg
cd $HOME/footloose_genome/S_parvus_wtdbg

# run assembler (genome size is a total guess here, but placed on the high end)
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 \
-o S_parvus_wtdbg \
-t 40 \
-x sq \
-g 8.5g \
-X 30 \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_1.fasta

# run consensus
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns \
-t 40 \
-i S_parvus_wtdbg.ctg.lay.gz \
-fo S_parvus_wtdbg.ctg.fa
```

After running this, I ran my genome metrics script with this command:

```bash
genomemetrics.job S_parvus_wtdbg/S_parvus_wtdbg.ctg.fa S_parvus_1.0
```

It turns out, this assembly is not good. `Contig N50 = 28kb and BUSCO (tetropoda) = 0.6% (?!)`. On a positive note, the genome size is a bit under 2.2 GB!

### Why the heck is this such a terrible assembly?! 

It could be a number of things. First, the parameters could be poor for this genome. So I tried some iterations of the assembly in which I modify the assembly parameters to fix this. I tried the parameters that were used to assemble an axolotl genome, and specified a smaller genome size:

```bash
#!/bin/bash
#SBATCH --partition=shared,macmanes
#SBATCH -J footloosegenome
#SBATCH --output wtdbg2_2.2.axolotlparameters.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

# convert subreads.bam to fasta file
# samtools fasta raw_data/1_A01/m64019_200422_015111.subreads.bam > raw_data/S_parvus_smrtcell_1.fasta


module purge

mkdir S_parvus_wtdbg_axolotlparameters
cd $HOME/footloose_genome/S_parvus_wtdbg_axolotlparameters

# run assembler (genome size is a total guess here, but placed on the high end)
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 \
-x sq \
-o S_parvus_wtdbg_axolotlparameters \
-t 40 \
-g 2.2g \
-L 5000 \
-p 21 \
-S 2 \
--aln-noskip \
--rescue-low-cov-edges \
--tidy-reads 2500 \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_1.fasta

# run consensus
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns \
-t 40 \
-i S_parvus_wtdbg_axolotlparameters.ctg.lay.gz \
-fo S_parvus_wtdbg_axolotlparameters.ctg.fa
```

This was also poor with a `contig N50 = 28kb and BUSCO = 0.5%`. Next I attempted to trim "short" reads even more aggressively, removing anything below 7.5k (which drops us down to ~40x coverage of the genome).

```bash
#!/bin/bash
#SBATCH --partition=shared,macmanes
#SBATCH -J footloosegenome
#SBATCH --output wtdbg2_2.2.axolotlparameters.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

# convert subreads.bam to fasta file
# samtools fasta raw_data/1_A01/m64019_200422_015111.subreads.bam > raw_data/S_parvus_smrtcell_1.fasta


module purge

mkdir S_parvus_wtdbg_axolotlparameters_7.5k
cd $HOME/footloose_genome/S_parvus_wtdbg_axolotlparameters_7.5k

# run assembler
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 \
-x sq \
-o S_parvus_wtdbg_axolotlparameters \
-t 40 \
-g 2.2g \
-L 5000 \
-p 21 \
-S 2 \
--aln-noskip \
--rescue-low-cov-edges \
--tidy-reads 7500 \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_1.fasta

# run consensus
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns \
-t 40 \
-i S_parvus_wtdbg_axolotlparameters.ctg.lay.gz \
-fo S_parvus_wtdbg_axolotlparameters.ctg.fa
```

**This is awful. My assembly is worse**

Next I tried....

In addition to changing the parameters, I also delved into the raw data to figure out what was going on. I wrote a python script to examine total number of reads and calculate and plot N50. Doing this, I found that the good subreads of our raw data have an N50 of  13,626. 

**insert figure here**

This is quite a bit lower than we expected given library prep aimed for 20-25 kb inserts. However, not horrible. I would expect to get an assembly that is better than roughly double the length of the raw reads!!!

Where to go from here? There are a few things that I need to examine, in no particular order.

1. Sequencing depth (align reads to ref with bwa-mem then use samtools depth **bwa-mem currently running**)
2. Is it a proliferation of recent repeats???? If this is the case, then it would be hard/impossible to assemble these repeat regions that have not had time to diverge. (Repeat Masker and/or Repeat Modeler **RepeatModeler currently running**)
3. Contamination? Could bacterial contamination and/or human contamination be driving this? Seems unlikely, but definitely need to do this sanity check.


### Sequencing depth, coverage, and mapping rate

Results from `samtools flagstat`:

```
10323702 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
5865329 + 0 supplementary
0 + 0 duplicates
9856255 + 0 mapped (95.47% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```


Ran `samtools depth`: `samtools depth SMRTcell_1.sorted.bam > depth.txt` to caluclate depth at every base.

Get total number of bases (`samtools view -H SMRTcell_1.sorted.bam  | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'` should be the same as all other estimates...): 2143282740

Calculate average coverage per base pair (`awk '{sum+=$3} END { print "Average = ",sum/2143282740}' depth.txt`): Average =  NUMBER

Next up: plot the distribution of sequencing depth.
