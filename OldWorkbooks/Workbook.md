### Genome assembly

We have gotten PacBio reads from Mt. Sinai. The data are from 3 SMRTcells of the PacBio Sequel II and are placed into a directory called `raw_PacBio_data/`. Data come in a `.bam.` format, and only the `subreads.bam` are good for genome assembly. The first task is to convert these `.bam` files into `.fasta` files for assembly. I have a small script setup to run a conversion from `.bam` to `.fasta.`. This is entitled `pbsamtools.job`:

```#!/bin/bash
#SBATCH --job-name=pbsamtools
#SBATCH --output=pbsamtools.log
#SBATCH --cpus-per-task=24
#SBATCH --partition=macmanes,shared
# echo commands to stdout
set -x

DIR=$(pwd)

# load environments
module purge
module load linuxbrew/colsa

# input variables
INPUT=$1

# script
samtools fasta ${INPUT}.bam > ${INPUT}.fasta
```

And command line to submit a job for each `subreads.bam` from the PacBio Sequel II:

```bash
### create fasta files
cd ${DIR}/raw_PacBio_data/
# list subread bam files
PBdat=$(ls */*subreads.bam | sed "s/.bam//g")

# quick for loop
for data in $PBdat
do
BAMDIR=$(dirname $data)
BAM=$(basename $data)
# change directory
cd $BAMDIR
# submit job
sbatch pbsamtools.job $BAM
cd ..
done
```

Next, combine all the fasta files into a single concatenated file:

```bash
cat <(ls */*subreads.fasta) > concatenated_PB_subreads.fasta
```

Finally, run an inital assembly with the assembler `wtdbg2` using this submission script.

```#!/bin/bash

#SBATCH --partition=macmanes
#SBATCH -J footloosegenome
#SBATCH --output wtdbg2.log
#SBATCH --cpus-per-task=40
#SBATCH --mem 700Gb
#SBATCH --exclude=node117,node118

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
-i $HOME/footloose_genome/raw_PacBio_data/concatenated_PB_subreads.fasta

# run consensus
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns \
-t 40 \
-i S_parvus_wtdbg.ctg.lay.gz \
-fo S_parvus_wtdbg.ctg.fa
```
