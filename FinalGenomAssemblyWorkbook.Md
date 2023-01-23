## Genome Assembly for the Borneon foot flagging frog

This is our approach to assembling, error correcting, and annotating this genome.

### Data overview

Downloaded data from the 5 PacBio SMRTcells and checked md5sums. Converted good reads (subreads) to fasta files using `samtools fasta`.

Quick diagnostics on read data:

```bash
cd $HOME/footloose_genome/raw_data

# use my python script
readlengths.py S_parvus_smrtcell_4.fasta S_parvus_smrtcell_5.fasta.txt
```

#### Total data:

SMRTcell | Reads | Bases | Average read length | Read N50
-- | -- | -- | -- | --
1 | 10,075,593| 104,320,373,130 | 10,353.8 | 13,626
2 | 6,467,364 | 61,060,153,324 | 9,441.3 | 17,321
3 | 6,595,111 | 62,098,555,689 | 9,415.8 | 17,826
4 | 8,972,917 | 84,993,649,263 | 9,472.2 |  14,344
4 | 7,636,119 | 71,903,661,306 | 9,416.3 |  14,381
Total | ‬39,747,104 | 384,376,392,712 | 9,670.6 | Calculate...

Ran an assembly using wtdbg2 on the combined datasets. 

```
#!/bin/bash
#SBATCH --partition=shared,macmanes
#SBATCH -J footloosegenome
#SBATCH --output Sparvus3.0.axolotlparameters.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700000

module load linuxbrew/colsa

# convert subreads.bam to fasta file
# samtools fasta raw_data/1_A01/m64019_200422_015111.subreads.bam > raw_data/S_parvus_smrtcell_1.fasta


module purge

mkdir S_parvus_wtdbg_axolotlparameters
cd $HOME/footloose_genome/S_parvus_wtdbg_axolotlparameters

# run assembler (genome size is a total guess here, but placed on the high end)
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtdbg2 \
-x sq \
-o S_parvus_wtdbg \
-t 40 \
-g 8.0g \
-L 5000 \
-p 21 \
-S 2 \
--aln-noskip \
--rescue-low-cov-edges \
--tidy-reads 2500 \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_1.fasta \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_2.fasta \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_3.fasta \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_4.fasta \
-i $HOME/footloose_genome/raw_data/S_parvus_smrtcell_5.fasta

# run consensus
/mnt/lustre/macmaneslab/macmanes/wtdbg2/wtpoa-cns \
-t 40 \
-i S_parvus_wtdbg.ctg.lay.gz \
-fo S_parvus_wtdbg.ctg.fa
```

### Genome assembly specs

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.4.0 | 3,985,797,087 | 608,715 | 23,418 | 0.0 | C:68.3%[S:67.2%,D:1.1%],F:12.6%,M:19.1%,n:3950

### Polishing the assembly with Racon 

The genic content of the assembly is not particularly good. I think this might be due to poor consensus calling. I am going to polish the assembly with the readsets that we have in order to see if this helps recover genic content that we are not seeing. It should.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=700000
set -x

DIR=$(pwd)
ASSEMBLY="$HOME/footloose_genome/S_parvus_wtdbg_axolotlparameters/S_parvus_wtdbg.ctg.fa"
genome=$(basename $ASSEMBLY)
READS=$"$HOME/footloose_genome/raw_data/S_parvus_all_smrtcells.fasta"

# racon requires a single readfile, concatenate reads:
cat raw_data/S_parvus_smrtcell_*.fasta > raw_data/S_parvus_all_smrtcells.fasta

module load linuxbrew/colsa

# preparation
mkdir racon
cd racon

#cp $ASSEMBLY .

### First align reads with minimap2
#echo aligning with minimap2
#minimap2 -I10G -t 40 -xmap-pb $genome $READS | gzip -c - > Footlose.PB.paf.gz

### Run racon
echo Polishing with racon
racon -t 240 $READS Footlose.PB.paf.gz $genome > S_parvus_wtdbg.ctg.polished1.fa
```

I ran this polishing procedure twice, which made the assembly slightly worse overall.

### Genome assembly comparisons

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.4.0 | 3,985,797,087 | 608,715 | 23,418 | 0.0 | C:68.3%[S:67.2%,D:1.1%],F:12.6%,M:19.1%,n:3950
S_parvus.4.0.polished1 | 3,982,189,551 | 611,229 | 22,402 | 0.0 | C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950
S_parvus.4.0.polished2 | 3,918,212,060 | 609,425 | 20,281 | 0.0 | C:73.3%[S:72.0%,D:1.3%],F:11.5%,M:15.2%,n:3950
S_parvus.4.0.preadspolished1 | 3,963,679,315 | 605,739 | 23,418 |  C:72.9%[S:71.6%,D:1.3%],F:11.5%,M:15.6%,n:3950

I then used the program `purge_dups` to...purge dups.

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J purge_dups
#SBATCH --output purge_dups_slurm.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118
#SBATCH --mem=300G
set -x

ASSEMBLY="$HOME/footloose_genome/genome_versions/S_parvus_wtdbg.ctg.polished1.fa"
READS=$"$HOME/footloose_genome/raw_data/S_parvus_all_smrtcells.fasta"
PAF="FootloosePB.paf.gz"

DIR=$(pwd)
genome=$(basename $ASSEMBLY)
PURGED=$(echo $genome | sed "s/.fa//" | sed "s/.fasta//")


module load linuxbrew/colsa


ln -s $ASSEMBLY

### run pipeline step by step
echo aligning with minimap2
minimap2 -I50G -t 40 -xmap-pb $genome $READS | gzip -c - > $PAF

echo initial calculations
pbcstat $PAF
calcuts PB.stat > cutoffs 2>calcults.log

echo Splitting assembly
split_fa $genome > $assembly.split
minimap2 -I50G -t 40 -xasm5 -DP $assembly.split $assembly.split | gzip -c - > $assembly.split.self.paf.gz

echo Purging haplotigs with overlaps
purge_dups -2 -T cutoffs -c PB.base.cov $assembly.split.self.paf.gz > dups.bed 2> purge_dups.log

echo Getting purged sequences...
get_seqs dups.bed $genome

# rename assembly
mv purged.fa $PURGED.purged.fa
```


Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.4.0 | 3,985,797,087 | 608,715 | 23,418 | 0.0 | C:68.3%[S:67.2%,D:1.1%],F:12.6%,M:19.1%,n:3950
S_parvus.4.0.polished1 | 3,982,189,551 | 611,229 | 22,402 | 0.0 | C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950
S_parvus.4.0.polished1.purged | 3,931,525,283 | 620,561 | 18,095 |  C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950


### Scaffolding with RNA data

What if I scaffold the assembly with bullfrog RNA seq data? Will this help recover missing genic content? I am going to run this on the `S_parvus.4.0.polished1.purged.fa` version of the genome (the one that was annotated). Data used for this are only the reads used to produce a trascriptome (see below).

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J PRNAscaf
#SBATCH --cpus-per-task=40
#SBATCH --output PRNAscaf.log
#SBATCH --mem 300Gb
#SBATCH --exclude=node117,node118


# input files
ASSEMBLY="$HOME/footloose_genome/genome_versions/S_parvus_wtdbg.ctg.polished1.fa.gz"
READ1="$HOME/footloose_genome/RNAdata/S_parvus_sub.1.fq"
READ2="$HOME/footloose_genome/RNAdata/S_parvus_sub.2.fq"
mkdir rna_scaffold
cd rna_scaffold

# This uses only the RNA seq reads that went in to the transcriptome assembly (60M total reads)

##### Maping with HiSat2
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n################## Indexing genome with HiSat2 ######################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"

hisat2-build -p 40 $ASSEMBLY ${ASSEMBLY}_hisat


printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n################# Alighing RNA reads with HiSat2 ####################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"

hisat2 --threads 40 -x ${ASSEMBLY}_hisat -1 $READ1 -2 $READ2 -k 3 -p 40 --pen-noncansplice 1000000 -S input.sam

printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#################### Scaffolding with P_RNA #######################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"
printf "\n\n\n#####################################################################\n\n\n"

P_RNA_scaffolder.sh -d $HOME/software/P_RNA_scaffolder -i input.sam -j $ASSEMBLY -F $READ1 -R $READ2  -t 40

echo Donezos
```

Results:

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | Scaffold N50 | Number of Scaffolds | BUSCO 
--- | --- | --- | --- | --- | --- | --- 
S_parvus.4.0.polished1 | 3,982,189,551 | 611,229 | 22,402 | 0.0 | C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950
S_parvus.4.0.polished1.RNAscaffolded (HiSat2) | 3,982,265,665 | 611,229 | 22,402 | 657,682 | 21,709 | C:82.5%[S:81.5%,D:1.0%],F:6.6%,M:10.9%,n:5310 
S_parvus.4.0.polished1.RNAscaffolded (BWA) | 3,982,222,951 | 611,229 | 22,402 | 627,725 | 22,068 | C:84.6%[S:83.6%,D:1.0%],F:5.8%,M:9.6%,n:5310
S_parvus.4.0.polished1.purged | 3,931,525,283 | 620,561 | 18,095 | NA | NA | C:74.5%[S:73.7%,D:0.8%],F:6.4%,M:19.1%,n:5310
S_parvus.4.0.polished1.purged.RNAscaffolded (HiSat2) | 3,931,525,283 | 620,561 | 18,095 | 667,429 | 17,410 | C:82.4%[S:81.5%,D:0.9%],F:6.6%,M:11.0%,n:5310
S_parvus.4.0.polished1.purged.RNAscaffolded (BWA) | 3,931,525,283 | 620,561 | 18,095 | 638,549 | 17,758 | C:81.7%[S:80.8%,D:0.9%],F:7.0%,M:11.3%,n:5310


I will be using this assembly moving forward: `S_parvus.4.0.polished1.RNAscaffoldedBWA.fa`

### Transcriptome assembly for foot-flaggin frog

The reads are DEEP. SO I subsmapled reads before doing an assembly.

```bash
seqtk sample -s13 BR2_CKDL210014905-1B_HJLVNDSX2_L1_1.fq.gz 20000000 > BR2_sub_1.fq
seqtk sample -s13 BR2_CKDL210014905-1B_HJLVNDSX2_L1_2.fq.gz 20000000 > BR2_sub_2.fq

seqtk sample -s13 Sp22_CKDL210014904-1B_HJLVNDSX2_L1_1.fq.gz 20000000 > SP22_sub_1.fq
seqtk sample -s13 Sp22_CKDL210014904-1B_HJLVNDSX2_L1_2.fq.gz 20000000 > SP22_sub_2.fq

seqtk sample -s13 LM2_CKDL210014907-1B_HJLVNDSX2_L1_1.fq.gz 20000000 > LM2_sub_1.fq
seqtk sample -s13 LM2_CKDL210014907-1B_HJLVNDSX2_L1_2.fq.gz 20000000 > LM2_sub_2.fq

# merge reads
cat *sub_1.fq > S_parvus_sub.1.fq
cat *sub_2.fq > S_parvus_sub.2.fq
```


WE used the Oyster River Protocol (a multi-kmer, multi-assembler approach) to produce a transcriptome across tissue types.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J ORP
#SBATCH --output ORP.log
#SBATCH --mem 600Gb
#SBATCH --cpus-per-task=40
module purge
module load anaconda/colsa

source activate orp-20191014

mkdir ORP_parvus
cd ORP_parvus

oyster.mk main \
TPM_FILT=1 \
MEM=575 \
CPU=40 \
READ1=$HOME/footloose_genome/RNAdata/S_parvus_sub.1.fq \
READ2=$HOME/footloose_genome/RNAdata/S_parvus_sub.2.fq \
RUNOUT=S_parvus_ORP
```



## Transcriptome results:

Transrate score: 0.41503
Transrate optimal score: 0.45624
BUSCO (metazoa??): C:93.1%[S:84.2%,D:8.9%],F:5.3%,M:1.6%,n:303
**NOTE:** run BUSCO with tetrapoda lineage.

OVerall, a good assembly. 

I want to look at only those assembled transcripts that have evidence suporting them. To do this, I used the `TransDecoder` pipeline with homoly evidence from `blastp` (with the uniprot database) and Pfam.

```

#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J TransDecoder
#SBATCH --cpus-per-task=40
#SBATCH --output TransDecoder.log
#SBATCH --mem 300Gb
#SBATCH --exclude=node117,node118


REF="$HOME/footloose_genome/ORP_parvus/assemblies/S_parvus_ORP.ORP.fasta"
REFNAME=$(basename $REF)
UNIPROT="$HOME/footloose_genome/uniprot_sprot.fasta"
PFAM_HMM="$HOME/footloose_genome/Pfam-A.hmm"


TransDecoder.LongOrfs -t $REF

# make blast database
makeblastdb -in $UNIPROT -parse_seqids -dbtype prot

# run blastp
blastp -query ${REFNAME}.transdecoder_dir/longest_orfs.pep  \
-db $UNIPROT  -max_target_seqs 1 \
-outfmt 5 -evalue 1e-5 -num_threads 40 > blastp.outfmt5

# build pfam binaries with hmmpress for use by hmmscan
hmmpress $PFAM_HMM

hmmscan --cpu 40 --domtblout pfam.domtblout $PFAM_HMM ${REFNAME}.transdecoder_dir/longest_orfs.pep

# TransDecoder.Predict -t target_transcripts.fasta [ homology options ]

TransDecoder.Predict -t $REF --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt5

```

I wanted to see how this transcriptome looks, so I ran BUSCO on it. 

```
cp S_parvus_ORP.ORP.fasta.transdecoder.cds S_parvus_ORP.ORP.fasta.transdecoder.cds.fasta
sbatch ~/scripts/txmetrics.job S_parvus_ORP.ORP.fasta.transdecoder.cds.fasta S_parvus_ORP.ORP.fasta.transdecoder
```

Total transcripts: 36217

Finally, I used these to annotate our assembly.

### Repeat modeling

Prior to annotating the assembly, it is useful to predict repeat regions. Especially in frogs since we know basically nothing. FOr this, I will use `RepeatModeler`.

```repeatmodeler.job
#!/bin/bash
#SBATCH --job-name=repeatmod
#SBATCH --output=repeatmodeler.log
#SBATCH --partition=macmanes
#SBATCH --cpus-per-task=40
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118

ASSEMBLY="$HOME/footloose_genome/genome_versions/S_parvus.4.0.polished1.RNAscaffoldedBWA.fa"
genome=$(basename $ASSEMBLY)

# env
module load linuxbrew/colsa

mkdir repeat_modeler
cd repeat_modeler

ln -s $ASSEMBLY

BuildDatabase -name $genome.repeatmodeler_db -engine ncbi $genome

RepeatModeler -database $genome.repeatmodeler_db -pa 40
```
Repeat database: `/mnt/oldhome/macmaneslab/ams1236/footloose_genome/repeat_modeler/RM_158667.ThuAug261505522021/consensi.fa.classified`

### Repeat masking


```repeatmasker.developedrepeatlibrary.job
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J repeatmask
#SBATCH --output repeatmasker.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

DIR=$(pwd)
ASSEMBLY="$HOME/footloose_genome/genome_versions/S_parvus.4.0.polished1.RNAscaffoldedBWA.fa"
genome=$(basename $ASSEMBLY)
PATH=/mnt/lustre/macmaneslab/macmanes/ncbi-blast-2.7.1+/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/shared/augustus_config/config


### input requires two fasta files. 1) the masked fasta file from RepeatMasker; 2) the genome assembly

# get most masked fasta file from RepeatModeler. This assumes that the most recent run is the correct one!
masked=$(ls repeat_modeler2/*db-families.fa)
maskedfile=$(basename $masked)

echo This is the masked fasta file from RepeatModeler $masked

## prep directory
mkdir repeatmasker
cd repeatmasker

# symlink everything here
ln -s ${DIR}/$masked
ln -s $ASSEMBLY

# sanity check
which perl

# make repeat library for parvus:
# first, I need to modify the database of repeats so that RepeatMasker picks up classes (I think, this is a test)
sed "s/\t/#/g" $HOME/repbase/RepBase25.05.fasta/vrtrep.ref | cut -d# -f 1-2 > $HOME/repbase/RepBase25.05.fasta/modifiedvertref.fa
echo Using repeat libary from parvus $masked
cat $HOME/repbase/RepBase25.05.fasta/modifiedvertref.fa ${DIR}/${masked} > $HOME/footloose_genome/vertrepeats_parvusRepModel.fa

# run repeatmasker
$HOME/software/RepeatMasker/RepeatMasker -pa 40 -gff -lib $HOME/footloose_genome/vertrepeats_parvusRepModel.fa -q $genome
```
