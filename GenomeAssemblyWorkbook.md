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

### Genome assembly comparisons

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.1.0 | 2,143,282,740 | 28,208 | 10,174 | 0.00 | C:0.6%[S:0.6%,D:0.0%],F:1.0%,M:98.4%,n:3950
S_parvus.2.0 | 3,782,941,384 | 156,120 | 49,172 | 0.00 | C:28.7%[S:28.6%,D:0.1%],F:14.2%,M:57.1%,n:3950
S_parvus.2.0_lowcoverage  | 3,756,961,657 | 197,812 | 42,434 | 0.00 | C:39.8%[S:39.3%,D:0.5%],F:15.2%,M:45.0%,n:3950
S_parvus.3.0 | 3,947,455,185 | 401,758 | 27,869 | 0.0 | C:56.5%[S:55.7%,D:0.8%],F:14.6%,M:28.9%,n:3950
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

I then reran this racon pipeline with the polished genome to produce a genome that had been polished 2x. This made the genome negligibly worse. 

We have gotten preads from PacBio. These are error corrected consensus reads that are produced by the FALCON assembler. PacBio ran FALCON for us using the data from our first 4 SMRTcells. This produced 75 TB of intermediate files, which...fuck that is so much. Anyway, I am going to use these error corrected reads to polish again see if I can't improve genic content a bit more.

```bash
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J racon
#SBATCH --output racon2.log
#SBATCH --cpus-per-task=24
#SBATCH --exclude=node117,node118
#SBATCH --mem=300000
set -x

DIR=$(pwd)
ASSEMBLY="$HOME/footloose_genome/racon/S_parvus_wtdbg.ctg.polished1.fa"
genome=$(basename $ASSEMBLY)
READS=$"$HOME/footloose_genome/raw_data/FromPacBioNov2020/falcon_assembly_results/preads4falcon.unwrapped.fasta"


module load linuxbrew/colsa

# preparation
mkdir racon_preads
cd racon_preads

#cp $ASSEMBLY .

awk '{print $1}' $genome > new.fasta
mv new.fasta $genome

### First align reads with minimap2
echo aligning with minimap2
minimap2 -I10G -t 40 -x asm20 $genome $READS | gzip -c - > Footlose.PB.paf.gz

### Run racon
echo Polishing with racon
racon -t 40 $READS Footlose.PB.paf.gz $genome > S_parvus_wtdbg.ctg.polished2.fa

```



### Genome assembly comparisons

Assembly | Genome Size (GB) | Contig N50 | Number of Contigs | %Ns | BUSCO 
--- | --- | --- | --- | --- | --
S_parvus.1.0 | 2,143,282,740 | 28,208 | 10,174 | 0.00 | C:0.6%[S:0.6%,D:0.0%],F:1.0%,M:98.4%,n:3950
S_parvus.2.0 | 3,782,941,384 | 156,120 | 49,172 | 0.00 | C:28.7%[S:28.6%,D:0.1%],F:14.2%,M:57.1%,n:3950
S_parvus.2.0_lowcoverage  | 3,756,961,657 | 197,812 | 42,434 | 0.00 | C:39.8%[S:39.3%,D:0.5%],F:15.2%,M:45.0%,n:3950
S_parvus.3.0 | 3,947,455,185 | 401,758 | 27,869 | 0.0 | C:56.5%[S:55.7%,D:0.8%],F:14.6%,M:28.9%,n:3950
S_parvus.4.0 | 3,985,797,087 | 608,715 | 23,418 | 0.0 | C:68.3%[S:67.2%,D:1.1%],F:12.6%,M:19.1%,n:3950
S_parvus.4.0.polished1 | 3,982,189,551 | 611,229 | 22,402 | 0.0 | C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950
S_parvus.4.0.polished2 | 3,918,212,060 | 609,425 | 20,281 | 0.0 | C:73.3%[S:72.0%,D:1.3%],F:11.5%,M:15.2%,n:3950
S_parvus.4.0.preadspolished1 | 3,963,679,315 | 605,739 | 23,418 |  C:72.9%[S:71.6%,D:1.3%],F:11.5%,M:15.6%,n:3950

It seems like the second round of polishing did not help. We have been in contact with PacBio about all of this, and they have graciously provided some corrected reads (preads) for us. I used those preads to polish the assembly again, and that did not help. So I am choosing to move forward with the assembly after it has been polished 1x by Racon. I then used the program `purge_dups` to...purge dups.

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
S_parvus.1.0 | 2,143,282,740 | 28,208 | 10,174 | 0.00 | C:0.6%[S:0.6%,D:0.0%],F:1.0%,M:98.4%,n:3950
S_parvus.2.0 | 3,782,941,384 | 156,120 | 49,172 | 0.00 | C:28.7%[S:28.6%,D:0.1%],F:14.2%,M:57.1%,n:3950
S_parvus.2.0_lowcoverage  | 3,756,961,657 | 197,812 | 42,434 | 0.00 | C:39.8%[S:39.3%,D:0.5%],F:15.2%,M:45.0%,n:3950
S_parvus.3.0 | 3,947,455,185 | 401,758 | 27,869 | 0.0 | C:56.5%[S:55.7%,D:0.8%],F:14.6%,M:28.9%,n:3950
S_parvus.4.0 | 3,985,797,087 | 608,715 | 23,418 | 0.0 | C:68.3%[S:67.2%,D:1.1%],F:12.6%,M:19.1%,n:3950
S_parvus.4.0.polished1 | 3,982,189,551 | 611,229 | 22,402 | 0.0 | C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950
S_parvus.4.0.polished2 | 3,918,212,060 | 609,425 | 20,281 | 0.0 | C:73.3%[S:72.0%,D:1.3%],F:11.5%,M:15.2%,n:3950
S_parvus.4.0.preadspolished1 | 3,963,679,315 | 605,739 | 23,418 |  C:72.9%[S:71.6%,D:1.3%],F:11.5%,M:15.6%,n:3950
S_parvus.4.0.polished1.purged | 3,931,525,283 | 620,561 | 18,095 |  C:74.4%[S:73.2%,D:1.2%],F:11.3%,M:14.3%,n:3950




I ran the purged assembly through BUSCO v4.1., using tetrapoda_odb10 (eukaryota, 2020-09-10

C:74.5%[S:73.7%,D:0.8%],F:6.4%,M:19.1%,n:5310





### Prepping publicly available RNA seq data for genome annotation purposes

First up, downloading RNA seq reads to use for annotation. We don't have any ourselves, so we will use RNA seq reads from the American bullfrog.

Tissue | Treatment | Link
--- | --- | ---
Whole tadpole | NA | https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=DRR040619
Olfactory bulb | NA | https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR4048907
Back skin | NA | https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR4045340
Lung sample | NA | https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR4045335
Tail fin | NA | https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR4045338
Brain | NA | https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR4044081


Downloading reads:

```bash
fasterq-dump $ACCESSION NUMBER
```

Then I had to fix header files:

```bash
# add this code here

```

Finally, I assembled a _de novo_ transcriptome for individuals from both studies. I did this to avoid some heterozygosity issues between studies (though depending on the source of the different regions of the bullfrog from the second study this may still be a big issue).

Assembly for the whole bullfrog tadpole:

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J ORPtad
#SBATCH --output ORPtad.log
####SBATCH --mem 300Gb
#SBATCH --cpus-per-task=24
module purge
module load anaconda/colsa

source activate orp-20191014

mkdir ORPtad
cd ORPtad

oyster.mk main \
TPM_FILT=1 \
MEM=100 \
CPU=24 \
READ1=$HOME/footloose_genome/raw_data/RNAseqData/DRR040619_1.fastq \
READ2=$HOME/footloose_genome/raw_data/RNAseqData/DRR040619_2.fastq \
RUNOUT=Bullfrog_tad_ORP
```

Assembly for the various tissues from the other study:

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J ORPmulti
#SBATCH --output ORPmulti.log
#SBATCH --mem 300Gb
#SBATCH --cpus-per-task=40
module purge
module load anaconda/colsa

source activate orp-20191014

# cat reads
cat raw_data/RNAseqData/SRR*1.fastq > raw_data/RNAseqData/combined.R1.fq
cat raw_data/RNAseqData/SRR*2.fastq > raw_data/RNAseqData/combined.R2.fq

delete files for quota
rm raw_data/RNAseqData/SRR*

# subsample to 40 million reads
seqtk sample -s 13 raw_data/RNAseqData/combined.R1.fq.gz 40000000 > raw_data/RNAseqData/multitissue.subsamp.R1.fq
seqtk sample -s 13 raw_data/RNAseqData/combined.R2.fq.gz 40000000 > raw_data/RNAseqData/multitissue.subsamp.R2.fq

delete files for quota
# rm raw_data/RNAseqData/SRR*

# prep
mkdir ORPmulti
cd ORPmulti


# ahw ya from southie? ya want some oystas??
oyster.mk main \
TPM_FILT=1 \
MEM=100 \
CPU=24 \
SPADES1_KMER=61 \
SPADES2_KMER=41 \
READ1=$HOME/footloose_genome/raw_data/RNAseqData/multitissue.subsamp.R1.fq \
READ2=$HOME/footloose_genome/raw_data/RNAseqData/multitissue.subsamp.R2.fq \
RUNOUT=Bullfrog_multitissues_ORP

```

I then analyzed these with TransRate and BUSCO using the tetrapoda database.

## Transcriptome results:

**TABLE HERE ONCE COMPLETE**

Finally, I used these to annotate our assembly.

### Annotation

First, a rousing round of repeat modeler:

```
#!/bin/bash
#SBATCH --job-name=repeatmod
#SBATCH --output=repeatmodeler.log
#SBATCH --partition=macmanes
#SBATCH --cpus-per-task=40
#SBATCH --open-mode=append
#SBATCH --exclude=node117,node118

ASSEMBLY="$HOME/footloose_genome/genome_versions/S_parvus_wtdbg.ctg.polished1.purged.fa"
genome=$(basename $ASSEMBLY)

# env
module load linuxbrew/colsa

mkdir repeat_modeler
cd repeat_modeler

ln -s $ASSEMBLY

BuildDatabase -name $genome.repeatmodeler_db -engine ncbi $genome

RepeatModeler -database $genome.repeatmodeler_db -pa 40
```

After this, I ran RepeatMasker.

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J repeatmask
#SBATCH --output repeatmasker.log
#SBATCH --cpus-per-task=40
#SBATCH --exclude=node117,node118

module load linuxbrew/colsa

DIR=$(pwd)
ASSEMBLY="$HOME/footloose_genome/genome_versions/S_parvus_wtdbg.ctg.polished1.purged.fa"
genome=$(basename $ASSEMBLY)
PATH=/mnt/lustre/macmaneslab/macmanes/ncbi-blast-2.7.1+/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/mnt/lustre/macmaneslab/shared/augustus_config/config


### input requires two fasta files. 1) the masked fasta file from RepeatMasker; 2) the genome assembly

# get most masked fasta file from RepeatModeler. This assumes that the most recent run is the correct one!
masked=$(ls -ltr repeat_modeler/RM_*/consensi.fa.classified | tail -n1 | awk '{print $NF}')
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

# run repeatmasker
$HOME/software/RepeatMasker/RepeatMasker -pa 40 -gff -lib $maskedfile  -q $genome


```

```
Repeat Masker results:

file name: S_parvus_wtdbg.ctg.polished1.purged.fa
sequences:         18095
total length: 3931525283 bp  (3931525237 bp excl N/X-runs)
GC level:         42.49 %
bases masked: 2430076406 bp ( 61.81 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
Retroelements       857352    442929408 bp   11.27 %
   SINEs:                0            0 bp    0.00 %
   Penelope         132303     48927208 bp    1.24 %
   LINEs:           744936    349034810 bp    8.88 %
    CRE/SLACS            0            0 bp    0.00 %
     L2/CR1/Rex     380312    164515603 bp    4.18 %
     R1/LOA/Jockey       0            0 bp    0.00 %
     R2/R4/NeSL          0            0 bp    0.00 %
     RTE/Bov-B        1290       458127 bp    0.01 %
     L1/CIN4        218200    125501718 bp    3.19 %
   LTR elements:    112416     93894598 bp    2.39 %
     BEL/Pao          6278      6722189 bp    0.17 %
     Ty1/Copia        3510      2316667 bp    0.06 %
     Gypsy/DIRS1     92726     73059855 bp    1.86 %
       Retroviral     6993      8490370 bp    0.22 %

DNA transposons     865064    493406526 bp   12.55 %
   hobo-Activator    61742     50885299 bp    1.29 %
   Tc1-IS630-Pogo   703957    389733641 bp    9.91 %
   En-Spm                0            0 bp    0.00 %
   MuDR-IS905            0            0 bp    0.00 %
   PiggyBac          93686     47833838 bp    1.22 %
   Tourist/Harbinger     0            0 bp    0.00 %
   Other (Mirage,        0            0 bp    0.00 %
    P-element, Transib)

Rolling-circles          0            0 bp    0.00 %

Unclassified:       5255637   1401539996 bp   35.65 %

Total interspersed repeats:  2337875930 bp   59.46 %


Small RNA:               0            0 bp    0.00 %

Satellites:              0            0 bp    0.00 %
Simple repeats:     818707     77722625 bp    1.98 %
Low complexity:      97081     14477851 bp    0.37 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element


RepeatMasker Combined Database: Dfam_3.1

run with rmblastn version 2.2.27+
The query was compared to classified sequences in "consensi.fa.classified"
```



Next, a grueling round of Maker. For Maker config files please see the [maker_data/](https://github.com/AdamStuckert/Footloose_genome/tree/master/maker_data) directory. I ran a single round of Maker, without running RepeatMasker within Maker (as I had already run both RepeatModeler and RepeatMasker). I used transcript evidence from the American bullfrog (Ranidae, *Lithobates catesbeianus*) to annotate this genome. These transcripts were produced via the Oyster River Protocol (code above) from publically available reads. 

Maker script:

```
#!/bin/bash
#SBATCH --partition=macmanes,shared
#SBATCH -J maker
#SBATCH --ntasks=160
#SBATCH --mem 110Gb
#SBATCH --output maker.masked.log
#SBATCH --open-mode=append

set -x

module purge
module load anaconda/colsa
source activate maker-3.01.02
module purge
module load linuxbrew/colsa


species="Staurois_parvus"
genome="$HOME/footloose_genome/repeatmasker/S_parvus_wtdbg.ctg.polished1.purged.fa.masked"
gbase=$(basename $genome | sed "s/.fasta//g" | sed "s/.fa.masked//g")
PREFIX="S_parv"
OUTPUT="maker"
MAKER="/mnt/lustre/macmaneslab/macmanes/test/maker/bin/"

mkdir -p $HOME/footloose_genome/${OUTPUT}
cd $HOME/footloose_genome/${OUTPUT}

echo Running maker on $genome
echo Maker output going to $OUTPUT

echo Using control files:
echo "$HOME/footloose_genome/maker_data/maker_opts_uniprot_NoRepeatMaskerWithinMaker.ctl"
echo "$HOME/footloose_genome/maker_data/maker_bopts.ctl"
echo "$HOME/footloose_genome/maker_data/maker_exe.ctl"

mpiexec -n 160 /mnt/lustre/macmaneslab/macmanes/test/maker/bin/maker \
-fix_nucleotides -base "$species" -quiet \
-genome "$genome" \
$HOME/footloose_genome/maker_data/maker_opts_uniprot_NoRepeatMaskerWithinMaker.ctl \
$HOME/footloose_genome/maker_data/maker_bopts.ctl \
$HOME/footloose_genome/maker_data/maker_exe.ctl

$MAKER/fasta_merge -d "$species".maker.output/"$species"_master_datastore_index.log -o "$species"."$gbase"
$MAKER/gff3_merge -d "$species".maker.output/"$species"_master_datastore_index.log -o "$species"."$gbase".gff3 -n
lastal -P22 /mnt/lustre/macmaneslab/macmanes/transporters/swissprot "$species"."$gbase".all.maker.proteins.fasta -f BlastTab > blast.out
$MAKER/maker_functional_fasta /mnt/lustre/macmaneslab/macmanes/transporters/uniprot_sprot.fasta blast.out "$species"."$gbase".all.maker.proteins.fasta > "$species"."$gbase".functional.proteins.fasta
$MAKER/maker_functional_fasta /mnt/lustre/macmaneslab/macmanes/transporters/uniprot_sprot.fasta blast.out "$species"."$gbase".all.maker.transcripts.fasta > "$species"."$gbase".functional.transcripts.fasta
MAKER/maker_functional_gff /mnt/lustre/macmaneslab/macmanes/transporters/uniprot_sprot.fasta blast.out "$species"."$gbase".gff3 > "$species"."$gbase".functional.gff3
MAKER/maker_map_ids --prefix "$PREFIX"_ --justify 6 "$species"."$gbase".functional.gff3 > "$species"."$gbase".genome.all.id.map
MAKER/map_fasta_ids "$species"."$gbase".genome.all.id.map  "$species"."$gbase".functional.proteins.fasta
MAKER/map_gff_ids "$species"."$gbase".genome.all.id.map  "$species"."$gbase".functional.gff3
$MAKER/map_fasta_ids "$species"."$gbase".genome.all.id.map  "$species"."$gbase".functional.transcripts.fasta

# get annotation information for RNAseq analyses
grep "^>" "$species"."$gbase".functional.transcripts.fasta | tr -d ">" > headers.txt
awk '{print $1}' headers.txt  > transcripts.txt
cut -f 2 -d '"' headers.txt  | sed "s/Similar to //g" > annotations.txt
paste transcripts.txt annotations.txt > "$species"."$gbase".annotations.tsv

```

Note: for now I'm actually annotating without repeatmodeler/repeatmasker (using maker_opts_uniprot.ctl) until I convert repeat annotations from gff2 to gff3. Conda install of the software I want has been slow.

Trying to get a functioning gff3 file:

```bash
# create GFF3 (from Daren Card's GitHub **add link**)
./rmOutToGFF3custom.sh -o S_parvus_wtdbg.ctg.polished1.purged.fa.out > S_parvus_wtdbg.ctg.polished1.purged.fa.out.gff3
# isolate complex repeats
grep -v -e "Satellite" -e "Low_complexity" -e "Simple_repeat" S_parvus_wtdbg.ctg.polished1.purged.fa.out.gff3 \
  > S_parvus_wtdbg.ctg.polished1.purged.fa.out.complex.gff3
# reformat to work with MAKER
cat S_parvus_wtdbg.ctg.polished1.purged.fa.out.complex.gff3 | \
  perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
  > S_parvus_wtdbg.ctg.polished1.purged.fa.out.complex.reformat.gff3
```

Quick note: Newer versions of Uniprot databases break this pipeline for some reason...which is quite annoying. Anyway, this works if I use an older version! I will update this with the path to that fasta file later.
