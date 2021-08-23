
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
