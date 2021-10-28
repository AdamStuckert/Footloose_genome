# Genome repeat masking pipeline 

The general approach:

1. RepeatModeler to produce species-specific TE libraries
2. Softmask the genome with RepeatMasker


## General information about downloads, etc

I downloaded `RepeatModeler2` and relevant programs in a conda environment (boringly called `repeat_modeler2`). This has everything to run these scripts, except for the newest version of RepeatMasker. I manually downloaded this.

RepeatMasker needs to be configured before use.

```bash
# download DFAM
cd RepeatMasker/Libraries
wget https://www.dfam.org/releases/Dfam_3.4/families/Dfam.h5.gz
gunzip Dfam.h5.gz
```

I also included RepBase data (`RepBase25.05`) by running ` perl addRepBase.pl` from within the RepeatMasker folder (the RepBas embl files were in the `Libraries` folder).

Then I configured `RepeatMasker`

```bash
perl ./configure

# I configured RMblast to be the default, which uses the DFAM and RepBase datasets.
# I also configured HMMER
```

Get Dipteran-specific repeats only:

`./famdb.py -i Libraries/RepeatMaskerLib.h5 families --format fasta_name --ancestors --descendants "vertebrata" --include-class-in-name > Libraries/vertebrata.fa`


I then ran RepeatModeler2.

This runs RepeatModeler to make species specific (predicted) TE libraries.

```bash RepMod.job
#!/bin/bash
#SBATCH -p macmanes
#SBATCH -J RepMod
##SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH -t 4-0:00:00
#SBATCH --output RepMod_%j.log

# note: submit from repeat_modeler2 conda env
module purge

REF="$HOME/Lutzo_assemblies/LUT06.asm_210301.fasta"
SPECIES="Staurois_parvus"

BuildDatabase -name $SPECIES.repeatmodeler_db -engine ncbi $REF

RepeatModeler -database $SPECIES.repeatmodeler_db -pa 40

```

