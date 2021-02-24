# Comparative Genomics

This details our approach to conducting comparative genomics. Broadly speaking, we are examining sequences for a core group of genes that we think may be playing a role in the ability to do rapid foot flagging in these frogs. Many of these are from previous work by Mangiamele et al.

## Identifying sequences for candidate genes

We are going to extract the sequences of the genes we want. As an easy first pass, I am going to extract the full sequence that has been identified as a "gene" by Maker. To do this, and to cut down on some compute time, I subsetted the annotation file to extract only the lines labelled as "gene".

> awk '$3 == "gene" {print $0}' Staurois_parvus.S_parvus_wtdbg.ctg.polished1.purged.fa.functional.gff3 > Staurois_parvus.S_parvus_wtdbg.ctg.polished1.purged.genesonly.gff3

I then wrote a small script that finds all instances of the genes we want, and extracts the sequences into individual fasta files:

```bash
#!/bin/bash extractgenesofinterest.sh

usage ()
{

    Extract all the genes I want from a genome.

    Note: For this I am using a modified gff3 file which includes only what Maker has anotated as a "gene". I created this gene gff3 file using awk '$3 == "gene" {print $0}' Staurois_parvus.S_parvus_wtdbg.ctg.polished1.purged.fa.functional.gff3 > Staurois_parvus.S_parvus_wtdbg.ctg.polished1.purged.genesonly.gff3
}

GFF="$HOME/footloose_genome/maker_rerun/Staurois_parvus.S_parvus_wtdbg.ctg.polished1.purged.genesonly.gff3"
GENOME="$HOME/footloose_genome/genome_versions/S_parvus_wtdbg.ctg.polished1.purged.fa"

mkdir genesearchoutput
cd genesearchoutput/

# clean up previous run
rm geneinformation.txt

genes="ar
bdnf
ntf3
ntf4
ntrk2
ngf
gja1
gjb6
foxp2
foxp1
gli3
ctnnb1
wnt
hoxc10
efna1
wnt11
gsk3b
apc
lef1
lrp5
lrp6
wnt8a
fzd3
fzd7
wnt5a
wnt5b
wnt6
wnt4"

#### Find each gene
echo Searching through $genes


for gene in $genes
do
grep -i "Similar to $gene:" $GFF >> $gene.txt
done



##### Make a usable bed file for each gene
for gene in $genes
do
awk -v OFS="\t" -F "\t" '{print $1, $4, $5, $9}' $gene.txt > tmp
sed -i "s/ID\\d.*Note=//g" tmp    #### for some reason this is not working, but whatever
sed -i "s/\s(.*//g" tmp
mv tmp $gene.bed
done



#for gene in $genes
#do
#info=$(grep -i $gene $GFF)
#contig=$(grep -i $gene $GFF | cut -f1)
#start=$(grep -i $gene $GFF | cut -f4 -d " ")
#end=$(grep -i $gene $GFF | cut -f5 -d " ")
#geneinfo=$(grep -i $gene $GFF | cut -f4 -d ";" | sed "s/Note=Similar to //g" | cut -f1 -d "(")
#printf "%s\t%s\t%s\t%s\n" "$contig" "$start" "$end" "$geneinfo" > $gene.bed
#printf "%s\n" "$info" >> geneinformation.txt
#done




#for gene in $genes
#do
#info=$(grep -i $gene $GFF)
#contig=$(echo $info | cut -f1 -d " ")
#start=$(echo $info | cut -f4 -d " ")
#end=$(echo $info | cut -f5 -d " ")
#printf "%s\t%s\t%s\tS_parvus%s\n" "$contig" "$start" "$end" "$gene" > $gene.bed
#printf "%s\n" "$info" >> geneinformation.txt
#done



### Extract each gene sequence

printf "##################################\n"
printf "##################################\n"
printf "#### Extracting each gene seq ####\n"
printf "##################################\n"
printf "##################################\n"



for bed in $(ls *bed)
do
bedname=$(echo $bed | sed "s/.bed//")
bedtools getfasta -fi $GENOME -bed $bed -name > $bedname.fa
done

```

### Identifying taxa to include in our comparative genomic approach

I queried NCBI's genomes on Feb 24, 2021 for all amphibians. Results:

```
Ambystoma mexicanum
Bufo bufo
Bufo gargarizans
Geotrypetes seraphini
Leptobrachium leishanense
Limnodynastes dumerilii
Lithobates catesbeianus
Microcaecilia unicolor
Nanorana parkeri
Oophaga pumilio
Pyxicephalus adspersus
Rana temporaria
Rhinatrema bivittatum
Rhinella marina
Scaphiopus couchii
Scaphiopus holbrookii
Spea bombifrons
Spea multiplicata
Xenopus laevis
Xenopus tropicalis
```

We will use the green anole and chicken as outgroups. 

I also queried all of the manakin genera (Antilophia, Ceratopipra, Chiroxiphia, Chloropipo, Corapipo, Cryptopipo, Heterocercus, Ilicura, Lepidothrix, Machaeropterus, Manacus, Masius, Neopelma, Pipra, Pseudopipra, Tyranneutes, Xenopipo). JK. The website is shit and you can't search just by genus name! Nonsense.
## Pulling out sequences from other organisms.


This code will pull out our genes of interest from a variety of taxa. First up, using an R script to identify the taxa ID of all of our groups.

```R
library(rentrez)
library(XML)
# If you are using many many searches and want to go faster than 3/second, you will need a key. https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/   
#set_entrez_key("YOURKEY") # I did not do this because I wrote a short pause in my function

## function to pull out taxa IDs
fetch_taxid = function(x){
  ret = NA
  tryCatch({
    search = entrez_search(db="taxonomy",term = x,)
    fetch = xmlToDataFrame(entrez_fetch(db="taxonomy",id=search$ids,rettype="xml"))
    ret=fetch$TaxId 
  }, error = function(e) {
    ret = NA
  })
  return(ret)
}

taxas=c(
"Ambystoma mexicanum",
"Bufo bufo",
"Bufo gargarizans",
"Geotrypetes seraphini",
"Leptobrachium leishanense",
"Limnodynastes dumerilii",
"Lithobates catesbeianus",
"Microcaecilia unicolor",
"Nanorana parkeri",
"Oophaga pumilio",
"Pyxicephalus adspersus",
"Rana temporaria",
"Rhinatrema bivittatum",
"Rhinella marina",
"Scaphiopus couchii",
"Scaphiopus holbrookii",
"Spea bombifrons",
"Spea multiplicata",
"Xenopus laevis",
"Xenopus tropicalis",
"anolis carolinensis", 
"chicken")

ID.DF <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c("taxa", "entrez_ID"))))
for (tax in taxas){
  ID <- fetch_taxid(tax)
  tmp <- c(tax, ID)
  ID.DF <- rbind(ID.DF, tmp)
  Sys.sleep(1) # this is so that we don't get interrupted due to too many searches
}

# rename colnames in df
colnames(ID.DF) <- c("taxa", "entrez_ID")

# save the datas
write.table(ID.DF, "taxaIDs.tsv", sep = "\t", row.names = FALSE)
```

### Downloading sequences, aligning, and creating gene trees

I'll now switch to bash to use eutils to actually pull out sequences. Currently, my approach is a bit rough, but works. I'll probably need to do a bit of a cleanup later. The general approach is to download all the genes of interest from our taxa of choice using `eutils`, pull out only the first sequence, create a single fasta file per gene (**I probably need to add a line to rename the sequences to something easier to read/interpret**), run mafft to align, and iqtree to create a tree.

Script:

```bash
#!/bin/bash
#SBATCH --job-name=compgen
#SBATCH --output=compgen.log
#SBATCH --cpus-per-task=24
#SBATCH --partition=macmanes,shared
#SBATCH --exclude=node117,node118


# prep
module purge
ml anaconda/colsa
conda activate busco-5.beta

# Now to do a big, systematic search.
mkdir comparativegenomics
mkdir comparativegenomics/taxaspecificfastas

printf "########################################\n###### Extracting gene sequences ######\n########################################\n"
# for loop for the genes we want to search:
for gene in $(cat comparativegenomics/Genes2Search.txt)
do

# inner for loop to pull in taxonomy ids
for taxa in $(cat comparativegenomics/taxaIDs.tsv) # note the way this markdown is written it isn't technically correct because I've just extracted the TaxaID column from the R script
do

# search
printf "Searching for %s in %s\n" "$gene" "$taxa"
esearch -db nucleotide -query "txid$taxa AND $gene[gene]" | efetch -format fasta > comparativegenomics/taxaspecificfastas/${gene}_${taxa}.fa

done
done

# Extract one sequence per gene x taxa combo:
printf "########################################\n###### Extracting one sequence per gene ######\n########################################\n"

```bash
for fasta in $(ls comparativegenomics/taxaspecificfastas/*fa)
do
# unwrap fastas
unwrap_fasta $fasta tmp.fa
# extract only first seq
head -n2 tmp.fa > $fasta
done

# Merge all the fasta files:

printf "########################################\n###### Merging fasta files ######\n########################################\n"

mkdir comparativegenomics/collatedgenefastas
for gene in $(cat comparativegenomics/Genes2Search.txt)
do
cat genesearchoutput/${gene}.fa > comparativegenomics/collatedgenefastas/${gene}.fa
cat $(ls comparativegenomics/taxaspecificfastas/${gene}*fa) >> comparativegenomics/collatedgenefastas/${gene}.fa



done

# Align with mafft.

printf "########################################\n###### Aligning sequences with mafft ######\n########################################\n"

module purge
ml linuxbrew/colsa
mkdir comparativegenomics/alignments
for fasta in $(ls comparativegenomics/collatedgenefastas/*fa)
do
# get gene name
gene=$(basename $fasta | sed "s/.fa//")
# run mafft to align
mafft --reorder --localpair --thread 24 $fasta  > comparativegenomics/alignments/$gene.aln
done


# View alignments:

# Navigate here: https://www.ebi.ac.uk/Tools/msa/mview/

# Make a tree for each gene:

printf "########################################\n###### Making gene trees ######\n########################################\n"

mkdir comparativegenomics/treefiles
for aln in $(ls comparativegenomics/alignments/*aln)
do
gene=$(basename $aln | sed "s/.aln//")
printf "########################################\n###### Producing tree file on $gene ######\n########################################\n"
iqtree -s $aln -st DNA -m HKY -nt 24 -bb 1000 -pre comparativegenomics/treefiles/${gene} # add -o for outgroup!
done
```

Ok, all of this preliminarily works! Next steps:

1. Identify all the taxa to include in the real version
2. Make it a slurm script, add in threads (24)
3. Verify best practices/choices in this pipeline
4. What are the next steps? How do we make inferences from this, how do we visualize these data?



