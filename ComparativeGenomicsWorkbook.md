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
