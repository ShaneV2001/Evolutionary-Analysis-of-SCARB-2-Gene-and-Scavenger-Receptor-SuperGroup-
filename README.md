# Project Repo: Analysis on SCARB-2 Gene (XP_021375881)



## Preamble
SCARB-2 Gene is named "project" for simplicity in editing command lines. All references and version numbers will be located in the Papers Citations. 

*WORK ON ADJUSTING CONTENTS HYPERLINKS*

# 1. Homolog Collection with BLAST

Used NCBI to obtain the BLAST data to collect preliminary data using the accession reference sequence XP_021373098.1

## Commands to Aligning Homologs
Here we downloaded our protein in FASTA format
```
ncbi-acc-download -F fasta -m protein XP_021373098.1
```
Here performed a BLAST search using our protein and are creating a table output for our data
```
blastp -db ../allprotein.fas -query XP_021373098.1.fa -outfmt 0 -max_hsps 1 -out project.blastp.typical.out
blastp -db ../allprotein.fas -query XP_021373098.1.fa -outfmt "6 sseqid pident length mismatch gappen evalue bitscore pident stitle"  -max_hsps 1 -out project.blastp.detail.out
```
Here we used selected our e value cut off. We used an e-value cutoff of 1e-35. We decided to use this e-value to have a high sequence identity, and a high alignment score as well as minimizes false homologs. This results in a filter BLAST sequences
```
awk '{if ($6< 1e-35 )print $1 }' project.blastp.detail.out > project.blastp.detail.filtered.out
```
This was used to count the number of hits after filtering 
```
wc -l project.blastp.detail.filtered.out
```
________________________________________________________________________________________________________________________________________________________________________________
# 2. Sequence Alignment 

Using the filtered BLAST results, we can use MUSCLE to align all the sequences with each other. MUSCLE will create multiple sequence alignments based on different alignment accuracy benchmarks. We used Alv as a visual aid for our alignment. We used alignbuddy to assist with calculating the lengths of the alignments. We used T-Coffee to calculate percent identities among the sequences

## Commands to Align Sequences
We use Seqkit to obtain the BLAST sequences 
```
seqkit grep --pattern-file ~/labs/lab3-$MYGIT/project/project.blastp.detail.filtered.out ~/labs/lab3-$MYGIT/allprotein.fas > ~/labs/lab4-$MYGIT/project/project.homologs.fas
```
This is the muscle command used to create the multiple sequence alignments 
```
muscle -in ~/labs/lab4-$MYGIT/project/project.homologs.fas -out ~/labs/lab4-$MYGIT/project/project.homologs.al.fas
```
Here we use the Alv command to view the alignment and save the alignment in a pdf format 
```
alv -kli --majority ~/labs/lab4-$MYGIT/project/project.homologs.al.fas | less -RS
```
This is how we formatted our alignments on alv. 
```
muscle -in ~/labs/lab4-$MYGIT/project/project.homologs.fas -html -out ~/labs/lab4-$MYGIT/project/project.homologs.al.html
sed -i 's/<PRE>/<pre style="font-family: 'FreeMono', monospaced;">/g' project.homologs.al.html
wkhtmltopdf ~/labs/lab4-$MYGIT/project/project.homologs.al.html --print-media-type ~/labs/lab4-$MYGIT/project/project.homologs.al.pdf
```
This is we begin to use alignbuddy. The First command line provides the width(length) of the alignment. The second command line calculates the alignment after removing any gaps. The last command line calculates the alignment after removing invariant positions.
```
alignbuddy  -al  ~/labs/lab4-$MYGIT/project/project.homologs.al.fas
alignbuddy -trm all  ~/labs/lab4-$MYGIT/gqr/gqr.homologs.al.fas | alignbuddy  -al
alignbuddy -dinv 'ambig' ~/labs/lab4-$MYGIT/gqr/gqr.homologs.al.fas | alignbuddy  -al
```
Here we calculated the average percent identity using the t_coffee program 
```
 alignbuddy -pi ~/labs/lab4-$MYGIT/project/project.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
________________________________________________________________________________________________________________________________________________________________________________

# 3. Gene Family Sequence Alignment

For the Gene Family Alignment, we used IQ tree to create the optimal phylogenetic trees based on our previously obtained alignments. 

## Commands for the Gene Family Alignments
Here we used IQ-tree with ultra-fast bootstrap support levels using the -bb 1000. This procedure takes about 10 minutes to run.
```
iqtree -s ~/labs/lab5-$MYGIT/project/project.homologs.al.fas -bb 1000 -nt 2 
```
Here we used gotree to create a midpoint rooted tree, since iq tree created an unrooted tree.
```
gotree reroot midpoint -i ~/labs/lab5-$MYGIT/project/project.homologs.al.fas.treefile -o ~/labs/lab5-$MYGIT/project/project.homologs.al.mid.treefile
```
Here we can see the rooted tree using the nw_display program
```
nw_order -c n ~/labs/lab5-$MYGIT/project/project.homologs.al.mid.treefile  | nw_display -
```

# 4. Tree Reconciliation and Rearrangement 

We used the midpoint-rooted gene tree produced in the last section to reconcile the genes using the program Notung. Notung will estimate duplication events and loss events. This allows us to differentiate between orthologs and paralogs. 

## Commands to Reconcile gene and species trees 
This command performs the reconciliation
```
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ~/labs/lab5-$MYGIT/species.tre -g ~/labs/lab6-$MYGIT/project/project.homologs.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/labs/lab6-$MYGIT/project/
```
We used this command to assign the internal nodes 
```
grep NOTUNG-SPECIES-TREE ~/labs/lab6-$MYGIT/project/project.homologs.al.mid.treefile.reconciled | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
```
We used this command to generate a RecPhyloXML to be able to view the gene-within-species using the ThirdKind program
```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/labs/lab6-$MYGIT/project/project.homologs.al.mid.treefile.reconciled --include.species
thirdkind -Iie -D 40 -f ~/labs/lab6-$MYGIT/project/project.homologs.al.mid.treefile.reconciled.xml -o  ~/labs/lab6-$MYGIT/project/project.homologs.al.mid.treefile.reconciled.svg
```
_________________________________________________________________________________________________________________________________________________________________________

# 5. Protein Domain Prediction

Protein Domain Prediction

## Commands to Aligning Homologs
Here we used RPS-BLAST to identify Pfam Domains without our specific protein sequences. We will be using the project.homolog.fas file containing the original, unaligned sequences from Section 2.
```
sed 's/*//' ~/labs/lab4-$MYGIT/project/project.homologs.fas > ~/labs/lab8-$MYGIT/project/project.homologs.fas
```
Following this, we download the Pfam database from NCBI. Then run the RPS blast command.
```
rpsblast -query ~/labs/lab8-$MYGIT/project/project.homologs.fas -db ~/data/Pfam -out ~/labs/lab8-$MYGIT/project/project.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
Here we used this Rscript to run R without opening an R console. plotTreeandDomain is a program developed by Dr. Rest. We used the data from the tree file and the rps-BLAST output.
```
sudo /usr/local/bin/Rscript  --vanilla ~/labs/lab8-$MYGIT/plotTreeAndDomains.r ~/labs/lab5-$MYGIT/project/project.homologs.al.mid.treefile ~/labs/lab8-$MYGIT/project/project.rps-blast.out ~/labs/lab8-$MYGIT/project/project.tree.rps.pdf
```
We used this command to examine the distribution of Pfam domains across the proteins.
```
cut -f 1 ~/labs/lab8-$MYGIT/project/project.rps-blast.out | sort | uniq -c
cut -f 6 ~/labs/lab8-$MYGIT/project/project.rps-blast.out | sort | uniq -c
```
_______________________________________________________________________________________________________________________________________________________________________________
