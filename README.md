Gene loss among Eukarya that led to disjoint EFL1 ortholog function among species
===
In examining the possible gene losses that may have shifted homolog function within orthologs of the EFL1 gene family, a procedural series of bioinformatic techniques have been utilized to further understand its evolutionary history. By following this annotated repository, it is possible to replicate the findings of this study, as well as possibly apply its principles to other gene families of interest.

## Table of Contents

- [Gene loss among Eukarya that led to disjoint EFL1 ortholog function among species](#gene-loss-among-eukarya-that-led-to-disjoint-efl1-ortholog-function-among-species)
  * [Table of Contents](#table-of-contents)
  * [Introduction](#introduction)
  * [Methods and Procedure](#methods-and-procedure)
  * [Conclusion](#conclusion)

## Introduction

Among all three domains of life, protein synthesis is highly conserved and is directed by the complex protein-assembling machinery of ribosomes.

The EFL1, or the elongation factor-like GTPase 1 protein, is among one of the key players in ribosome synthesis, being inherently crucial in removing an anti-association factor eFI6 that leads to the full maturation of the large 60S subunit of the ribosomal body. EFL1 contains five structural domains, with domain II providing an insertion site for a variable length of amino acids that separates it from EF-2 translocases -- an important distinction for function.

![](https://i.imgur.com/zWanYoD.png)
(Asano et al. 2013)

It plays a major role in the etiology of SDS, or Schwann-Diamond syndrome, and is of interest to elucidate the phylogenetic relationship of EFL1 amongst orthologs within varying species, in order to deduce why some orthologs of EFL1 are species-specific at domain II within the protein’s architecture. In other words, did significant gene loss contribute to the loss of SBDS interaction among EFL1 orthologs?

In order to explore this question, the procedure below was used to uncover where in this gene's evolutionary history did gene loss among orthologs contribute to its functionality. 

Methods and Procedure
---
**1. Using BLAST to find EFL1 homologs**

First, create a database in which the BLAST search with be performed and download the NCBI FASTA formatted file of target gene ```XP_021362732.1```

```
mkdir /home/ec2-user/labs/lab3-L4-jle/efl1
pwd
ncbi-acc-download -F fasta -m protein XP_021362732.1
```
Then, run the BLAST search, making note of the top scoring hits. Use the command ```-less``` to peer into the file.

```
blastp -db …/allprotein.fas -query XP_021362732.1.fa -outfmt 0 -max_hsps 1 -out elf1.blastp.typical.out
less elf1.blastp.typical.out
```
It is important to then begin filtering the BLAST results from the query and view it as well.
```
blastp -db …/allprotein.fas -query XP_021362732.1.fa -outfmt “6 sseqid pident length mismatch gapopen evalue bitscore pident stitle” -max_hsps 1 -out efl1.blastp.detail.out
less -S efl1.blastp.detail.out
```
Next, it is important to begin filtering out only for high scoring, putative homologs. Be sure to use an e-value of 1e-35, using a program called "awk" that allows particular records in a file to be modified. Here, ```grep``` will search the query for all files named "efl1.blastp.detail.out"
```
grep -c Hsapiens efl1.blastp.detail.out
awk ‘{if ($6< 1e-35 )print $1 }’ efl1.blastp.detail.out > efl1.blastp.detail.filtered.out
wc -l efl1.blastp.detail.filtered.out
grep -o -E “[1][a-z]+.” efl1.blastp.detail.filtered.out | sort | uniq -c
```
The output of the BLAST search and filtering of the retrieved NCBI file should look as such:
```
A. planci 3
A. vaga 6
H. belcheri 3
C. intestinalis 4
D. melanogaster 3
E. granulosus 3
H. sapiens 3
L. anatina 4
M. yessoensis 5
```
Be sure to check if the display shows the number of similar homologs among the nine species explored in this study, as presented above.

**2. Using MUSCLE alignment to deduce sequence percent identity**

Using [seqkit](https://bioinf.shenwei.me/seqkit/), which is an ultrafast FASTA toolkit, to obtain the sequences in the selected FASTA protein file, and then view the output.

```
seqkit grep --pattern-file ~/labs/lab3-$MYGIT/efl1/efl1.blastp.detail.filtered.out ~/labs/lab3-$MYGIT/allprotein.fas > ~/labs/lab4-$MYGIT/efl1/efl1.homologs.fas
```
Following the retrieval of the target sequences, using [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/), align the sequences selected from the BLAST query:
``` 
muscle -in ~/labs/lab4-$MYGIT/efl1/efl1.homologs.fas -out ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas
```
Using [alv](https://github.com/arvestad/alv), a command-line alignment viewing tool, view the MUSCLE output and adjust the screen output using ```-majority``` 

```
alv -kli  ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas | less -RS
alv -kli --majority ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas | less -RS
```
Also, create the output as a pdf, named ```efl1.homologs.al.pdf``` for further viewing.
```
alv -ki -w 100 ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas | aha > ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.html
a2ps -r --columns=1 ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.html -o ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.ps
ps2pdf ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.ps ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.pdf
```
Using [alignbuddy](https://github.com/biologyguy/BuddySuite/blob/master/buddysuite/AlignBuddy.py), retrieve more information (width, length, gaps) regarding the aligned sequence:

```
alignbuddy  -al  ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas
alignbuddy -trm all  ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas | alignbuddy  -al
```
Then, remove and calculate the length of the alignment after removing all invariant positions in the sequence:

```
alignbuddy -dinv 'ambig' ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas | alignbuddy  -al
```

Another method of confirming the robustness of the alignment is to use [t_coffee](https://www.ebi.ac.uk/Tools/msa/tcoffee/) to calculate the percent identity.

```
t_coffee -other_pg seq_reformat -in ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas -output sim
```
It is also wise to compare that percent identity with another calculation from alignbuddy, in order to see which value sequence alignment is better aligned.

```
alignbuddy -pi ~/labs/lab4-$MYGIT/efl1/efl1.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```
**3. Constructing the Phylogenetic tree for EFL1 homologs Using Sequence Data**

[IQ-TREE](http://www.iqtree.org/) is a phylogenetic tree building algorithm that infers trees using maximum-likelihood. It will compile this data by choosing the most appropriate model, composed of three parts.

1) The amino-acid exchange rate matrix
2) The model-based frequencies, in this case the amino-acid frequency
3) The rate heterogeneity along places on the sequence, denoted by an integer specifying the number of categories


Using this program, the file ```efl1.homologs.al.fas``` will be used as the base data for this inference. The flag ```-bb 1000``` will be used as the base bootstrap support level for the tree. The most optimal score will be found here ```efl1.homologs.al.fas.iqtree```
```
iqtree -s ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.fas -bb 1000 -nt 2 
```
> Another interesting implication of using IQ-TREE is that is allows for the user to see rate homegeneity values, for which some sequences are invariant and have evolved very slowy, or others that have evolved quickly. 

To view it, use this command with the appropriate directory to the EFL1 homolog file:

```
nw_display ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.fas.treefile
```

Since the output from using IQ-TREE will be in the ```.iqtree``` format, the tree will be pasted in ASCII, to view in terminal. Therefore, to view fully, use [Rscript](https://www.r-project.org/) to export the .iqtree file as a viewable pdf. *Note: 0.4* in the code signifies the text size in the output pdf.
```
Rscript --vanilla ~/labs/lab5-$MYGIT/plotUnrooted.R  ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.fas.treefile ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.fas.treefile.pdf 0.4
```
**4. Midpoint root the optimal phylogenetic tree using gotree**

Oftentimes, it is sensible to root the earliest divergence point of the homologs of interest by using [gotree](https://github.com/evolbioinfo/gotree), a set a commandline tools written in [Go](https://go.dev/).

```
gotree reroot midpoint -i ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.fas.treefile -o ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.mid.treefile
```

Then, view it in the command line:
```
nw_order -c n ~/labs/lab5-$MYGIT/efl1/efl1.homologs.al.mid.treefile  | nw_display -
```
> Note: ```nw_order``` will order the clades by the descendent number, allowing the output tree to look cleaner and easily readable. 

**5. Using NOTUNG to reconcile the EFL1 gene family**

Sometimes, it is necessary to reconcile gene and species trees, due to the nature of many duplication and loss events that may have happened in its evolutionary history. [NOTUNG](http://www.cs.cmu.edu/~durand/Notung/) is a gene tree-species reconciliation software that estimates these duplication and loss events using a parsimony-based approach. 

The version that will be used is NOTUNG 3.0.

Using [java](https://www.java.com/en/download/help/whatis_java.html), install NOTUNG using this command line:

```
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar --help  
```
Then, use this following command to initialize NOTUNG and import the following target gene file ``` efl1.homologs.al.mid.treefile``` in order to reconcile:

```
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ~/labs/lab5-$MYGIT/species.tre -g ~/labs/lab6-$MYGIT/efl1/efl1.homologs.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/labs/lab6-$MYGIT/efl1/
```
View the output using this command, where ```-s``` cuts the long lines that are formed in the output:

```
less -s efl1.homologs.al.mid.treefile.reconciled.events.tx
```
Also, be sure to view and note any internal nodes that do not have formal taxonomic names using this command line:

```
grep NOTUNG-SPECIES-TREE ~/labs/lab6-$MYGIT/gqr/gqr.homologs.al.mid.treefile.reconciled | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
```
Another method to view the reconciled tree from NOTUNG is by viewing it in thirdkind.

In [thirdkind](https://academic.oup.com/bioinformatics/article/38/8/2350/6525213), by creating a RecPhyloXML object, a graphical network of the gene loss/duplication lineages will be mapped on top of each other, in an concise manner:

```
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/labs/lab6-$MYGIT/efl1/efl1.homologs.al.mid.treefile.reconciled --include.species
thirdkind -Iie -D 40 -f ~/labs/lab6-$MYGIT/gqr/gqr.homologs.al.mid.treefile.reconciled.xml -o  ~/labs/lab6-$MYGIT/gqr/gqr.homologs.al.mid.treefile.reconciled.svg
```

**6. Using RPS-BLAST to identify protein domains with Pfam**

[RPS-BLAST](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBWhat) is a specialized version of BLAST, in which it stands for "Reverse Position-Specific BLAST." Essentially, RPS-BLAST will use the query sequence to search backwards from the query to the subject, as opposed to BLAST which does the opposite. 

First, it is necessary to retrive the [Pfam ](https://www.ebi.ac.uk/interpro/)database, in which it is possible to conduct functional analysis of target proteins.

```
wget -O ~/data/Pfam_LE.tar.gz ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz && tar xfvz ~/data/Pfam_LE.tar.gz  -C ~/data
```

To run RPS-BLAST on the target protein ddomain, use this command. Note: The e-value is the width of the search of which certain domains will be returned. ```-evalue``` is set at 0.0000000001, which is appropriate for the EFL1 protein family; however, it may need to be adjusted for other families.

```
rpsblast -query ~/labs/lab8-$MYGIT/efl1/efl1.homologs.fas -db ~/data/Pfam -out ~/labs/lab8-$MYGIT/efl1/efl1.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```

In order to view the predicted Pfam domains of the phylogeny, it is necessary to incorporate two programs, [ggtree](https://guangchuangyu.github.io/software/ggtree/) and [drawProteins](https://www.bioconductor.org/packages/devel/bioc/vignettes/drawProteins/inst/doc/drawProteins_BiocStyle.html), in order to plot the pfam domains that was searched using RPS-BLAST.

The following command line has an element written in [R](https://www.r-project.org/) and involves several key elements.

1. ```sudo```: gives permission as computer administrator to install packages
2. ```--vanilla```: a commmand flag used to not restore a previous workspace environment
3. ```plotTreeAndDomains.r```: cited from Dr. Joshua Rest, PhD. It defines the boundaries of the output.

**7. Viewing the predicted Pfam domains at length**

Using this command, it is possible to view more detailed annotations of the target gene:

```
mlr --inidx --ifs "\t" --opprint  cat ~/labs/lab8-$MYGIT/efl1/efl1.rps-blast.out | tail -n +2 | less -S
```
In order to view the distribution of Pfam domains, use this command, that counts the domains automatically:

```
cut -f 1 ~/labs/lab8-$MYGIT/efl1/efl1.rps-blast.out | sort | uniq -c
```

If needed, in order to view the domain with the greatest number of residues, use this command line:

```
awk '{a=$4-$3;print $1,'\t',a;}' ~/labs/lab8-$MYGIT/gqr/gqr.rps-blast.out |  sort  -k2nr
```

Using the ```sort``` command, which will order the results depending on the user-set values, use to command line to find the domain with the most optimal e-value:

```
sort  -k5rg ~/labs/lab8-$MYGIT/efl1/efl1.rps-blast.out | less -S
```

## Conclusion
If followed thoroughly, the work flow outline is repeatable and may be applied to other gene families of interest. 

Some considerations to keep in mind:

1. This is a user-profile centered work flow; to repeat this procedure, be sure to appropriately create directories, as shown above, in a personal GitHub account.

2. In further viewing the Pfam domains, an R script written by Joshua Rest, PhD, was used. 

3. E-values for this experiment are only specific to the EFL1 protein family. Any other target family will require adjusted, if needed, in order to capture the scope of the study.


