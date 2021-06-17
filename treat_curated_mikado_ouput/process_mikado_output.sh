#!/bin/bash 

# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-03-27
# @Description: Bash script to process the MIKADO curated output. This script take as input the curated mikado file 
#(transcriptome_assembly.complete.gff3) as a reference instead of the usual mikado output

## Outputs: (review)
## 1) Gff3 file with just the primary transcript per locus
## 2) Gff3 file with all the transcipts (including isoforms) for coding genes + nc genes
## 3) Gff3 file with all the transcipts (including isoforms) for coding genes (which will be used in quantification) 
## 4) Gff3 file with all the transcipts (including isoforms) for nc genes
## 5) Gff3 compressed files and index for genome browser 
## 6) Fasta files (.mRNA, .cds .pep)
## 7) Gene to transcript table (with the additional collum with the classication into ncRNA and coding RNA)
## 8) List of class codes of ncRNA compared with the coding regions
## 9) Output statistics related to the transcriptome (N50 and distributions)
## 10) GTF file to be used in the featureCounts quantification

transcriptome_assembly.complete.gff=$1
genome.fasta=$2

#Ouput assembly stats 
mikado util stats transcriptome_assembly.complete.gff3 transcriptome_assembly.complete.stats.tsv

#run the treat_mikado_output.py to split mikado output into several files
python3 scripts/treat_mikado_out.py -i transcriptome_assembly.complete.gff3 -o transcriptome_assembly

#sort the output file and make the files for jbrowser
#complete file 
gt gff3 -sortlines -tidy -retainids transcriptome_assembly.complete.gff3 > transcriptome_assembly.complete.sorted.gff3
bgzip transcriptome_assembly.complete.sorted.gff3
tabix transcriptome_assembly.complete.sorted.gff3.gz

#mRNA file 
gt gff3 -sortlines -tidy -retainids transcriptome_assembly.mRNA.gff3 > transcriptome_assembly.mRNA.sorted.gff3
bgzip transcriptome_assembly.mRNA.sorted.gff3
tabix transcriptome_assembly.mRNA.sorted.gff3.gz

#convert to fasta files transcripts #-y options to output the proteins, -x to output the cds and -w to output the mRNA
gffread transcriptome_assembly.complete.gff3 -g genome.fasta -y transcriptome_assembly.complete.pep.fasta
gffread transcriptome_assembly.mRNA.gff3 -g genome.fasta -y transcriptome_assembly.mRNA.pep.fasta
gffread transcriptome_assembly.primary.gff3 -g genome.fasta -y transcriptome_assembly.primary.pep.fasta

gffread transcriptome_assembly.complete.gff3 -g genome.fasta -x transcriptome_assembly.complete.cds.fasta
gffread transcriptome_assembly.mRNA.gff3 -g genome.fasta -x transcriptome_assembly.mRNA.cds.fasta
gffread transcriptome_assembly.primary.gff3 -g genome.fasta -x transcriptome_assembly.primary.cds.fasta

gffread transcriptome_assembly.complete.gff3 -g genome.fasta -w transcriptome_assembly.complete.transcripts.fasta
gffread transcriptome_assembly.mRNA.gff3 -g genome.fasta -w transcriptome_assembly.mRNA.transcripts.fasta
gffread transcriptome_assembly.primary.gff3 -g genome.fasta -w transcriptome_assembly.primary.transcripts.fasta

#tranform the mRNA gff3 output into gtf (this file will be used for quantification) 
gffread transcriptome_assembly.mRNA.gff3 -T > transcriptome_assembly.mRNA.gtf

#get a gene to transcript tsv file (with a third collum which either show mRNA or ncRNA)
python3 scripts/gene2transcript_table.py -i transcriptome_assembly.complete.gff3 -o g2t.tsv

#Compute N50 for the complete and coding sequences
sh scripts/N50.sh transcriptome_assembly.complete.transcripts.fasta > n50_results.txt
sh scripts/N50.sh transcriptome_assembly.mRNA.transcripts.fasta >> n50_results.txt

## Compare the ncRNA with the mRNA to get the corresponding class (-T option is to ensure there is no rmap and tmap files created). This will give us an idea of the ncRNA classes we have 
gffcompare -r transcriptome_assembly.mRNA.gff3 -T transcriptome_assembly.ncRNA.gff3 -o gffcompare_ncRNA

#use greps to find the number of each class code in the file 
for code in = c k m n j e o s x i y p u;  
do  
	echo -e "$code\t$(grep "class_code \"$code\"" gffcompare_ncRNA.annotated.gtf | wc -l)" >> ncRNA_class_codes_stats.tsv;
done