# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-02-17
# @Description: Parse a Gff file and get the transcript to gene table 

import argparse
import csv, re 

def get_gene2transcript_table(inpt,output): 
    """Main function. Takes as input the gff3 file and outputs a tsv file with gene*tab*transcript"""
    file_int = csv.reader(open(inpt, "r"), delimiter = "\t")
    transcript_list = [] #list of list, where the transcript are saved as [0]gene [1]transcript [2]mRNA/nCRNA
    for row in file_int: 
        if row[0].startswith("#"): 
            continue
        if row[2] == "mRNA" or row[2] == "ncRNA": 
            gene_id = re.search("Parent=(CATHA_\d+)", row[8]).group(1)
            transcript_id = re.search("ID=(CATHA_\d+\.\d+)", row[8]).group(1)
            transcript_list.append((gene_id, transcript_id, row[2]))
    
    with open(output, "w") as ofh: 
        for t in transcript_list:  
            ofh.write("\t".join(t))
            ofh.write("\n")
        


if __name__ == "__main__": 
    ##Parse arguments
    parser = argparse.ArgumentParser(description= "Read a gff3 file from Mikado and output a tsv table with gene\ttranscript")
    parser.add_argument("-i", "--input", help = "Input gff3 file")
    parser.add_argument("-o", "--output", help = "output file")
    args = parser.parse_args()

    get_gene2transcript_table(args.input, args.output)