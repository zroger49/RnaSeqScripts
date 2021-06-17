# @Author: Rogerio Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2021-03-12
# @Description: Extract sequences from fasta files. Ouput is a fasta filw with the name of the list.fasta

import argparse, csv
from Bio import SeqIO


def extract(gene_list, fasta_file): 
    """Main"""
    #Get the input genes
    with open(gene_list, "r") as ofh:
        my_gene_list = csv.reader(ofh)
        extract_gene_list = [gene[0] for gene in my_gene_list]
    #parse fasta file and get interesting sequences
    fasta_results = []
    fasta = SeqIO.parse(fasta_file, "fasta")
    for sequence in fasta: 
        gene_id = sequence.id.split(".")[0]
        if gene_id in extract_gene_list: 
            fasta_results.append(sequence)

    #output the sequence
    SeqIO.write(fasta_results, "output.fasta", "fasta")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= "Output a fasta file with the sequences requested")
    parser.add_argument("-l", "--gene-list", help = "Input Gene list")
    parser.add_argument("-fasta_file", "--fasta", help = "Input Gene list")
    args = parser.parse_args()

    extract(args.gene_list, args.fasta)