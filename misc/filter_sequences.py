# @Author: Rogerio Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2021-04-14
# @Description: Filter sequences based on their size 


import argparse, csv
from Bio import SeqIO


def filter(size, fasta_file, output): 
    """"Filter sequences based on min size. Print sequences that were excluded""" 
    #parse fasta file and get interesting sequences
    fasta_results = []
    fasta = SeqIO.parse(fasta_file, "fasta")
    for sequence in fasta: 
        if len(sequence.seq) >= size: 
            fasta_results.append(sequence)
        else: 
            print("{}, len = {}".format(sequence.id, len(sequence.seq)))

    SeqIO.write(fasta_results, output, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description= "Filter sequences based on size limit")
    parser.add_argument("-s", "--size", help = "min size")
    parser.add_argument("-fasta_file", "--fasta", help = "Input Gene list")
    parser.add_argument("-output", "--out", help = "Output file with filtered sequences")
    args = parser.parse_args()

    filter(int(args.size), args.fasta, args.out) 
