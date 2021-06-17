# @Author: Rogerio Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-02-05
# @Description: Replace non ATGCN nucleotides for n (this makes it compatible with dammit)



from Bio import SeqIO
from Bio.Seq import Seq 
import argparse


def run(input_file,output_file): 
    fasta_parser = SeqIO.parse(input_file, "fasta")
    new_seq_list = []
    for SeqObj in fasta_parser: 
        new_seq = "" 
        for nucl in SeqObj.seq: 
            if nucl.upper() in ["A", "T", "C", "G", "N"]: 
                new_seq += nucl 
            else: 
                new_seq += "N"
        SeqObj.seq = Seq(new_seq)
        new_seq_list.append(SeqObj)
    SeqIO.write(new_seq_list, output_file, "fasta")


if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description= "Test the transcriptome for fused ref transcript")
    parser.add_argument("-f", "--fasta", help = "Input transcriptome fasta file")
    parser.add_argument("-o", "--output", help = "Output fasta file")
    arg = parser.parse_args()

    run(arg.fasta, arg.output)


