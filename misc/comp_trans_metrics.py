"""Simple script to compute transcriptome metrics

Computes: 

- N50 
- Median transcript len 
- Average transcript len
- Number of sequences 
- Number of Genes (optional) 
- Busco analysis (optional) 

"""
import os
import argparse
import re 
import numpy as np
import pandas as pd 

from Bio import SeqIO



def run_busco_analysis(fasta_file, args): 
    """Runs Busco analysis"""
    #busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]
    busco_comand = f"busco -i {fasta_file} -l {args.linage} -o {fasta_file.split('.fa')[0]} -m transcriptome".split(" ")
    os.subprocess(busco_comand) 

    def parse_busco_file(short_summary_file) -> dict: 
        """Parse Busco resuls"""
        content = open(short_summary_file)
        for line in content:
            if "Complete and single-copy BUSCOs" in line:
                comp = int(line.split("\t")[1])
            elif "Complete and duplicated BUSCOs" in line:
                dupl = int(line.split("\t")[1])
            elif "Fragmented BUSCOs" in line:
                frag = int(line.split("\t")[1])
            elif "Missing BUSCOs" in line:
                miss = int(line.split("\t")[1])
        return {'Complete': comp, 'Duplicate': dupl, 'Fragmented': frag, 'Missing': miss}


    short_summary_file = f"{fasta_file.split('.fa')[0]}_short_summary.specific{args.lineage}.txt"
    if not short_summary_file in os.listdir(): 
        print(f"No BUSCO results found for {fasta_file}")
        return

    data = parse_busco_file(short_summary_file)
    return data


def compute_n50(unsorted_array): 
    """Compute N50"""
    sorted_array = np.sort(unsorted_array)
    csum = np.cumsum(sorted_array)
    half = int(np.sum(sorted_array)/2)
    csumn2 = min(csum[csum >= half])
    ind = np.where(csum == csumn2)
    n50 = sorted_array[ind[0]]
    return n50
    
def compute_basic_metrics(sequences, args) -> list:
    """Compute basic metrics for a fasta file
    N50 
    Median transcript len 
    Average transcript len
    Number of transcripts 
    Number of genes
    """
    len_list = [] 
    gene_counts = set()
    gene_regex = re.compile(f"({args.geneNamePattern})") 
    
    for seq in sequences: 
        len_list.append(len(seq.seq))
        if args.geneNamePattern != None:
            result = gene_regex.match(f"{seq.id}").groups()
            gene_counts.add(result[0])
    len_array = np.array(len_list)
    mean = np.mean(len_array)
    median = np.median(len_array)
    number_of_transcipts = len(len_list)
    if args.geneNamePattern != None:
        number_of_genes = len(gene_counts)
    else: 
        number_of_genes = number_of_transcipts
    n50 = compute_n50(len_array)
    return [mean, median, number_of_genes, number_of_transcipts, n50[0]]
    

def get_mode(fasta_file) -> str: 
    """Check if script should be ran on a folder or a file"""
    if os.path.isfile(fasta_file): 
        return 'single_file'
    return 'folder'


def main(args):
    """Main function""" 
    #Determine run mode and file list
    if get_mode(args.f) == 'folder': 
        files_in_folder = os.listdir(args.f)
        files = [os.path.join(args.f, file) for file in files_in_folder if (file.endswith(".fa") or file.endswith(".fasta"))]
    else: 
        files = [args.f]

    #Loop through files and run anaylsis
    row_names = []
    result_map = {'Mean': [], 
                'Median': [], 
                'Nº genes': [], 
                'Nº transcripts': [], 
                'N50': []}
    
    if args.busco:
        result_map['Busco Complete Singe-Copy'] = []
        result_map['Busco Complete Duplicates'] = []
        result_map['Busco Fragmented'] = []
        result_map['Busco Missing'] = []

    for fasta_file in files: 
        sequences = SeqIO.parse(fasta_file, "fasta")
        mean, median, number_of_genes, number_of_transcipts, n50 = compute_basic_metrics(sequences, args)
        result_map['Mean'].append(mean)
        result_map['Median'].append(median)
        result_map['Nº genes'].append(number_of_genes)
        result_map['Nº transcripts'].append(number_of_transcipts)
        result_map['N50'].append(n50)
        if args.busco: 
            data = run_busco_analysis(fasta_file, args) #Quick and dirty implementation
            result_map['Busco Complete Singe-Copy'].append(data['Complete'])
            result_map['Busco Complete Duplicates'].append(data['Duplicate'])
            result_map['Busco Fragmented'].append(data['Fragmented'])
            result_map['Busco Missing'].append(data['Missing'])
        row_names.append(fasta_file.split('.fa')[0])
    result_dataframe = pd.DataFrame(data=result_map)
    result_dataframe.index = row_names

    result_dataframe.to_csv(args.out, sep='\t')

def argument_parser(): 
    parser = argparse.ArgumentParser(description = "Compute transcriptome metrics (N50, median transcript len, avg transcript len, number of sequences, number of Genes, Busco analysis (optional)")
    parser.add_argument("-f", help = "Input (either specify a fasta file or a folder with fasta files (ending in either .fa or .fasta))", dest = "f", type = str, required=True)
    parser.add_argument("-out", help = "Output file with transcriptome metrics", dest = "out", type = str, default = "out.tsv") 
    parser.add_argument('-b', '--busco', help='Run busco analysis (optional, time consuming)', default = False)
    parser.add_argument('-l', '--lineage', help='Specify lineage to use in BUSCO', default = False)
    parser.add_argument("-geneNamePattern", help = 'Gene Panttern, provided in Regex. Used to count genes the number of genes in a file (regex). E.g TRINITY_DN1000_c115_g5 -> TRINITY_DN\d+_c\d+_g\d+', default=None)
    args = parser.parse_args()
    main(args) 

if __name__ == '__main__': 
    argument_parser() 
