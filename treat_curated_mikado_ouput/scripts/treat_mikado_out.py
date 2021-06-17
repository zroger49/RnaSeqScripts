# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date: 2020-03-27
# @Description: Parse and process the output of mikado pick
# 1)Convert the gff3 file into: one file with just the primary transcript, one file with all the transcript and one file with just the coding transcript
# (removes Lines with Super Locus, as well as any atributes referencing it)
# 2) renames all the transcripts (from mikado.something.XGXXX to prefix_XXX)

import argparse, csv, re

def rename_atribute(row, rename_dictionary, prefix, name_counter): 
    """Rename the atributes in the Gff3 file, based on the ID atribute of the gff3 file. This function expects to find the "gene" row before the child atributes"""
    gene_id =  row[8].split(";")[0].split("=")[1].split(".")[1]
    if gene_id not in rename_dictionary.keys(): #new entry in the dictionary if the gene id has not been found yet
        name_counter += 1
        rename_dictionary[gene_id] = f"{prefix}_{name_counter}"


    row[8] = row[8].replace(f"mikado.{gene_id}", rename_dictionary[gene_id])
    return row, rename_dictionary, name_counter


def is_primary(row, primary_list):
    """Helper function to determine if the row being parsed belong to a primary or non primary transcript"""
    #primary status is returned as a string and not as a boolean 
    if row[2] == "ncRNA" or row[2] == "mRNA":  #Primary atribute are in these rows
        primary = re.search("primary=(True|False)", row[8]).group(1)
        if primary == "True":
            transcript_id = row[8].split(";")[0].replace("ID=", "") #Transcript name is in these rows ID field
            primary_list.append(transcript_id)
        return primary, primary_list
    else: #parse other rows
        transcript_id = row[8].split(";")[1].replace("Parent=", "")
        if transcript_id in primary_list:
            return "True", primary_list
        else:
            return "False", primary_list
        
def get_id(col9): 
    """Helper function to get the ID of a feature"""
    return col9.split(";")[0].replace("ID=", "")


def parse_gff3_file(input, output_prefix, header): 
    if header == True: 
        header_list = [] #Keep a list with the header lines
    primary_list = []#List for primary transcript
    mRNA_locus_list = [] #List with name of the mRNA loci 
    gff3_primary = [] #List with the lines relating to the primary transcripts in each locus 
    gff3_all = [] #List with all the transcript
    gff3_mRNA = [] #List with mRNA transcripts
    gff3_ncRNA = [] #List with the ncRNA transcripts
    with open(input, "r") as gff3fh:
        mikado_loci_gff3 = csv.reader(gff3fh, delimiter = "\t")
        for row in mikado_loci_gff3: 
            if row[0].startswith("#"):
                if header == True and row[0] != "###": 
                    header_list.append(row[0])
                else: #Ignores lines with "###"
                    continue 
            
            #elif row[2] == "superlocus": #Ignore lines with SuperLocus 
             #   continue
            
            elif row[2] == "gene" or row[2] == "ncRNA_gene": #remove the super locus information from the "genes" lines and append them only to the gff3_all list (it is important to have the gene line since there are multiple transcript per gene)
                gff3_all.append(row)
                if row[2] == "gene": #append mRNA to the respective list 
                    gff3_mRNA.append(row)
                    mRNA_locus_list.append(get_id(row[8])) #append the name of the locus to the list
                else: #implicit to be a row related to ncRNA 
                    gff3_ncRNA.append(row)

            else: 
                primary, primary_list = is_primary(row,primary_list) #returns if the row belong to a primary transcript and the updated primary_list
                if primary == "True": ##Only append rows from primary transcripts
                    gff3_primary.append(row)
                
                if ".".join(get_id(row[8]).split(".")[0:1]) in mRNA_locus_list: #check if the gene is in the mRNA gene list
                    gff3_mRNA.append(row)
                else: #implicit to be a ncRNA genes
                    gff3_ncRNA.append(row)
                
                gff3_all.append(row)
    
    #Ouput files
    with open(f"{output_prefix}.primary.gff3", "w") as pofh: 
        if header == True: 
            for hr in header_list: 
                pofh.write(hr)
                pofh.write("\n")
        
        for row in gff3_primary: 
            pofh.write("\t".join(row))
            pofh.write("\n")

    with open(f"{output_prefix}.mRNA.gff3", "w") as mofh:
        if header == True: 
            for hr in header_list: 
                mofh.write(hr)
                mofh.write("\n")
        
        for row in gff3_mRNA: 
            mofh.write("\t".join(row))
            mofh.write("\n")
    

    with open(f"{output_prefix}.ncRNA.gff3", "w") as ncofh:
        if header == True: 
            for hr in header_list: 
                ncofh.write(hr)
                ncofh.write("\n")
        
        for row in gff3_ncRNA: 
            ncofh.write("\t".join(row))
            ncofh.write("\n")

    
if __name__ == "__main__": 
    ##Parse arguments
    parser = argparse.ArgumentParser(description= "Parse Mikado gff3 output")
    parser.add_argument("-i", "--input", help = "Input gff3 file")
    parser.add_argument("-o", "--output", help = "prefix for the output files")
    parser.add_argument("--header", help = "Keep the file in the output header?", action = 'store_true') 
    args = parser.parse_args()

    if args.header: 
        header = True
    else: 
        header = False

    parse_gff3_file(args.input, args.output, header)