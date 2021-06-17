# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-03-09
# @Description: Parses the Interpro, pannzer2 and the manual annotation file into a single annotation file 
###Attention: This only works if the primary transcript per gene as annotated 
##Gene without a DE are marked as hypothetical proteins

#general libraries
import argparse, csv

#functions from personal script 
import parse_interpro, parse_pannzer2_output

#libraries that need instalation
from goatools import obo_parser


def manual_annotation_parser(manual_file): 
    """Parses the manual annotation """
    manual_map = {}
    manual_annotation = csv.reader(open(manual_file, "r"), delimiter = "\t")
    for row in manual_annotation: 
        manual_map[row[1]] = f"{row[0]}*"
    return manual_map

def manual_go_terms_parser(manual_go): 
    """Parses information about the manual annotated genes GO terms. 
    The transcriptome gene id is in the first row, the annotation in the third. GO terms are seperated by "; \" """
    manual_go_map = {}
    manual_annotation = csv.reader(open(manual_go, "r"), delimiter = "\t")
    for row in manual_annotation: 
        manual_go_map[row[1]] = [term for term in row[2].split("; ")]
    return manual_go_map

def output_final_file(manual, manual_go, pannzer, interpro, go, g2t):
    """Ouput the final annotation file (named annotation_output.tsv)"""
    t2g = {} #saves a dictionary with transcript as a key and gene as a items 
    #I could replace this with a simples read from a list of genes but that would require changes in the rest of the code AND I might want to use the g2t eventually
    g2thandler = csv.reader(open(g2t, "r"), delimiter = "\t")
    for grow in g2thandler: 
        if grow[1].endswith(".1") and grow[2] == "mRNA": #Since I only annoatated the primary transcript per locus,of ncRNA genes, which always ends with .1 I only need to parse one these transcripts
            t2g[grow[1]] = grow[0]
    
    output_row_list = [["Gene id", "DE", "GOBP", "GOBP_terms", "GOMF", "GOMF_terms", "GOCC", "GOCC_terms", "CATH-Gene3D", "CDD", "PANTHER", "PFAM", "PIRSF", "PRINTS", "Prosite Patterns", "Prosite profiles", "SMART", "SFLD", "SUPERFAMILY", "TIGRFAM", "KEGG", "EC"]] #list of output rows. Each transcript in the t2g dictionary is processed one at a time
    #each output row has the following collums
    #Gene id DE (or manual annotation) GOBP_codes, GOBP_terms, GOMF_codes, GOMF_terms, GOCC_codes, GOCC_terms, CATH-Gene3D, CDD, PANTHER, PFAM, PIRSF, PRINTS, Prosite Patterns, Prosite profiles, SMART, SFLD, SUPERFAMILY, TIGERFAMS, KEGG, EC

    for transcript, gene in t2g.items(): 
        current = [""] * 22
        current[0] = gene

        #Parse Pannzer annotation
        if transcript in pannzer.keys():
            if "Description" in pannzer[transcript].keys(): #Parse decription annotation from pannzer
                current[1] = pannzer[transcript]["Description"][0] #0 is necessary as there is only one DE descriptor but it is still saved as a list 
            else: 
                current[1] = "hypothetical protein"
            
            if "go" in pannzer[transcript].keys(): ##Parse GO annotation from pannzer2
                for go_annotation in pannzer[transcript]["go"]: 
                    go_term_name = go[go_annotation].name

                    #for BP
                    if go[go_annotation].namespace == "biological_process": 
                        if current[2] != "": 
                            current[2] += ";"
                        current[2] += go_annotation
                        if current[3] != "": 
                            current[3] += ";"
                        current[3] += go_term_name
                    
                    #for MF
                    if go[go_annotation].namespace == "molecular_function": 
                        if current[4] != "": 
                            current[4] += ";"
                        current[4] += go_annotation
                        if current[5] != "": 
                            current[5] += ";"
                        current[5] += go_term_name

                    if go[go_annotation].namespace == "cellular_component": 
                        if current[6] != "": 
                            current[6] += ";"
                        current[6] += go_annotation
                        if current[7] != "": 
                            current[7] += ";"
                        current[7] += go_term_name
        
            #Parse KEGG and EC annotation
            if "EC" in pannzer[transcript].keys(): 
                for ec_annotation in pannzer[transcript]["EC"]: 
                    if current[21] != "": 
                        current[21] += ";"
                    current[21] += ec_annotation

            if "KEGG" in pannzer[transcript].keys(): 
                for ec_annotation in pannzer[transcript]["KEGG"]: 
                    if current[20] != "": 
                        current[20] += ";"
                    current[20] += ec_annotation
        
        ###Parse interpro annotation (for each gene)
        if transcript in interpro.keys(): 
            if "Gene3D" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["Gene3D"]: 
                    if current[8] != "": 
                        current[8] += ";"
                    current[8] += annotation

            if "CDD" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["CDD"]: 
                    if current[9] != "": 
                        current[9] += ";"
                    current[9] += annotation

            if "PANTHER" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["PANTHER"]: 
                    if current[10] != "": 
                        current[10] += ";"
                    current[10] += annotation

            if "Pfam" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["Pfam"]: 
                    if current[11] != "": 
                        current[11] += ";"
                    current[11] += annotation

            if "PIRSF" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["PIRSF"]: 
                    if current[12] != "": 
                        current[12] += ";"
                    current[12] += annotation

            if "PRINTS" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["PRINTS"]: 
                    if current[13] != "": 
                        current[13] += ";"
                    current[13] += annotation

            if "ProSitePatterns" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["ProSitePatterns"]: 
                    if current[14] != "": 
                        current[14] += ";"
                    current[14] += annotation

            if "ProSiteProfiles" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["ProSiteProfiles"]: 
                    if current[15] != "": 
                        current[15] += ";"
                    current[15] += annotation

            if "SMART" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["SMART"]: 
                    if current[16] != "": 
                        current[16] += ";"
                    current[16] += annotation

            if "SFLD" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["SFLD"]: 
                    if current[17] != "": 
                        current[17] += ";"
                    current[17] += annotation

            if "SUPERFAMILY" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["SUPERFAMILY"]: 
                    if current[18] != "": 
                        current[18] += ";"
                    current[18] += annotation

            if "TIGRFAM" in interpro[transcript].keys(): 
                for annotation in interpro[transcript]["TIGRFAM"]: 
                    if current[19] != "": 
                        current[19] += ";"
                    current[19] += annotation


        #replace the DE annotation if there is a manual annotation
        if gene in manual.keys():
            current[1] = manual[gene] 

        #supplement the GO annotation of the manual annotation
        if gene in manual_go.keys(): 
            #parse the manual annotation
            for term in manual_go[gene]: 
                go_term_name = go[term].name
                if go[term].namespace == "biological_process": 
                    if term not in current[2].split(";"): 
                        if current[2] != "": 
                            current[2] += ";"
                        current[2] += term

                        if current[3] != "": 
                            current[3] += ";"
                        current[3] += go_term_name

                if go[term].namespace == "molecular_function": 
                    if term not in current[4].split(";"): 
                        if current[4] != "": 
                            current[4] += ";"
                        current[4] += term

                        if current[5] != "": 
                            current[5] += ";"
                        current[5] += go_term_name

                if go[term].namespace == "cellular_component": 
                    if term not in current[6].split(";"): 
                        if current[6] != "": 
                            current[6] += ";"
                        current[6] += term

                        if current[7] != "": 
                            current[7] += ";"
                        current[7] += go_term_name

        if current[1] == "": #if there is no description append hypothetical protein
            current[1] = "hypothetical protein"

        #finally append the current list to the final list 
        output_row_list.append(current)
    
    #output the file 
    with open("annotation_output.tsv", "w") as fofh:
        for row in output_row_list: 
            for i in range(len(row)): 
                if row[i] == "": 
                    row[i] = "Na"
            if len(row) != 22: #for debug  
                print(row)
            fofh.write("\t".join(row))
            fofh.write("\n")


if __name__ == "__main__": 
    ##Parse arguments
    parser = argparse.ArgumentParser(description= "Parse the Intepro annotation file")
    parser.add_argument("-m", "--input-manual", help = "Input tsv file of the manual annotation (where the first collum is the gene name and second colum is the gene id")
    parser.add_argument("-p", "--input-pannzer", help = "Input tsv file of the pannzer input")
    parser.add_argument("-i", "--input-interpro", help = "Input tsv file of the interpro input")
    parser.add_argument("-go", "--input-go", help = "Input go obo file")
    parser.add_argument("-g2t", "--transcript2gene", help = "transcript to gene list of the transcriptome")
    parser.add_argument("-mgo", "--manual-go", help = "Input manual GO term annotation file")
    args = parser.parse_args()

    manual_map = manual_annotation_parser(args.input_manual)
    manual_go_map = manual_go_terms_parser(args.manual_go)

    interpro_map = parse_interpro.parse_Interpro_annotation(args.input_interpro)
    pannzer_map = parse_pannzer2_output.parse_pannzer_anotation(args.input_pannzer)
    go_tree = obo_parser.GODag(args.input_go)

    output_final_file(manual_map, manual_go_map, pannzer_map, interpro_map, go_tree, args.transcript2gene)