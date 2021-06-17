# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-03-04
# @Description: Parse the InterPro output in tsv format

import argparse,csv
from goatools import obo_parser

def compute_statistics(map, go_map):
    """Compute statistics for interpro annotation
    In this case, there o point computing the number of genes annotated with each db, since redudant 
    annotation are removed"""
    statistics_map = {"Annotated": 0, "Total_GO_BP":0,"Total_GO_MF":0, "Total_GO_CC":0}
    for value in main_map.values(): 
        statistics_map["Annotated"] += 1
        if "go" in value.keys():
            for go in value["go"]: 
                if go_map[go].namespace == "biological_process": 
                    statistics_map["Total_GO_BP"] += 1
                if go_map[go].namespace == "molecular_function": 
                    statistics_map["Total_GO_MF"] += 1
                if go_map[go].namespace == "cellular_component": 
                    statistics_map["Total_GO_CC"] += 1
            

    return statistics_map 
    

def parse_Interpro_annotation(input_file): 
    """Parse the InterPro annotation. Results is a dictionary of dictionaries. Keys are the sequence name and item is a dictionary where keys are the db and items are the description. Gene Ontology information is also stored"""
    main_map = {}
    inter_pro_annotation = csv.reader(open(input_file, "r"), delimiter = "\t")
    for row in inter_pro_annotation: 
        if row[0] not in main_map.keys(): #check if the protein is in the map
            main_map[row[0]] = {}
        if row[5] != "-": #sometimes there is no information on the domain
            information = f"{row[5]} ({row[4]})" #get the description as well as the acession number]
        else: 
            information = f"{row[4]}"
            
        if row[3] in main_map[row[0]].keys(): #check if there is an entry for that proteins and that analysis type
            if information not in main_map[row[0]][row[3]]: #check if the parsed information is not the dictionary (avoiding duplicate entries)
                main_map[row[0]][row[3]].append(information)
            #else skip adding the information 
        else: #if there is no information for that protein from this analysis
            main_map[row[0]][row[3]] = [information]
        
        ##parse the gene ontology information
        if len(row) >= 14:
            if row[13] != "-": #if there is gene ontology information
                go_list = row[13].split("|")
                for go in go_list:
                    if "go" in main_map[row[0]].keys() and go not in main_map[row[0]]["go"]:
                        main_map[row[0]]["go"].append(go)
                    else: 
                        main_map[row[0]]["go"] = [go]
            
    return main_map



if __name__ == "__main__": 
    ##Parse arguments
    parser = argparse.ArgumentParser(description= "Parse the Intepro annotation file")
    parser.add_argument("-i", "--input-interpro", help = "Input tsv file")
    parser.add_argument("-go", "--input-go", help = "Input go obo file")
    
    args = parser.parse_args()

    
    main_map = parse_Interpro_annotation(args.input_interpro)

    go = obo_parser.GODag(args.input_go)
    stats = compute_statistics(main_map, go)
    print(stats)