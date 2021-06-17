# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-03-07
# @Description: Parse the Pannzer2 output in tsv format. The output is a dictionary of dictionaries
#also outputs stats 

import argparse,csv
from goatools import obo_parser


def parse_pannzer_anotation(input_file): 
    """Parse the Pannzer2 input. Results is a dictionary of dictionaries"""
    main_map = {}
    pannzer2_annotation = csv.reader(open(input_file, "r"), delimiter = "\t")
    next(pannzer2_annotation) #skip the header
    for row in pannzer2_annotation: 
        if row[1] == "DE" or "BP_" in row[1] or "MF_" in row[1] or "CC_" in row[1] or "EC_" in row[1] or "KEGG_" in row[1]: #ignore the other rows
            if row[0] not in main_map.keys(): #check if the protein is in the map
                main_map[row[0]] = {}
            if row[1] == "DE": #parse decription line
                information = row[5]
                key = "Description"
            if "BP" in row[1] or "MF" in row[1] or "CC" in row[1]: #parses the GO lines. I am assuming there can not be repeated terms
                information = f"GO:{row[4]}"
                key = "go"
            
            if "EC_" in row[1]: #Parses the EC line along with the GO term
                information = row[4]
                key = "EC"
            
            if "KEGG_" in row[1]: #parses the Kegg line along with the GO term
                information = row[4]
                key = "KEGG"

            #append the information to the dictionary
            if key in main_map[row[0]].keys(): #check if there is already information about this type of annotation
                main_map[row[0]][key].append(information)
            else: 
                main_map[row[0]][key] = [information]

            if key == "KEGG" or key == "EC": ###These line also contain information about Gene Ontology
                information_2 = row[5]
                if "go" in main_map[row[0]].keys(): #check if there is already information about this type of annotation. 
                    if information_2 not in main_map[row[0]]["go"]: #Also need to check if the GO terms is already in the list 
                        main_map[row[0]]["go"].append(information_2)
                else: 
                    main_map[row[0]]["go"] = [information_2]
                

    return main_map
        
def compute_statistics(main_map, go): 
    """Return a map with the statistics of the annotation"""
    statistics_map = {"Annotated": 0, "DE": 0, "KEGG":0, "EC": 0, "GO":0, "GO_BP": 0, "GO_MF": 0, "GO_CC":0, "Total_GO_BP":0,"Total_GO_MF":0, "Total_GO_CC":0}
    for value in main_map.values(): 
        statistics_map["Annotated"] += 1
        for annotation_type, annotations in value.items():
            if annotation_type == "Description": 
                statistics_map["DE"] += 1
            if annotation_type == "KEGG": 
                statistics_map["KEGG"] += 1
            if annotation_type == "EC": 
                statistics_map["EC"] += 1
            if annotation_type == "go":
                statistics_map["GO"] += 1 
                #set flags to false. If one terms from one of these ontology is found then add to the final output
                flag_bp = False 
                flag_mf = False
                flag_cc = False 
                for go_term in annotations:
                    ontology = go[go_term].namespace
                    if ontology == "biological_process": 
                        statistics_map["Total_GO_BP"] += 1 
                        flag_bp = True
                    if ontology == "molecular_function": 
                        statistics_map["Total_GO_MF"] += 1 
                        flag_mf = True
                    if ontology == "cellular_component": 
                        statistics_map["Total_GO_CC"] += 1 
                        flag_cc = True
                if flag_bp == True: 
                    statistics_map["GO_BP"] += 1 
                if flag_mf == True: 
                    statistics_map["GO_MF"] += 1 
                if flag_cc == True: 
                    statistics_map["GO_CC"] += 1 
    
    #Compute the average number of BP terms per annotated gene with any type of GO terms (total GO_BP/GO)
    statistics_map["Average_BP"] = statistics_map["Total_GO_BP"] / statistics_map["GO"]
    statistics_map["Average_MF"] = statistics_map["Total_GO_MF"] / statistics_map["GO"]
    statistics_map["Average_CC"] = statistics_map["Total_GO_CC"] / statistics_map["GO"]
    
    #return the dictionary
    return statistics_map



if __name__ == "__main__": 
    ##Parse arguments
    parser = argparse.ArgumentParser(description= "Parse the Intepro annotation file")
    parser.add_argument("-i", "--input-pannzer", help = "Input tsv file")
    parser.add_argument("-go", "--input-go", help = "Input go obo file")
    args = parser.parse_args()

    #output the dictionary
    main_map = parse_pannzer_anotation(args.input_pannzer)
    
    #compute statistics of the annotation 
    #go = obo_parser.GODag(args.input_go)
    #statistics = compute_statistics(main_map, go)
    #print(statistics)