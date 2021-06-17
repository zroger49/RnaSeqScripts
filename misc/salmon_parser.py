# @Author: Rogerio Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-03-19
# @Description: Replace non ATGCN nucleotides for n (this makes it compatible with dammit)


import argparse
import os

def get_log_files_path(input_folder):
    file_map = {} 
    for sample in os.listdir(input_folder): 
        file_list = os.walk("{}/{}".format(input_folder, sample))
        for path,folder,files in file_list: 
            for file_log in files: 
                if file_log == "salmon_quant.log": 
                    file_map[sample] ="{}/{}".format(path, file_log).replace("/", "\\")
    return file_map


def parse(input_folder): 
    data_map = {}
    file_map = get_log_files_path(input_folder)
    for sample, file_path in file_map.items(): 
        with open(file_path, "r") as fh: 
            content = fh.readline()
            while not content.startswith("Total"): 
                content = fh.readline()
            mapped_reads = content.split(": ")[1].replace("\n", "")
            content = fh.readline()
            unique_mapped = content.split(": ")[1].replace("\n", "")
            content = fh.readline()
            ambi_mapped = content.split(": ")[1].replace("\n", "")
            data_map[sample] = [mapped_reads, unique_mapped, ambi_mapped]
    return data_map

def output_table(data_map, output_file): 
    with open(output_file, "w")as ofh:
        ofh.write("\t".join(["sample", "Aligned", "unique", "ambigous"]))
        ofh.write("\n")
        for key, values in data_map.items():
            ofh.write("\t".join([key, values[0], values[1], values[2]]))
            ofh.write("\n")



def argument_parser(): 
    parser = argparse.ArgumentParser(description="Parse salmon")
    parser.add_argument("-f", help = "folder containing subfolders (one subfolder per sample)", dest = "f", type = str, required=True)
    parser.add_argument("-out", help = "Output table", dest = "out", type = str, default = "out.tsv") 
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args) 

def main(args): 
    data_map = parse(args.f)
    output_table(data_map,args.out)

if __name__ == "__main__":
    argument_parser()