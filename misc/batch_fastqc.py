# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-01-30
# @Description: Run fastqc on the folder


import os
import subprocess


if __name__ == "__main__":
    output_dir = "01_fastqc"
    os.mkdir(output_dir)
    os.mkdir(os.path.join(output_dir, "illumina"))
    os.mkdir(os.path.join(output_dir, "nanopore"))
    home_dir = os.getcwd()

    for sample in os.listdir(home_dir):
        file_list = os.walk("{}/{}".format(home_dir, sample))
        for path, folder, files in file_list:
            print(files)
            if len(files) == 1: 
                R1_file = os.path.join(home_dir, sample, files[0])
                fastqc_call = ["fastqc", R1_file, "-o", os.path.join(output_dir, "nanopore"), "-f", "fastq"]
                subprocess.call(fastqc_call)
				
            elif len(files) == 2: 
                R1_file = os.path.join(home_dir, sample, files[0])
                R2_file = os.path.join(home_dir, sample, files[1])

                fastqc_call = ["fastqc", R1_file, R2_file, "-o", os.path.join(output_dir, "illumina"), "-f", "fastq"]
                subprocess.call(fastqc_call)