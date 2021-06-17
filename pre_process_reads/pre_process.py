# @Author: Rogério Ribeiro
# @E-mail: zroger499@gmail.com
# @Date:   2020-01-29
# @Description: RNA-seq data pre-processing script

# TODO
# Replace hardcorded dirs e cleanup for an argparser/or just use Snakemake

import os
import subprocess

CLEANUP = True #clean intermediary files 

# raw data dir
raw_data_dir = "00_raw_data"
trimomatic_path = "/home/up201505534/software/Trimmomatic-0.39"

##Different sequencing libraries were processed in different ways
smart_seq = ["idio_1","idio_2", "idio_3", "meso_1", "meso_2", "meso_3", "leaf_1", "leaf_2", "totalptt_1", "totalptt_2", "totalptt_3"]
true_seq = ["P1F1Int","P2F1Int","P3F1Int","P1F1Ext","P2F1Ext","P3F1Ext","P1F4Int","P2F4Int","P3F4Int","P1F4Ext","P2F4Ext","P3F4Ext",]
nano = ["P1F2Int", "P2F2Int", "P3F2Int", "P1F2Ext", "P2F2Ext", "P3F2Ext", "P1F3Int", "P2F3Int", "P3F3Int", "P1F3Ext", "P2F3Ext","P3F3Ext",]

# Make directories for fastqc and processed reads
#Fastqc with raw data
os.makedirs("01_fastqc/illumina")
os.makedirs("01_fastqc/nanopore")

#Processed samples
os.makedirs("01_processed/illumina")
os.makedirs("01_processed/nanopore/reports")
os.makedirs("01_processed/nanopore/pychopper")


#Fastqc of the processed data
os.makedirs("02_fastqc_processed/illumina")
os.makedirs("02_fastqc_processed/nanopore")


for sample in os.listdir(raw_data_dir):
    file_list = os.walk("{}/{}".format(raw_data_dir, sample))
    for path, folder, files in file_list:
        
        #Code for processing Illumina samples
        if sample in smart_seq or sample in true_seq:
            # Set the R1 and the R2 file
            if "R1" in files[0]:
                R1_file = os.path.join(raw_data_dir, sample, files[0])
                R2_file = os.path.join(raw_data_dir, sample, files[1])
            else:
                R1_file = os.path.join(raw_data_dir, sample, files[1])
                R2_file = os.path.join(raw_data_dir, sample, files[0])

            # Fastqc comand:
            # fastqc R1_file R2_file -o 01_fastqc -f fastq
            fastqc_call = ["fastqc", R1_file, R2_file, "-o", "01_fastqc/illumina", "-f", "fastq"]
            subprocess.call(fastqc_call)

            ##trimmomatic comand : java -jar $trimomaticPath/trimmomatic-0.39.jar PE -threads 10 -phred33 R1_raw.fastq.gz R2_raw.fastq.gz
            ##R1_cut_paired.fastq.gz  R1_cut_unpaired.fastq.gz  R2_cut_paired.fastq.gz  R2_cut_unpaired.fastq.gz
            ##ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
            
            #For smart_seq also add ILLUMINACLIP:NexteraPE-PE.fa:2:30:1
            trimmomatic_out = os.path.join("01_processed/illumina", sample)
            os.mkdir(trimmomatic_out)

            if sample in smart_seq:
                trimmomatic_call = ["java","-jar",os.path.join(trimomatic_path, "trimmomatic-0.39.jar"), 
                "PE", "-threads", "10", "-phred33", 
                R1_file, R2_file, 
                os.path.join(trimmomatic_out, "{}_R1_cut_paired.fastq.gz".format(sample)),
                os.path.join(trimmomatic_out, "{}_R1_cut_unpaired.fastq.gz".format(sample)),
                os.path.join(trimmomatic_out, "{}_R2_cut_paired.fastq.gz".format(sample)),
                os.path.join(trimmomatic_out, "{}_R2_cut_unpaired.fastq.gz".format(sample)),
                "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
                "ILLUMINACLIP:NexteraPE-PE.fa:2:30:10",
                "LEADING:3",
                "TRAILING:10",
                "SLIDINGWINDOW:4:15",
                "MINLEN:50"]
            
            else:
                trimmomatic_call = ["java","-jar",os.path.join(trimomatic_path, "trimmomatic-0.39.jar"), 
                "PE", "-threads", "10", "-phred33", 
                R1_file, R2_file, 
                os.path.join(trimmomatic_out, "{}_R1_cut_paired.fastq.gz".format(sample)),
                os.path.join(trimmomatic_out, "{}_R1_cut_unpaired.fastq.gz".format(sample)),
                os.path.join(trimmomatic_out, "{}_R2_cut_paired.fastq.gz".format(sample)),
                os.path.join(trimmomatic_out, "{}_R2_cut_unpaired.fastq.gz".format(sample)),
                "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
                "LEADING:3",
                "TRAILING:10",
                "SLIDINGWINDOW:4:15",
                "MINLEN:50"]
            
            subprocess.call(trimmomatic_call)

            
            # Fastqc comand 2
            # fastqc R1_cut_paired_file R2_cut_paired_file -o 02_fastqc_processed -f fastq
            fastqc_call = ["fastqc", os.path.join(trimmomatic_out, "{}_R1_cut_paired.fastq.gz".format(sample)), 
                os.path.join(trimmomatic_out, "{}_R2_cut_paired.fastq.gz".format(sample)),
                "-o", "02_fastqc_processed/illumina",
                "-f","fastq"]

            subprocess.call(fastqc_call)

        #Code for processing nanopore samples
        if sample in nano: 
            read_file = os.path.join(raw_data_dir, sample, files[0])
            
            # Fastqc comand:
            # fastqc nano_sample.fastq -o fastqc/ -f fastq
            fastqc_call = ["fastqc", read_file, "-o", "01_fastqc/nanopore", "-f", "fastq"]
            subprocess.call(fastqc_call)

            #Orient_nanopore reads
            #cdna_classifier.py P1F3Ext.fastq P1F3Ext.oriented.fastq -r P1F3Ext_report.pdf 
            #-w P1F3Ext_rescued_reads.fastq -S P1F3Ext_stat.tsv -u P1F3Ext_unclass_reads.fastq 
            #-z 50 -l P1F3Ext.small.fastq
            pychopper_intermediary_folder = "01_processed/nanopore/pychopper"
            pychopper_report_folder = "01_processed/nanopore/reports" 

            pychopper_call = ["cdna_classifier.py",  read_file, 
            os.path.join(pychopper_intermediary_folder, "{}.oriented.fastq".format(sample)), #output oriented reads
            "-r", os.path.join(pychopper_report_folder, "{}_report.pdf".format(sample)),  #output pdf report 
            "-w", os.path.join(pychopper_intermediary_folder, "{}.rescued.reads.fastq".format(sample)),  #output rescued reads (reads with internal adaptors)
            "-S", os.path.join(pychopper_report_folder, "{}_stat.tsv".format(sample)), #output stat file
            "-u", os.path.join(pychopper_intermediary_folder,"{}_unclass_reads.fastq".format(sample)), #output non-oriented reads 
            "-z", "50", "-l", os.path.join(pychopper_intermediary_folder,"{}.small.fastq".format(sample))] #output small reads
            
            subprocess.call(pychopper_call)

            #rescue reads might give valuable information. Concatenate the oriented and rescue reads: 
            #cat reads.oriented reads.rescued > reads.oriented.all
            output_reads_file = open(os.path.join("01_processed/nanopore", "{}.oriented.all.fastq".format(sample)), "w")
            concatenate_call = ["cat", os.path.join(pychopper_intermediary_folder, "{}.oriented.fastq".format(sample)), 
            os.path.join(pychopper_intermediary_folder, "{}.rescued.reads.fastq".format(sample))]

            subprocess.call(concatenate_call, stdout=output_reads_file)
            output_reads_file.close()

            
            #call fastqc on the concatenated files: 
            #fastqc reads.oriented.all.fastq -o 02_fastqc_processed/nanopore
            fastqc_call = ["fastqc", os.path.join("01_processed/nanopore", "{}.oriented.all.fastq".format(sample)), 
            "-o", "02_fastqc_processed/nanopore"]
            
            subprocess.call(fastqc_call)


if CLEANUP == True:
    subprocess.call(["rm", "-r" "01_processed/nanopore/pychopper"]) #remove the intermediaries pychopper files
