# Misc 
This folder will be home to all the script that do not fall under a specified type of analysis step

### **filter_sequences.py**

Filter sequences based on size. Filtered out sequences are printed in the terminal

Example usage: 

``` python3 filter_sequences.py -s minSize -fasta_file InputFastaFile -output OutputFastaFile ``` 


### **replace_non_atgcn_nucl**

Replace all the non ATGCN nucleotides in a file for 'n'

Example usage: 

```  python3 replace_non_atgcn_nucl.py -f inputFastaFile -o outputFastaFile ``` 

### **salmon_parser.py** 

Parses salmon_quant.log files. Recursively look for subfolders with salmon_quant.log under the root directory

Example usage: 

```  python3 salmon_parser.py -f rootDirectory -out outputTable ``` 