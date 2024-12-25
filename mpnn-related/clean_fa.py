#!/home/amax/.conda/envs/mlfold/bin/python
from Bio import SeqIO
import sys, os

def filter_sequence(seq):
    return "".join([base for base in seq if base in "ACGT"])

def fasta_name_modify(input_file):
    with open(input_file,'r') as file:
        lines = file.readlines()
    with open(input_file,'w') as file:
        for line in lines:
            if line.startswith('>'):
                line_blank_remove = line.replace(" ","")
                line_name_modify = line_blank_remove.replace(",","_")
                line = line_name_modify
            file.write(line)

def remove_duplicates_and_save(input_fasta, output_unique_fasta, output_duplicates_fasta):
    fasta_name_modify(input_fasta)
    
    sequences = SeqIO.parse(input_fasta, "fasta")
    seen_sequences = set()
    unique_sequences = []
    duplicate_sequences = []

    for seq_record in sequences:
        seq_record.description = ""
        #seq_record.seq = filter_sequence(seq_record.seq)
        seq_str = str(seq_record.seq)
        if seq_str not in seen_sequences:
            seen_sequences.add(seq_str)
            unique_sequences.append(f">{seq_record.id}\n{seq_str}")
        else:
            duplicate_sequences.append(f">{seq_record.id}\n{seq_str}")

    with open(output_unique_fasta, "w") as unique_fasta_file:
        unique_fasta_file.write("\n".join(unique_sequences))

    with open(output_duplicates_fasta, "w") as duplicates_fasta_file:
        duplicates_fasta_file.write("\n".join(duplicate_sequences))

    print(f"Removed duplicates. The unique sequences are saved in {output_unique_fasta}")
    print(f"Duplicate sequences are saved in {output_duplicates_fasta}")

input_fasta = sys.argv[1]
output_unique_fasta = str(input_fasta).split('.')[0] + "_cleaned.fasta"
output_duplicates_fasta = str(input_fasta).split('.')[0] + "_duplicates_output.fasta"

remove_duplicates_and_save(input_fasta, output_unique_fasta, output_duplicates_fasta)