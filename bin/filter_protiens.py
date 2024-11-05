#!/usr/bin/env python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import logging 


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

invalid_amino_acids = set("J")

def valid_seq(sequence):
    
    for aa in sequence:
        if aa in invalid_amino_acids:
            return (False,aa)
        
    return (True, None)


def clean_sequences(input_file, output):
    
    with open(output, "w") as out_handle:
        
        for record in SeqIO.parse(input_file, "fasta"):
            valid, residue =  valid_seq(record.seq)
            if valid:
                SeqIO.write(record, out_handle, "fasta")
            else:
                logging.info(f"Invalid sequence: {record.id}\tResidue: {residue}")
                

def main():
    
    parser = argparse.ArgumentParser(description="Filter and clean amino acid sequences in a FASTA file.")
    parser.add_argument('input_file', type=str, help='Path to the input FASTA file containing sequences.')
    parser.add_argument('-o','--output', type=str, required = True,  help='Path to the output FASTA file to save cleaned sequences.')
    args = parser.parse_args()
    clean_sequences(args.input_file, args.output)

if __name__ == "__main__":
    main()
