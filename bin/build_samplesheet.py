#!/usr/bin/env python
from argparse import RawTextHelpFormatter
import os
import sys
import glob
import pprint
import argparse

def build_samplesheet(read_dir, spliton, samplesheet, paired):

    search_dir = os.path.join(read_dir,"*.f*q*")
    read_files  =  list(map(os.path.abspath, sorted(glob.glob(search_dir))))

    if not read_files:raise Exception("No reads found check read directory path. Do the read extensions match the pattern '*.f*q*' ")
    
    label  = lambda fname:  os.path.basename(fname.rsplit(spliton, 1)[0])
    
    n = len(read_files) 

    samplesheet_fp = open(samplesheet, 'w')
    
    if paired:
        read_data = {}
        for i in range(0, n, 2):
            paired_reads = read_files[i:i+2]
            sample_IDs = list(map(label,  paired_reads))
            if (sample_IDs[0]) != sample_IDs[1]:
                raise Exception("Reads appear to be unpaired based on sample IDs")
            read_data[list(set(sample_IDs))[0]] = list(map(os.path.basename, paired_reads))
            fmt = "{},{},{}"
            
        print("sample,fastq_1,fastq_2", file = samplesheet_fp)
        
    else:
        read_data = dict((label(seq_read), os.path.basename(seq_read)) for  seq_read in  read_files)
        fmt = "{},{}"
        print("sample,fastq_1", file = samplesheet_fp)

     
        
    for k,v in read_data.items():

        if paired:
            print(fmt.format(k,v[0], v[1]), file = samplesheet_fp)
        else:
            print(fmt.format(k,v), file = samplesheet_fp)
            
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Generates a samplesheetfile for paired-end data https://docs.qiime2.org/2021.8/tutorials/importing/#fastq-samplesheet-formats \nAssumes that reads have have the '*.fastq' | *.fq | '*.fastqz' extension.", formatter_class=RawTextHelpFormatter)

    parser.add_argument('-r','--read-dir', dest='read_dir',
                        action='store', required=True, type=str)
    parser.add_argument('-s','--spliton', dest='spliton', default="_", help ="For sample name detection. The character in all read filename immediately before readnumber. e.g for 'RemovePrimer_Final.LB6_1.fq.gz' spliton = '_'",
                        action='store', required=False, type=str)
    parser.add_argument('-f','--samplesheet', dest='samplesheet',
                        action='store', required=True, type=str)
    parser.add_argument('-p','--paired', default = False,action='store_true')
    
    args = parser.parse_args()
    
    build_samplesheet(args.read_dir, args.spliton, args.samplesheet, args.paired)
