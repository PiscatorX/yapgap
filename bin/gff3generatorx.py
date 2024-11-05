#!/usr/bin/env python
import itertools as it
import csv
import copy
import pprint
import argparse
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation



class GenGFF(object):
   
    def __init__(self):
   
        self.gff_headers = ['seqid' , 'source' , 'type' , 'start' , 'end', 'score', 'strand', 'phase', 'attributes']


    def ParseFile(self):

        self.f_obj = open(self.input)
        self.tRNA_fpass_out()
        self.extract_feat(self.f_obj, self.header_map)


        
    def tRNA_fpass_out(self, header_lines = [0, 1]):

        self.head_map = {'Sequence Name':'seqid',
                  'tRNAscan-SE':'source',
                  'tRNA Bounds Begin':'start',
                  'tRNA End':'end',
                  'Score':'score'}

        pair_zip = lambda l: it.zip_longest(l[0],l[1], fillvalue = '')
        header_tuples = []

        for n_line, line in  enumerate(self.f_obj):
           if n_line in header_lines:
               header = line.strip().split("\t")
               if line.startswith("Sequence"):
                   pass
               header_tuples.append(header)
           if (n_line == max(header_lines)):
               break

        headers_raw = [' '.join(pair).strip() for pair in  pair_zip(header_tuples) ]
        
        self.header_map =  [ self.head_map.get(key, key) for key in  headers_raw]      

        print(self.header_map)
        
      
    def extract_feat(self, f_obj, header_map):

       
        self.trna_df = []
        
        reader = csv.reader(self.f_obj, delimiter='\t')

        for row in reader:
            if row[0].startswith('-'):
                continue
           
            row_dict = dict(zip(header_map, row))
            qualifiers = dict((k,v ) for k,v in row_dict.items() if k not in self.gff_headers)
            row_dict = dict((k.lower(), str(v).lower() ) for k,v in row_dict.items() if k in self.gff_headers)
            row_dict.update({'source': 'tRNAscan-SE', 'type': 'tRNA', 'phase': '.'})
            row_dict['strand'] = '-1' if int(row_dict['start'])  >  int(row_dict['end']) else '1'
            row_dict['qualifiers'] = qualifiers
            self.trna_df.append(row_dict)
            

    def write_features(self):

       self.gff_records = []
       for row in self.trna_df:
           seq = Seq('', length = None)
           rec = SeqRecord(seq, row.pop('seqid'))
           if not self.source:
               row.update({ 'source': self.source })
               
           start = int(row['start'])
           end = int(row['end'])
           if (start > end):
               start, end = (end, start)
               
           top_feature = SeqFeature(FeatureLocation(start, end),
                                    type = self.type,
                                    qualifiers = row['qualifiers'] )
           rec.features = [top_feature]
           
           self.gff_records.append(rec)
           
       with open(self.outfile, "w") as out_handle:
           GFF.write(self.gff_records, out_handle)



    def cmd_arguments(self):
      
        parser = argparse.ArgumentParser(description="Parse an output file and generate a GFF3 file")
        parser.add_argument("input", help="Path to the file to be parsed")
        parser.add_argument('-s','--source', action='store', required=False, default='tRNAscan',
                            help="type of file to converted ")
        parser.add_argument('-t','--type',  action='store', required=False, default= 'SO:0001410',
                            help="GFF3 entry for 'type'")
        parser.add_argument('-o','--outfile',  action='store', required=True, help="GFF3 output file")


        args = parser.parse_args()
        self.__dict__.update(args.__dict__)
        


if __name__ == "__main__":
   
    gff3_writer = GenGFF()
    gff3_writer.cmd_arguments()
    gff3_writer.ParseFile()
    gff3_writer.write_features()
