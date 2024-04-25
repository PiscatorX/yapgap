#!/usr/bin/env python
import itertools as it
import pandas as pd
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
               header_tuples.append(line.strip().split("\t"))
           if (n_line == max(header_lines)):
               break

        headers_raw = [' '.join(pair).strip() for pair in  pair_zip(header_tuples) ]
        self.header_map =  [ self.head_map.get(key, key) for key in  headers_raw]

      
    def combine_columns(self, row):
   
        attributes_cols = [ col for col in self.trna_df.columns if col not in self.gff_headers ]  

        return dict([(col, row[col]) for col in attributes_cols ])
        
      
    def extract_feat(self, f_obj, header_map):

        self.trna_df = pd.read_csv(self.f_obj, delimiter="\t", names = header_map, encoding="ascii")
        self.trna_df = self.trna_df[self.trna_df.seqid != "--------"]
        
        #seqid  tRNA #  start  end Anti Type  Codon  SeqID  SeqLen  score
        self.trna_df = self.trna_df.assign(source = 'tRNAscan-SE', type = 'tRNA', strand = 1, phase = '.')
        self.trna_df.loc[self.trna_df.loc[:,'start'] > self.trna_df.loc[:,'end'], 'strand']  = -1
        self.trna_df['attributes'] = self.trna_df.apply(self.combine_columns, axis=1)
        


    def write_features(self):

       self.gff_records = []
       for index, row in self.trna_df.iterrows():
           seq = Seq('', length = None)
           rec = SeqRecord(seq, row.pop('seqid'))
           row.attributes.update({ 'source': self.source })
           qualifiers = row.attributes
           start = int(row['start'])
           end = int(row['end'])
           if (start > end):
               start, end = (end, start)
           top_feature = SeqFeature(FeatureLocation(start, end),
                                    type = self.type,
                                    strand = row.strand,
                                    qualifiers=qualifiers)
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
    
    


# from BCBio import GFF
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.SeqFeature import SeqFeature, FeatureLocation

# # ##gff-version 3
# # ##sequence-region ID1 1 20
# # ID1     prediction      gene    1       20      10.0    +       .       ID=gene1;other=Some,annotations
# # ID1     prediction      exon    1       5       .       +       .       Parent=gene1
# # ID1     prediction      exon    16      20      .       +       .       Parent=gene1

# out_file = "your_file.gff"
# seq = Seq("GATCGATCGATCGATCGATC")
# rec = SeqRecord(seq, "ID1")
# qualifiers = {
#     "source": "prediction",
#     "score": 10.0,
#     "other": ["Some", "annotations"],
#     "ID": "gene1",
# }

# sub_qualifiers = {"source": "prediction"}
# top_feature = SeqFeature(
#     FeatureLocation(0, 20), type="gene", strand=1, qualifiers=qualifiers
# )
# top_feature.sub_features = [
#     SeqFeature(FeatureLocation(0, 5), type="exon", strand=1, qualifiers=sub_qualifiers),
#     SeqFeature(
#         FeatureLocation(15, 20), type="exon", strand=1, qualifiers=sub_qualifiers
#     ),
# ]
# rec.features = [top_feature]



