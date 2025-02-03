#!/usr/bin/env python
import argparse
import gzip
from more_itertools import chunked

def demux_and_trim(r1_path, r2_path, out_path, umilen=12):

    with gzip.open(r1_path,'rt') as r1,\
        gzip.open(r2_path,'rt') as r2,\
        gzip.open(out_path,'wt',compresslevel=2) as out:
        
        for read1, read2 in zip(chunked(r1,4), chunked(r2,4)):
            cb, umi = read1[1][:16],read1[1][16:16+umilen].rstrip()
            header = read2[0].rstrip().split()[0]
            
            sequence = read2[1].rstrip()
            comment = read2[2].rstrip()
            qualities = read2[3].rstrip()
            
            # Remove poly A if present
            poly_a_pos = sequence.find('AAAAAAAAAAAA')
            if poly_a_pos>0:
                sequence = sequence[:poly_a_pos]
                qualities = qualities[:poly_a_pos]
                
            if len(sequence)>8:
                out.write(
                    f'{header} RX:Z:{umi} CB:Z:{cb} MI:Z:{cb+umi}\n{sequence}\n{comment}\n{qualities}\n' 
                )

def main_cli():
    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Parse fastq files, demultiplex and write read 2 reads to disk""")
    argparser.add_argument(
        'r1',
        type=str,
        help='Read 1')
    argparser.add_argument(
        'r2',
        type=str,
        help='Read 2')
    argparser.add_argument(
        'output',
        type=str,
        help='Output file')
    args = argparser.parse_args()

    demux_and_trim(args.r1, args.r2, args.output)
    
