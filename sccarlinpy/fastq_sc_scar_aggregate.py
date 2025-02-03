#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple
from more_itertools import chunked
from singlecellmultiomics.utils import hamming_distance
import pandas as pd
from collections import Counter, defaultdict
import argparse
import numpy as np
import gzip
import argparse
from .scar_functions import k_v_splitter
import pysam

def primer_filter_and_deduplicate(fastq_path_or_bam: str, 
                                out_path: str, 
                                primer: str, 
                                max_primer_hamming_distance:int = 2,
                                minlength:int=None,
                                target_contig=None,
                                ) -> Tuple[ np.array, defaultdict ]:
    """
    Read fastq file and demultiplex and deduplicate, returns a histogram of hamming distances to the specified primer,
    and a Counter of cell to sequences
    """
    def split_tag(s):
        tag, tagtype,value = s.split(':')
        return tag, value
        

    def iter_records(bam_or_fastq, target_contig):
        dropped_wrong_location = 0
        if bam_or_fastq.endswith('.bam'):
            bam = pysam.AlignmentFile(bam_or_fastq, 'rb')
            for read in bam:
                if read.is_secondary or read.is_supplementary:
                    continue
                if target_contig is not None and read.is_mapped and read.reference_name != target_contig:
                    dropped_wrong_location+=1
                    continue

                t = read.get_tag('RX')
                # "RX:Z:TGCTCCGGAAGT CB:Z:ACGTACAGTGACTCGC MI:Z:ACGTACAGTGACTCGCTGCTCCGGAAGT"
                parts = t.split(' ')
                tags = dict( (split_tag(part) for part in parts[1:]))
                tags['RX'] = parts[0]
                yield read.query_name, read.query_sequence, None, read.query_qualities, tags
            print(f'Excluded {dropped_wrong_location} reads mapping at unexpected location')
        else:
            with gzip.open(bam_or_fastq,'rt') as h:
                for i, (query_name, query, comment, query_qual ) in enumerate(chunked(h, 4)):
                    tags = dict( (k_v_splitter(kv)  for kv in query_name.split()[1:]))
                    yield query_name, query, comment, query_qual, tags


    cell_to_sequences = defaultdict(Counter) # Cell to (carlin array) sequences
    h_dists = pd.Series(np.zeros(len(primer)+1), index=range(len(primer)+1)) # Histogram of Hamming distances
    for query_name, query, comment, query_qual, tags in iter_records(fastq_path_or_bam, target_contig):
        prefix = query[:len(primer)]
        h_dist = hamming_distance(primer, prefix)
        h_dists[h_dist] += 1
        if h_dist <= max_primer_hamming_distance:
            
            seq = primer + query[len(primer):-1] # minus one to strip the newline
            if minlength is not None and len(seq)<minlength:
                continue
            cell_to_sequences[tags['CB']][seq] += 1
            
    with open(out_path, 'w') as out:
        for cell, cell_data in cell_to_sequences.items():
            for seq, count in cell_data.most_common():
                out.write(f'{cell}\t{seq}\t{count}\n')
                
    return h_dists, cell_to_sequences




def main_cli():

    argparser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Identify reads starting with the expected carlin primer, demultiplex and count unique sequence occurrences')

    argparser.add_argument('fastq', type=str, help='Can also use BAM to exclude reads which map to the genome')
    argparser.add_argument('-o',type=str)
    argparser.add_argument('-primer',type=str, default='CGGATTAACGTGTAAGCGGC')
    argparser.add_argument('-max_primer_hamming_distance',type=int, default=2, help='Maximum hamming distance between primer and beginning of read')
    argparser.add_argument('-minlength',type=int, default=None, help='Minimum read length')
    argparser.add_argument('-target_contig',type=str, help='Keep reads which map to this contig, or which do not map at all.')
    
    args = argparser.parse_args()

    primer_filter_and_deduplicate(args.fastq, args.o, primer=args.primer, max_primer_hamming_distance=args.max_primer_hamming_distance, minlength=args.minlength, target_contig=args.target_contig)