
from typing import Tuple
from singlecellmultiomics.utils import hamming_distance

def k_v_splitter(kv: str) -> Tuple[str, str]:
    parts = kv.split(':')
    return parts[0], parts[-1]



def add_to_readset(path, readset, sequences_per_cell):
    
    
    name = path.split('/')[-1].split('.')[0].replace('dCarlnday','').replace('dCarlndy','')
    day_n = name[0]
    if day_n=="5" or day_n=="7":
        suffix = " abc"[int(name[2])]
    else:
        suffix = ""
        
    r_list_index = f'day{day_n}{suffix}'
    r_cell_suffix = f'_d{day_n}{suffix}'

    library = r_list_index
    with open(path) as h:
        for line in h:
            cell, read, n = line.rstrip().split('\t')
            cell = f'{cell}{r_cell_suffix}'
            if 'N' in read:
                continue
            if len(read)<20:
                continue
            readset[read] += int(n)
            sequences_per_cell[library][cell][read] += int(n)
            
            

def trim_last_match_from_cigar(cigartuples):
    if len(cigartuples)==1 and  cigartuples[-1][1] == 'M': # If all matches, just return 'WT'
        return [(1,'WT'),]
    
    if cigartuples[-1][1] == 'M':
        return cigartuples[:-1]
    return cigartuples

def tuples_to_cigar_str(cigartuples):
    return ''.join([str(x[0])+x[1] for x in cigartuples])

def get_alignment(read, ref, sw, ref_name):
    return read, trim_last_match_from_cigar( sw.align(ref, read, query_name="nvt", ref_name=ref_name, q_qual=None).cigar )

