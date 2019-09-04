from pyfaidx import Fasta
import re
import random
import pandas as pd

def main():
    ref_file = "reference_data/human_g1k_v37/human_g1k_v37.fasta"
    singleton_file = "testData.tsv"
    output_list = []
    # Create fasta object
    fasta_obj = Fasta(ref_file, read_ahead=10000, as_raw=True)
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    with open(singleton_file) as fp:
        line = fp.readline()
        line = fp.readline()
        while line:
            output_list.extend(process_line(line, fasta_obj))
            line = fp.readline()
    final = pd.DataFrame(output_list)
    final.to_csv("data/sample.csv", index = None, header=True)
    print("Done!")


def process_line(x, fsObj):
    new_pair = []
    content = x.split("\t")
    chrom = content[0]
    pos = int(content[1])
    ref = content[2]
    motif = content[6][:9]
    new_pair.append({
        'chrom':chrom,
        'pos':pos,
        'motif': motif,
        'mutation' : 1
    })
    new_pair.append(sample_control(chrom, pos, ref, fsObj))
    return new_pair

def sample_control(chrom, pos, ref, fsObj, window=300, bp=4):
    seq = fsObj['{}'.format(chrom)][(pos-1-window):(pos-1+window)]
    sites = [m.start() for m in re.finditer(ref, seq)]
    ix = random.choice(sites) + pos - 1 - window
    while ix == pos:
        ix = random.choice(sites) + pos - 1 - window
    newSeq = fsObj['{}'.format(chrom)][(ix - bp):(ix+1+bp)]
    while 'N' in newSeq:
        ix = random.choice(sites) + pos - 1 - window
        newSeq = fsObj['{}'.format(chrom)][(ix - bp):(ix+1+bp)]
    entry = {
        'chrom' : chrom,
        'pos' : ix,
        'motif' : newSeq,
        'mutation' : 0 
    }
    return entry



if __name__ == "__main__":
    main()