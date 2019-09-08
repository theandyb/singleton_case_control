from pyfaidx import Fasta
import re
import random
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Sample control distribution.")
    parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
    parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
    parser.add_argument("-o", "--output", help="Path to output", required=True)
    parser.add_argument("chrom", help="Chromosome we are sampling from")
    args = parser.parse_args()

    ref_file = args.fasta #"reference_data/human_g1k_v37/human_g1k_v37.fasta"
    singleton_file = args.singleton #"/net/snowwhite/home/beckandy/research/smaug-redux/summaries/singletons.full.summary"
    chrom = args.chrom
    output_list = []
    # Create fasta object
    fasta_obj = Fasta(ref_file, read_ahead=10000, as_raw=True)
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    with open(singleton_file) as fp:
        line = fp.readline()
        line = fp.readline()
        while line:
            output_list.append(process_line(line, chrom, fasta_obj))
            line = fp.readline()
    final = pd.DataFrame(output_list)
    final.to_csv(args.output, index = None, header=True)
    print("Done!")


def process_line(x, chrom, fsObj):
    content = x.split("\t")
    pos = int(content[1])
    ref = content[2]
    motif = content[3][:9]
    cat = content[4]
    return sample_control(chrom, pos, ref, cat, fsObj)

def sample_control(chrom, pos, ref, cat, fsObj, window=300, bp=4):
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
        'cat': cat,
        'mutation' : 0
    }
    return entry



if __name__ == "__main__":
    main()
