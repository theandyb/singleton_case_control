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
    parser.add_argument("-x", "--skip", help = "Number of additional lines to skip", required=False, type=int)
    parser.add_argument("chrom", help="Chromosome we are sampling from")
    args = parser.parse_args()

    ref_file = args.fasta #"reference_data/human_g1k_v37/human_g1k_v37.fasta"
    singleton_file = args.singleton #"/net/snowwhite/home/beckandy/research/smaug-redux/summaries/singletons.full.summary"
    chrom = args.chrom
    output_list = []
    # Create fasta object
    fasta_obj = Fasta(ref_file)
    seq = fasta_obj["{}".format(chrom)]
    seqstr = seq[0:len(seq)].seq
    print("FASTA read!")
    # Iterate over singletons file
    print("Sampling control observations for singletons...")
    counter = 1
    with open(singleton_file) as fp:
        next(fp)
        if args.skip:
            cLine = 1
            while cLine <= args.skip:
                next(fp)
                cLine += 1
        line = fp.readline()
        while line:
            output_list.append(process_line(line, chrom, seqstr))
            line = fp.readline()
            counter += 1
            if counter % 10000 == 0:
                print(counter)
                pd.DataFrame(output_list).to_csv(args.output, index = None, header=False, mode='a')
                output_list = []
    print("Done sampling...")
    if output_list:
        pd.DataFrame(output_list).to_csv(args.output, index = None, header=False, mode='a')
    print("Done!")


def process_line(x, chrom, seq):
    content = x.strip().split("\t")
    pos = int(content[1])
    ref = content[2]
    motif = content[3][:9]
    cat = content[4]
    return sample_control(chrom, pos, ref, cat, seq)

def sample_control(chrom, pos, ref, cat, seq, window=150, bp=4):
    lowBound = max((pos-1-window-bp), 1)
    upBound = min(len(seq), pos + window+bp)
    subseq = seq[lowBound:upBound]
    sites2 = [m.start() for m in re.finditer(ref, subseq)]
    sites = [s for s in sites2 if (s >= bp and s < (len(subseq)-bp))]
    if (window - 1 + bp) in sites:
        sites.remove(window - 1 + bp)
    ix = random.choice(sites)
    newSeq = subseq[(ix - bp - 1):(ix+bp)]
    while not re.search("[ATCG]{9}", newSeq):
        print("IT HAPPENED!")
        ix = random.choice(sites)
        newSeq = subseq[(ix - bp - 1):(ix+bp)]
    entry = {
        'chrom' : chrom,
        'pos' : pos,
        'motif' : newSeq,
        'cat': cat,
        'ref': ref
    }
    return entry



if __name__ == "__main__":
    main()
