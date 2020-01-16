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
        print(line)
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
    motif = content[6][:9]
    cat = content[7]
    return sample_control(chrom, pos, ref, cat, seq)

def sample_control(chrom, pos, ref, cat, seq, window=150, bp=4):
    sites1 = []
    sites2 = []
    
    while(len(sites1) == 0 and len(sites2) == 0):
        lseg_lb = max((pos-1-window-bp), 0)
        lseg_ub = pos - bp - 1
        useg_lb = pos + bp
        useg_ub = upBound = min(len(seq), pos + window + bp)
        subseq1 = seq[lseg_lb:lseg_ub]
        subseq2 = seq[useg_lb:useg_ub]
        subseq1 = re.sub(r'^N+', '', subseq1)
        subseq2 = re.sub(r'N+$', '', subseq2)

        sites1 = [m.start() for m in re.finditer(ref, subseq1)]
        sites2 = [m.start() for m in re.finditer(ref, subseq2)]
        sites1 = [s for s in sites1 if (s >= bp and s < (len(subseq1)-bp))]
        sites2 = [s for s in sites2 if (s >= bp and s < (len(subseq2)-bp))]
        window += 50 #expand window in edge case where mut_site is only ref_allele in window
    window -= 50    
    flip = random.randint(0, 1)
    if ((flip == 0 and len(sites1)>0) or (len(sites2)==0)):
        subseq = subseq1
        sites = sites1
    else:
        subseq = subseq2
        sites = sites2
    if(len(sites)==0):
        print("Bad pos: {}".format(pos))
    ix = random.choice(sites)
    newSeq = subseq[(ix - bp):(ix+bp+1)]
    while not re.search("[ATCG]{9}", newSeq):
        #print(pos)
        sites.remove(ix)
        ix = random.choice(sites)
        newSeq = subseq[(ix - bp):(ix+bp+1)]
    entry = {
        'chrom' : chrom,
        'pos' : pos,
        'motif' : newSeq,
        'cat': cat,
        'ref': ref,
        'window': window
    }
    return entry



if __name__ == "__main__":
    random.seed( 8675 ) # threeeeee ohhhhh niiiiii-eee-iiiiine
    main()
