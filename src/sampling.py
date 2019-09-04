from pyfaidx import Fasta
import re

def main():
    ref_file = "reference_data/human_g1k_v37/human_g1k_v37.fasta"
    singleton_file = "testData.tsv"

    fasta_obj = Fasta(ref_file, read_ahead=10000, as_raw=True)


def process_line(x, fsObj):
    content = x.strip().split("\t")
    chrom = content[0]
    pos = content[1]
    ref = content[2]

def sample_control(chrom, pos, ref, fsObj, window=300):
    seq = fsObj['{}'.format(chrom)][(pos-1-window):(pos-1+window)]


if __name__ == "__main__":
    main()