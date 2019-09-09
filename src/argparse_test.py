import argparse
import os.path
from os import path

parser = argparse.ArgumentParser(description='Sample control distribution.')
parser.add_argument("-s", "--singleton", help="Location of singleton file", required=True)
parser.add_argument("-f", "--fasta", help="FASTA file to grab sequence from", required=True)
args = parser.parse_args()
print("{} - {}".format(args.singleton, args.fasta))