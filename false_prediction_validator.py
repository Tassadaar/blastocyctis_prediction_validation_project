import argparse
import re
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Validates false predictions flagged by liftoff")

parser.add_argument("-a", "--assembly", type=str, required=True, help="Input genome assembly in fasta")
parser.add_argument("-p", "--prediction", type=str, required=True, help="Input predictions in gff3")

args = parser.parse_args()

