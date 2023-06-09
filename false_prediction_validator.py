import argparse
import re
import gffutils
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Validates false predictions flagged by liftoff")

parser.add_argument("-a", "--assembly", type=str, required=True, help="Input path for genome assembly in fasta")
parser.add_argument("-p", "--prediction", type=str, required=True, help="Input path for predictions in gff3")

args = parser.parse_args()

# parse gff3 file and store in a sqlite3 database
db = gffutils.create_db(args.prediction, "test.db", force=True, merge_strategy="merge")

# print the database, write to gff3
for g in db.features_of_type('gene', order_by=('seqid', 'start')):
    print()
    print(g)
    for f in db.children(g, order_by='start'):
        print(f)
        