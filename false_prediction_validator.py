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

for gene in db.features_of_type("gene"):
    # this needs input validation: what if valid_ORF is null?
    valid_ORF = gene.attributes["valid_ORF"]  # for some reason, this returns a list with a single item, hence index 0

    if valid_ORF[0] == "False":
        missing_start_codon = gene.attributes["missing_start_codon"]

        if missing_start_codon[0] == "True":
            search_start = None
            search_end = None

            if gene.strand == "+":
                search_start = int(gene.end)
                search_end = search_start + 15
            elif gene.strand == "-":
                search_end = int(gene.start)
                search_start = search_end - 15

            print(f"{search_start}, {search_end}")
            
# print the database, write to gff3
# for g in db.features_of_type('gene', order_by=('seqid', 'start')):
#     print()
#     print(g)
#     for f in db.children(g, order_by='start'):
#         print(f)
