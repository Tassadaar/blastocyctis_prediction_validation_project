import argparse
import re
import gffutils
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Validates false predictions flagged by liftoff")

    parser.add_argument("-f", "--fasta_file", type=str, required=True, help="Input path for genome assembly in fasta")
    parser.add_argument("-p", "--prediction", type=str, required=True, help="Input path for predictions in gff3")
    parser.add_argument("-o", "--output", type=str, required=True, help="Input name for the output file")

    args = parser.parse_args()

    # load fasta
    fasta = SeqIO.index(args.fasta_file, "fasta")
    # parse gff3 file and store in a sqlite3 database
    db = gffutils.create_db(args.prediction, dbfn=":memory:", id_spec=id_func)
    # counter for number of motifs found
    motif_counter = 0
    output = args.output

    with open(output, "w") as out:
        for gene in db.features_of_type("gene"):
            # validation
            if "valid_orf" and "missing_start_codon" not in gene.attributes:
                out.write(f"\n{str(gene)}\n")
                for child in db.children(gene, order_by="start"):
                    out.write(f"{str(child)}\n")
                continue

            valid_orf = gene.attributes["valid_ORF"]  # for some reason this returns a list with a single item
            missing_start_codon = gene.attributes["missing_start_codon"]

            if valid_orf[0] == "False" and missing_start_codon[0] == "True":

                if gene.strand == "+":
                    search_start = int(gene.end)
                    search_end = search_start + 15
                    motif_range = get_seq(fasta, gene.seqid, search_start, search_end)

                    if re.search(r"T[ACTG]{4}TGTTTGTT", motif_range):  # check if the motif exists
                        gene.attributes["motif"] = "True"
                        motif_counter += 1
                    else:
                        gene.attributes["motif"] = "False"

                elif gene.strand == "-":
                    search_end = int(gene.start)
                    search_start = search_end - 15
                    motif_range = get_seq(fasta, gene.seqid, search_start, search_end)

                    if re.search(r"AACAAACA[ACTG]{4}A", motif_range):  # check if the motif exists
                        gene.attributes["motif"] = ["True"]
                        motif_counter += 1
                    else:
                        gene.attributes["motif"] = ["False"]

            out.write(f"\n{str(gene)}\n")
            for child in db.children(gene, order_by="start"):
                out.write(f"{str(child)}\n")

    # status report
    prepend_to_file(output, f"# There were {motif_counter} transferred predictions with the motif found")


# written by Joran(novigit)
# this function creates a custom ID for each CDS feature, using the ID, the start and end coordinates to
# make sure it is unique
def id_func(x):
    new_id = None

    if x.featuretype == 'gene':
        new_id = x.attributes['ID'][0]
    elif x.featuretype == 'CDS':
        new_id = '-'.join([x.attributes['ID'][0], x.seqid, str(x.start), str(x.end)])

    return new_id


# this function extracts genomic sequence for motif search
def get_seq(fasta, contig_id, start, end):
    contig_seq = fasta[contig_id].seq.upper()
    motif_range = str(contig_seq[start - 1:end - 1])  # python slicing is 0-based

    return motif_range


# this function adds a report as a comment to the beginning of the new gff3 files
def prepend_to_file(file_path, text):
    with open(file_path, 'r') as f:
        old_content = f.read()

    with open(file_path, 'w') as f:
        f.write(text + '\n')  # Write the new content first
        f.write(old_content)  # Write the original content afterwards


if __name__ == "__main__":
    main()
