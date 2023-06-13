import argparse
import re
import gffutils
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Validates false predictions flagged by liftoff")

    parser.add_argument("-f", "--fasta_file", type=str, required=True, help="Input path for genome assembly in fasta")
    parser.add_argument("-p", "--prediction", type=str, required=True, help="Input path for predictions in gff3")

    args = parser.parse_args()

    # load fasta
    fasta = SeqIO.index(args.fasta_file, "fasta")
    # parse gff3 file and store in a sqlite3 database
    db = gffutils.create_db(args.prediction, dbfn=":memory:", id_spec=id_func)
    # counter for number of motifs found
    motif_counter = 0

    for gene in db.features_of_type("gene"):
        # validation
        if "valid_orf" and "missing_start_codon" not in gene.attributes:
            continue

        valid_orf = gene.attributes["valid_ORF"]  # for some reason this returns a list with a single item
        missing_start_codon = gene.attributes["missing_start_codon"]

        if valid_orf[0] == "False":

            if missing_start_codon[0] == "True":

                if gene.strand == "+":
                    search_start = int(gene.end)
                    search_end = search_start + 15
                    motif_range = get_seq(fasta, gene.seqid, search_start, search_end)

                    if re.search(r"T[ACTG]{2}[ACTG]{2}TGTTTGTT", motif_range):  # check if the motif exists
                        gene.attributes["motif"] = "True"
                        motif_counter += 1
                        print(gene)
                    else:
                        gene.attributes["motif"] = "False"

                elif gene.strand == "-":
                    search_end = int(gene.start)
                    search_start = search_end - 15
                    motif_range = get_seq(fasta, gene.seqid, search_start, search_end)

                    if re.search(r"AACAAACA[ACTG]{2}[ACTG]{2}A", motif_range):  # check if the motif exists
                        gene.attributes["motif"] = ["True"]
                        motif_counter += 1
                        print(gene)
                    else:
                        gene.attributes["motif"] = ["False"]

        # gene.attributes["ID"] = gene.attributes["ID"][0] + "-1"



    # write to gff3
    # with open("test.gff3", "w") as out:
    #     out.write(f"# There were {motif_counter} transferred predictions with the motif found\n")
    #     for gene in db.features_of_type("gene", order_by=("seqid", "start")):
    #         out.write("\n")
    #         out.write(f"{str(gene)}\n")
    #         for child in db.children(gene, order_by="start"):
    #             out.write(f"{str(child)}\n")

    # status
    print(f"\nThere were {motif_counter} transferred predictions with the motif found")


def id_func(x):
    new_id = None

    if x.featuretype == 'gene':
        new_id = x.attributes['ID'][0]
    elif x.featuretype == 'CDS':
        new_id = '-'.join([x.attributes['ID'][0], x.seqid, str(x.start), str(x.end)])

    return new_id


def get_seq(fasta, contig_id, start, end):
    contig_seq = fasta[contig_id].seq.upper()
    motif_range = str(contig_seq[start - 1:end - 1])  # python slicing is 0-based

    return motif_range


if __name__ == "__main__":
    main()
