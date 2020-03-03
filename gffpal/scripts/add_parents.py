import sys
import argparse

from gffpal.gff import GFF, GFF3Record


def cli_add_parents(parser):
    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The gff file to add parents to.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the tidied gff to."
    )

    return


def add_parents(args: argparse.Namespace) -> None:

    gff = GFF.parse(args.infile)
    gff.infer_missing_parents()

    for f in gff.select_type("mRNA"):
        if len(f.parents) > 0:
            continue

        if f.attributes is None:
            continue

        if f.attributes.id is None:
            continue

        id_ = f.attributes.id
        gene_id = f"gene.{id_}"
        gene = GFF3Record.infer_from_children([f], id=gene_id, type="gene")
        f.add_parent(gene)
        gff.add_record(gene)

    print("##gff-version 3", file=args.outfile)
    for feature in gff.traverse_children(sort=True):
        print(feature, file=args.outfile)
    return
