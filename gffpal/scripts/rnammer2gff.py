import sys
import argparse
from copy import deepcopy

from gffpal.gff import GFF
from gffpal.gff import GFFRecord

import logging
logger = logging.getLogger(__name__)

TYPE_MAP = {
    "euk": {
        "5s_rrna": "rRNA_5S",
        "8s_rrna": "rRNA_5S",
        "18s_rrna": "rRNA_18S",
        "28s_rrna": "rRNA_28S",
    },
    "bac": {
        "5s_rrna": "rRNA_5S",
        "16s_rrna": "rRNA_16S",
        "23s_rrna": "rRNA_23S",
    },
    "arc": {
        "5s_rrna": "rRNA_5S",
        "16s_rrna": "rRNA_16S",
        "23s_rrna": "rRNA_23S",
    },
}


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        """
    )

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help=(
            "Input gff2 result file. "
            "Use '-' for stdin."
        ),
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output gff file path. Default stdout.",
    )

    parser.add_argument(
        "-s", "--source",
        default="RNAmmer",
        help=f"What to put in the source gff column.",
    )

    parser.add_argument(
        "-k", "--kingdom",
        default="euk",
        choices=["arc", "bac", "euk"],
        help="What kingdom was used to run rnammer?",
    )
    return parser.parse_args(args)


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])
    records = []
    for line in args.infile:
        if line.startswith("#"):
            continue

        sline = line.strip().split("\t")
        rrna_type = sline[8]
        new_type = TYPE_MAP[args.kingdom][rrna_type.lower()]
        sline[1] = args.source
        sline[2] = new_type
        sline[8] = "."

        rna_record = GFFRecord.parse("\t".join(sline))
        gene_record = deepcopy(rna_record)
        gene_record.type = "rRNA_gene"
        gene_record.add_child(rna_record)

        records.append(gene_record)
        records.append(rna_record)

    num = 0
    for record in GFF(records).traverse_children(sort=True):
        if record.type == "rRNA_gene":
            num += 1
            record.attributes.id = f"rRNA_gene{num}"
        else:
            record.attributes.id = f"rRNA{num}"
            record.attributes.parent = [
                p.attributes.id
                for p
                in record.parents
            ]
        print(record, file=args.outfile)

    return
