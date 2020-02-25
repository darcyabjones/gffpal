import sys
import argparse
from copy import deepcopy

from typing import cast
from typing import List

from gffpal.gff import GFF
from gffpal.gff import GFF3Record
from gffpal.attributes import GFF3Attributes

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


def cli_rnammer2gff(parser):
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
    return


def rnammer2gff(args: argparse.Namespace) -> None:
    records: List[GFF3Record] = []

    for line in args.infile:
        if line.startswith("#"):
            continue

        sline = line.strip().split("\t")
        rrna_type = sline[8]
        new_type = TYPE_MAP[args.kingdom][rrna_type.lower()]
        sline[1] = args.source
        sline[2] = new_type
        sline[8] = "."

        rna_record = cast(GFF3Record, GFF3Record.parse("\t".join(sline)))
        gene_record = deepcopy(rna_record)
        gene_record.type = "rRNA_gene"
        gene_record.add_child(rna_record)

        records.append(gene_record)
        records.append(rna_record)

    num = 0
    for record in GFF(records).traverse_children(sort=True):
        if record.attributes is None:
            attr = GFF3Attributes()
            record.attributes = attr
        else:
            attr = record.attributes

        if record.type == "rRNA_gene":
            num += 1
            attr.id = f"rRNA_gene{num}"
        else:
            attr.id = f"rRNA{num}"
            attr.parent = [
                p.attributes.id
                for p
                in record.parents
                if (p.attributes is not None
                    and p.attributes.id is not None)
            ]

        print(record, file=args.outfile)

    return
