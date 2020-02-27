import sys
import argparse
from collections import defaultdict

from typing import Dict, List

from gffpal.gff import GFF3Record
from gffpal.parsers.domtbl import DomTbl

import logging
logger = logging.getLogger(__name__)


def cli_add_antifam(parser):
    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help=(
            "Input gff3 result file. "
            "Use '-' for stdin."
        ),
    )

    parser.add_argument(
        "antifam",
        type=argparse.FileType('r'),
        help=(
            "Input domtbl result file. "
            "Use '-' for stdin."
        ),
    )

    parser.add_argument(
        "-f", "--field",
        type=str,
        default="ID",
        help="Align matches using this field in the attributes."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output gff file path. Default stdout.",
    )
    return


def add_antifam(args: argparse.Namespace) -> None:

    antifam_records: Dict[str, List[DomTbl]] = defaultdict(list)
    for rec in DomTbl.from_file(args.antifam):
        antifam_records[rec.target_name].append(rec)

    for line in args.infile:
        sline = line.strip()
        if sline.startswith("#") or sline == "":
            args.outfile.write(line)
            continue

        record = GFF3Record.parse(sline)
        if record.attributes is None:
            print(record, file=args.outfile)
            continue

        field = record.attributes.get(args.field, None)
        if field is None:
            print(record, file=args.outfile)
            continue

        if str(field) not in antifam_records:
            print(record, file=args.outfile)
            continue

        dbxrefs = []
        matches = []
        for antifam_record in antifam_records[str(field)]:
            dbxrefs.append(f"AntiFam:{antifam_record.query_acc}")
            matches.append(
                f"{antifam_record.query_acc} {antifam_record.full_evalue} "
                f"{antifam_record.full_score} {antifam_record.domain_score} "
                f"{antifam_record.hmm_from} {antifam_record.hmm_to} "
                f"{antifam_record.ali_from} {antifam_record.ali_to}"
            )

        record.attributes.dbxref.extend(dbxrefs)
        record.attributes.custom["antifam_match"] = ",".join(matches)
        print(record, file=args.outfile)

    return
