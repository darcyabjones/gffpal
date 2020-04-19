import sys
import argparse
from copy import deepcopy

from typing import cast
from typing import Set, Dict
from typing import Iterator

from gffpal.gff import GFF, GFF3Record


def cli_select(parser):

    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The gff file fix the CDSs of.",
    )

    parser.add_argument(
        "--so",
        type=str,
        default="https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/master/Ontology_Files/so.obo",  # noqa
        help="The sequence ontology OBO file.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the tidied gff to."
    )

    return


def prune_gff(records: Set[GFF3Record]) -> GFF:

    new_records: Dict[GFF3Record, GFF3Record] = dict()

    # Create a mapping from old to new objects
    # to preserve hashing/lookup capability
    for record in records:
        new_record = deepcopy(record)
        new_record.children = []
        new_record.parents = []

        new_records[record] = new_record

    for record in records:
        new_record = new_records[record]

        for parent in record.parents:
            # Don't add parents that shouldn't be in this set
            if parent not in records:
                continue

            new_parent = new_records[cast(GFF3Record, parent)]
            new_record.add_parent(new_parent)

        for child in record.children:
            # Don't add children that shouldn't be in this set
            if child not in records:
                continue

            new_child = new_records[cast(GFF3Record, child)]
            new_record.add_child(new_child)

        # Update the record parent IDS to reflect the new split set.
        if new_record.attributes is not None:
            new_record.attributes.parent = []
            for parent in new_record.parents:
                # This should always be true, as the ID is necessary
                assert parent.attributes is not None
                assert parent.attributes.id is not None
                new_record.attributes.parent.append(parent.attributes.id)
        else:
            # This necessarily should be true, since attributes define parent
            # child relationships.
            assert len(new_record.children) == 0
            assert len(new_record.parents) == 0

    return GFF(list(new_records.values()))


def select(args: argparse.Namespace) -> None:

    gff = GFF.parse(args.infile)
    ids = {l.strip() for l in args.ids}

    to_keep: Set[GFF3Record] = set()

    for record in gff:
        if record.attributes is not None and record.attributes.id in ids:
            to_keep.update(cast(
                Iterator[GFF3Record],
                record.traverse_parents()
            ))

            to_keep.update(cast(
                Iterator[GFF3Record],
                record.traverse_children()
            ))

    pruned = prune_gff(to_keep)

    print("#gff-version 3", file=args.outfile)
    for feature in pruned.traverse_children(sort=True):
        print(feature, file=args.outfile)
