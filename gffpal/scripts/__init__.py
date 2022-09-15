import sys
import argparse
import logging

from gffpal.exceptions import GPException, ECode
from gffpal.scripts.hints import hints, cli_hints
from gffpal.scripts.expandcds import expandcds, cli_expandcds
from gffpal.scripts.rnammer2gff import rnammer2gff, cli_rnammer2gff
from gffpal.scripts.trnascan2gff import trnascan2gff, cli_trnascan2gff
from gffpal.scripts.exonerate2gff import exonerate2gff, cli_exonerate2gff
from gffpal.scripts.add_antifam import add_antifam, cli_add_antifam
from gffpal.scripts.add_parents import add_parents, cli_add_parents
from gffpal.scripts.select import select, cli_select
from gffpal.scripts.coords_to_contigs import coords, cli_coords
from gffpal.scripts.antismash2gff import antismash2gff, cli_antismash2gff

logging.basicConfig(level=logging.ERROR)


def cli(prog, args):

    parser = argparse.ArgumentParser(
        prog=prog,
        description=""
    )

    subparsers = parser.add_subparsers(dest='subparser_name')

    hints_subparser = subparsers.add_parser(
        "hints",
        help="Construct augustus hints files from various gff formats."
    )

    cli_hints(hints_subparser)

    expandcds_subparser = subparsers.add_parser(
        "expandcds",
        help="Construct expand CDS features to include stop or start codons."
    )

    cli_expandcds(expandcds_subparser)

    rnammer2gff_subparser = subparsers.add_parser(
        "rnammer2gff",
        help="Convert RNAmmer output to a gff3 format."
    )

    cli_rnammer2gff(rnammer2gff_subparser)

    trnascan2gff_subparser = subparsers.add_parser(
        "trnascan2gff",
        help="Convert TRNAScan output to a gff3 format."
    )

    cli_trnascan2gff(trnascan2gff_subparser)

    exonerate2gff_subparser = subparsers.add_parser(
        "exonerate2gff",
        help="Convert exonerate output to a gff3 format."
    )

    cli_exonerate2gff(exonerate2gff_subparser)

    add_antifam_subparser = subparsers.add_parser(
        "add_antifam",
        help="Add antifam matches to gff3 attributes."
    )

    cli_add_antifam(add_antifam_subparser)

    add_parents_subparser = subparsers.add_parser(
        "add_parents",
        help="Add parents that should be there for genes to gff3."
    )

    cli_add_parents(add_parents_subparser)

    select_subparser = subparsers.add_parser(
        "select",
        help="Select a subset based on their ids."
    )

    cli_select(select_subparser)

    coord_subparser = subparsers.add_parser(
        "coord2contig",
        help="Get some contigs given an alignment"
    )

    cli_coords(coord_subparser)

    antismash2gff_subparser = subparsers.add_parser(
        "antismash2gff",
        help="Convert the AntiSMASH genbank output to gff3"
    )

    cli_antismash2gff(antismash2gff_subparser)

    parsed = parser.parse_args(args)

    if parsed.subparser_name is None:
        parser.print_help()
        sys.exit(0)

    return parsed


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])
    try:
        if args.subparser_name == "hints":
            hints(args)
        elif args.subparser_name == "expandcds":
            expandcds(args)
        elif args.subparser_name == "rnammer2gff":
            rnammer2gff(args)
        elif args.subparser_name == "trnascan2gff":
            trnascan2gff(args)
        elif args.subparser_name == "exonerate2gff":
            exonerate2gff(args)
        elif args.subparser_name == "add_antifam":
            add_antifam(args)
        elif args.subparser_name == "add_parents":
            add_parents(args)
        elif args.subparser_name == "select":
            select(args)
        elif args.subparser_name == "coord2contig":
            coords(args)
        elif args.subparser_name == "antismash2gff":
            antismash2gff(args)
        else:
            raise ValueError("I shouldn't reach this point ever")
    except GPException as e:
        print(f"Error: {str(e)}")
        sys.exit(e.ecode)
    except BrokenPipeError:
        # Pipes get closed and that's normal
        sys.exit(0)
    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(ECode.SIGINT)
    except EnvironmentError as e:
        print((
            "Encountered a system error.\n"
            "We can't control these, and they're usually related to your OS.\n"
            "Try running again."
        ), file=sys.stderr)
        raise e
    except Exception as e:
        print((
            "I'm so sorry, but we've encountered an unexpected error.\n"
            "This shouldn't happen, so please file a bug report with the "
            "authors.\nWe will be extremely grateful!\n\n"
        ), file=sys.stderr)
        raise e
    return
