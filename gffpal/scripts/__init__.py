import sys
import argparse


from gffpal.exceptions import GPException, ECode
from gffpal.scripts.hints import hints, cli_hints
from gffpal.scripts.expandcds import expandcds, cli_expandcds
from gffpal.scripts.rnammer2gff import rnammer2gff, cli_rnammer2gff
from gffpal.scripts.trnascan2gff import trnascan2gff, cli_trnascan2gff


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
        else:
            raise ValueError("I shouldn't reach this point ever")
    except GPException as e:
        print(f"Error: {e.msg}")
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
