import sys
import argparse

"""
def cli(prog, args):

    parser = argparse.ArgumentParser(
        prog=prog,
        description=""
    )

    subparsers = parser.add_subparsers(dest='subparser_name')

    encode_subparser = subparsers.add_parser(
        "encode",
        help="Encode a files ids using some strategy."
    )

    cli_encode(encode_subparser)

    decode_subparser = subparsers.add_parser(
        "decode",
        help=("Decodes a file given a mapping file.")
    )

    cli_decode(decode_subparser)

    parsed = parser.parse_args(args)

    if parsed.subparser_name is None:
        parser.print_help()
        sys.exit(0)

    return parsed


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])
    try:
        if args.subparser_name == "encode":
            encode(args)
        elif args.subparser_name == "decode":
            decode(args)
        else:
            raise ValueError("I shouldn't reach this point ever")
    except SRException as e:
        print(f"Error: {e.msg}")
        sys.exit(e.ecode)
    except BrokenPipeError:
        # Pipes get closed and that's normal
        sys.exit(0)
    except KeyboardInterrupt:
        print("Received keyboard interrupt. Exiting.", file=sys.stderr)
        sys.exit(1)
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
"""
