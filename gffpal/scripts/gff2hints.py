import sys
import argparse

from gffpal.gff import GFF
from gffpal.attributes import GFFAttributes

SOURCES = ["M", "E", "P", "RM", "W", "XNT", "C", "D", "T", "R", "PB"]

GFF_TYPE_MAP = {
    "gene": "genicpart",
    "mRNA": "genicpart",
    "transcription_start_site": "transcription_start_site",
    "TSS": "transcription_start_site",
    "SO:0000315": "transcription_start_site",
    "transcription_end_site": "transcription_end_site",
    "SO:0000616": "transcription_end_site",
    "exon": "exon",
    "SO:0000147": "exon",
    "interior_coding_exon": "exon",
    "SO:0000004": "exon",
    "coding_exon": "exon",
    "SO:0000195": "exon",
    "five_prime_coding_exon": "exon",
    "SO:0000200": "exon",
    "three_prime_coding_exon": "exon",
    "SO:0000202": "exon",
    "interior_exon": "exon",
    "SO:0000201": "exon",
    "UTR": "UTR",
    "SO:0000203": "UTR",
    "noncoding_exon": "UTR",
    "SO:0000198": "UTR",
    "three_prime_UTR": "three_prime_UTR",
    "SO:0000205": "three_prime_UTR",
    "five_prime_UTR": "five_prime_UTR",
    "SO:0000204": "five_prime_UTR",
    "three_prime_noncoding_exon": "three_prime_UTR",
    "SO:0000444": "three_prime_UTR",
    "five_prime_noncoding_exon": "five_prime_UTR",
    "SO:0000445": "five_prime_UTR",
    "intron": "intron",
    "SO:0000188": "intron",
    "five_prime_intron": "intron",
    "SO:0000190": "intron",
    "interior_intron": "intron",
    "SO:0000191": "intron",
    "three_prime_intron": "intron",
    "SO:0000192": "intron",
    "UTR_intron": "intron",
    "SO:0000446": "intron",
    "SO:0000447": "intron",
    "five_prime_UTR_intron": "intron",
    "SO:0000448": "intron",
    "three_prime_UTR_intron": "intron",
    "CDS": "CDS",
    "coding_sequence": "CDS",
    "SO:0000316": "CDS",
    "CDS_fragment": "CDSpart",
    "SO:0001384": "CDSpart",
    "CDS_supported_by_peptide_spectrum_match": "CDSpart",
    "SO:0002071": "CDSpart",
    "start_codon": "start_codon",
    "SO:0000318": "start_codon",
    "non_canonical_start_codon": "start_codon",
    "SO:0000680": "start_codon",
    "stop_codon": "stop_codon",
    "SO:0000319": "stop_codon",
    "five_prime_cis_splice_site": "five_prime_cis_splice_site",
    "SO:0000163": "five_prime_cis_splice_site",
    "donor splice site": "five_prime_cis_splice_site",
    "five prime splice site": "five_prime_cis_splice_site",
    "splice donor site": "five_prime_cis_splice_site",
    "canonical_five_prime_splice_site": "five_prime_cis_splice_site",
    "SO:0000677": "five_prime_cis_splice_site",
    "non_canonical_five_prime_splice_site": "five_prime_cis_splice_site",
    "SO:0000679": "five_prime_cis_splice_site",
    "three_prime_cis_splice_site": "three_prime_cis_splice_site",
    "SO:0000164": "three_prime_cis_splice_site",
    "acceptor splice site": "three_prime_cis_splice_site",
    "splice acceptor site": "three_prime_cis_splice_site",
    "three prime splice site": "three_prime_cis_splice_site",
    "canonical_three_prime_splice_site": "three_prime_cis_splice_site",
    "SO:0000676": "three_prime_cis_splice_site",
    "non_canonical_three_prime_splice_site": "three_prime_cis_splice_site",
    "SO:0000678": "three_prime_cis_splice_site",
    "match": "CDS",
    "SO:0000343": "CDS",
    "protein_match": "CDS",
    "SO:0000349": "CDS",
    "nucleotide_to_protein_match": "CDS",
    "translated_nucleotide_match": "CDS",
    "SO:0000181": "CDS",
    "nucleotide_match": "exon",
    "SO:0000347": "exon",
    "cDNA_match": "exon",
    "SO:0000689": "exon",
}

HINT_TYPE = [
    "start",
    "stop",
    "tss",
    "tts",
    "ass",
    "dss",
    "exonpart",
    "exon",
    "intronpart",
    "intron",
    "CDSpart",
    "CDS",
    "UTRpart",
    "UTR",
    "irpart",
    "nonexonpart",
    "genicpart",
]


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        """
    )

    parser.add_argument(
        "infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help="Input fasta files. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output fasta file path. Default stdout.",
    )

    parser.add_argument(
        "-s", "--source",
        default="M",
        help=f"The type of hint to create. Usually one of {SOURCES}.",
    )

    parser.add_argument(
        "-p", "--priority",
        default=1,
        type=int,
        help="The priority to give all hints.",
    )

    parser.add_argument(
        "-g", "--group-level",
        default="mRNA",
        help="The level to group features at.",
    )

    parser.add_argument(
        "-c", "--cds",
        default="CDSpart",
        choices=HINT_TYPE,
        help="The type to map CDS features to."
    )

    parser.add_argument(
        "-i", "--intron",
        default="intron",
        choices=HINT_TYPE,
        help="The type to map intron features to."
    )

    parser.add_argument(
        "-e", "--exon",
        default="exonpart",
        choices=HINT_TYPE,
        help="The type to map exon features to."
    )

    parser.add_argument(
        "-5", "--utr5",
        default="UTRpart",
        choices=HINT_TYPE,
        help="The type to map five_prime_UTR features to."
    )

    parser.add_argument(
        "-3", "--utr3",
        default="UTRpart",
        choices=HINT_TYPE,
        help="The type to map three_prime_UTR features to."
    )

    parser.add_argument(
        "-u", "--utr",
        default="UTRpart",
        choices=HINT_TYPE,
        help="The type to map UTR features to."
    )

    parser.add_argument(
        "-f", "--feature",
        nargs="+",
        help="Pairs to map between.",
    )

    parser.add_argument(
        "--cds-trim",
        default=6,
        type=int,
        help="Trim cds hints by this many basepairs.",
    )

    parser.add_argument(
        "--intron-trim",
        default=0,
        type=int,
        help="Trim intronpart hints by this many basepairs.",
    )

    parser.add_argument(
        "--exon-trim",
        default=0,
        type=int,
        help="Trim exon hints by this many basepairs.",
    )

    parser.add_argument(
        "--utr-trim",
        default=0,
        type=int,
        help="Trim utr hints by this many basepairs.",
    )

    parser.add_argument(
        "--ir-trim",
        default=50,
        type=int,
        help="Trim genepart hints by this many basepairs.",
    )

    parser.add_argument(
        "--nonexon-trim",
        default=0,
        type=int,
        help="Trim nonexonpart hints by this many basepairs.",
    )

    parser.add_argument(
        "--gene-trim",
        default=9,
        type=int,
        help="Trim genepart hints by this many basepairs.",
    )

    parser.add_argument(
        "--cds-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    parser.add_argument(
        "--intron-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    parser.add_argument(
        "--exon-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    parser.add_argument(
        "--utr-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    parser.add_argument(
        "--ir-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    parser.add_argument(
        "--nonexon-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    parser.add_argument(
        "--gene-priority",
        default=0,
        type=int,
        help="Give this hint a priority.",
    )

    return parser.parse_args(args)


def parse_custom_features(features):
    if features is None or len(features) == 0:
        return {}

    assert len(features) % 2 == 0
    features = dict(f.split("=", maxsplit=1) for f in features)

    for feature in features.values():
        if feature not in HINT_TYPE:
            raise ValueError("custom features must map to valid hint type.")
    return features


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    gff_to_hints = {
        "CDS": args.cds,
        "CDSpart": "CDSpart",
        "intron": args.intron,
        "intronpart": "intronpart",
        "exon": args.exon,
        "exonpart": "exonpart",
        "five_prime_UTR": args.utr5,
        "three_prime_UTR": args.utr3,
        "UTR": args.utr,
        "UTRpart": "UTRpart",
        "transcription_start_site": "tss",
        "transcription_end_site": "tts",
        "start_codon": "start",
        "stop_codon": "stop",
        "five_prime_cis_splice_site": "dss",
        "three_prime_cis_splice_site": "ass",
        "irpart": "irpart",
        "nonexonpart": "nonexonpart",
        "genicpart": "genicpart",
    }

    gff_to_hints.update(parse_custom_features(args.feature))

    gff = GFF.parse(args.infile)
    for parent in gff.select_type(args.group_level):
        group_name = parent.attributes.id
        if group_name is None:
            raise ValueError("Parent doesn't have an id.")

        for feature in gff.traverse_children([parent]):
            if feature.type in gff_to_hints:
                type_ = gff_to_hints[feature.type]
            else:
                type_ = GFF_TYPE_MAP.get(feature.type, None)
                type_ = gff_to_hints.get(type_, None)

            if type_ is None:
                continue

            feature.type = type_

            if feature.type == "exonpart":
                feature.trim_ends(args.exon_trim)
                priority = args.exon_priority
            elif feature.type == "CDSpart":
                feature.trim_ends(args.cds_trim)
                priority = args.cds_priority
            elif feature.type == "UTRpart":
                feature.trim_ends(args.utr_trim)
                priority = args.utr_priority
            elif feature.type == "intronpart":
                feature.trim_ends(args.intron_trim)
                priority = args.intron_priority
            elif feature.type == "genicpart":
                feature.trim_ends(args.gene_trim)
                priority = args.gene_priority
            elif feature.type == "irpart":
                feature.trim_ends(args.ir_trim)
                priority = args.ir_priority
            elif feature.type == "nonexonpart":
                feature.trim_ends(args.nonexon_trim)
                priority = args.nonexon_priority

            attr = GFFAttributes(
                custom=dict(
                    source=args.source,
                    group=group_name,
                    priority=args.priority + priority
                )
            )
            feature.attributes = attr
            print(feature, file=args.outfile)
    return


if __name__ == "__main__":
    main()
