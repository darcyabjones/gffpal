import sys
import argparse
from copy import copy

from typing import List
from typing import Dict
from typing import Optional

from gffpal.gff import GFF, GFFRecord
from gffpal.attributes import GFFAttributes
from gffpal.exceptions import GPCLIError, GPMissingID
from gffpal.higher import fmap, applicative

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


def cli_hints(parser):
    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="Input gff3 file. Use '-' for stdin.",
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output gff3 hints file path. Default stdout.",
    )

    parser.add_argument(
        "-s", "--source",
        type=str,
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
        type=str,
        help="The level to group features at.",
    )

    parser.add_argument(
        "-c", "--cds",
        default="CDSpart",
        choices=HINT_TYPE,
        type=str,
        help="The type to map CDS features to."
    )

    parser.add_argument(
        "-i", "--intron",
        default="intron",
        choices=HINT_TYPE,
        type=str,
        help="The type to map intron features to."
    )

    parser.add_argument(
        "-e", "--exon",
        default="exonpart",
        choices=HINT_TYPE,
        type=str,
        help="The type to map exon features to."
    )

    parser.add_argument(
        "-5", "--utr5",
        default="UTRpart",
        choices=HINT_TYPE,
        type=str,
        help="The type to map five_prime_UTR features to."
    )

    parser.add_argument(
        "-3", "--utr3",
        default="UTRpart",
        choices=HINT_TYPE,
        type=str,
        help="The type to map three_prime_UTR features to."
    )

    parser.add_argument(
        "-u", "--utr",
        default="UTRpart",
        choices=HINT_TYPE,
        type=str,
        help="The type to map UTR features to."
    )

    parser.add_argument(
        "-f", "--feature",
        nargs="+",
        type=str,
        default=None,
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
        "--start-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--stop-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--tss-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--tts-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--ass-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--dss-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--cds-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--intron-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--exon-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--utr-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--ir-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--nonexon-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    parser.add_argument(
        "--gene-priority",
        default=0,
        type=int,
        help="Give this hint a priority boost.",
    )

    return


def parse_custom_features(features: List[str]) -> Dict[str, str]:
    if features is None or len(features) == 0:
        return {}

    split_features: Dict[str, str] = dict()

    for feature in features:
        split_feature = feature.split("=", maxsplit=1)
        try:
            key = split_feature[0]
            val = split_feature[1]
        except IndexError:
            raise GPCLIError(
                "Custom features must have the format 'key=value'. "
                f"The offending feature was '{feature}'."
            )

        if val not in HINT_TYPE:
            raise GPCLIError(
                "Custom features must map to valid hint type. "
                f"The offending feature was '{feature[1]}'."
            )

        split_features[key] = val

    return split_features


def get_hints_map(args: argparse.Namespace) -> Dict[str, str]:
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
    return gff_to_hints


def get_trim_map(args: argparse.Namespace) -> Dict[str, int]:
    return {
        "exonpart": args.exon_trim,
        "CDSpart": args.cds_trim,
        "UTRpart": args.utr_trim,
        "intronpart": args.intron_trim,
        "genicpart": args.gene_trim,
        "irpart": args.ir_trim,
        "nonexonpart": args.nonexon_trim,
    }


def get_priority_map(args: argparse.Namespace) -> Dict[str, int]:
    return {
        "start": args.start_priority,
        "stop": args.stop_priority,
        "tss": args.tss_priority,
        "tts": args.tts_priority,
        "ass": args.ass_priority,
        "dss": args.dss_priority,
        "exonpart": args.exon_priority,
        "exon": args.exon_priority,
        "CDSpart": args.cds_priority,
        "CDS": args.cds_priority,
        "UTRpart": args.utr_priority,
        "UTR": args.utr_priority,
        "intronpart": args.intron_priority,
        "intron": args.intron_priority,
        "genicpart": args.gene_priority,
        "irpart": args.ir_priority,
        "nonexonpart": args.nonexon_priority,
    }


def transform_child(
    feature: GFFRecord,
    group_name: str,
    gff_to_hints: Dict[str, str],
    type_to_trim: Dict[str, int],
    type_to_priority: Dict[str, int],
    source: str,
    priority: int,
) -> Optional[GFFRecord]:
    """ Converts a regular feature to a hint record. """

    feature = copy(feature)
    if feature.type not in gff_to_hints:
        mapped_type = GFF_TYPE_MAP.get(feature.type, None)
    else:
        mapped_type = feature.type

    hint_type: Optional[str] = applicative(
        lambda t: gff_to_hints.get(t, None),
        mapped_type
    )

    if hint_type is None:
        return None

    feature.type = hint_type
    feature.trim_ends(type_to_trim[feature.type])
    priority_boost = type_to_priority[feature.type]

    attr = GFFAttributes(
        custom=dict(
            source=source,
            group=group_name,
            priority=str(priority + priority_boost)
        )
    )
    feature.attributes = attr
    return feature


def hints(args: argparse.Namespace) -> None:
    gff_to_hints = get_hints_map(args)

    type_to_trim = get_trim_map(args)
    type_to_priority = get_priority_map(args)

    gff = GFF.parse(args.infile)
    for parent in gff.select_type(args.group_level):
        group_name = fmap(lambda a: getattr(a, "id"), parent.attributes)

        if group_name is None:
            raise GPMissingID(
                "One of the selected records doesn't have an ID. "
                f"The offending line is {parent}."
            )

        for feature in gff.traverse_children([parent]):
            hint_feature = transform_child(
                feature,
                group_name,
                gff_to_hints,
                type_to_trim,
                type_to_priority,
                args.source,
                args.priority,
            )

            if hint_feature is not None:
                print(hint_feature, file=args.outfile)
    return
