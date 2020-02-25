import sys
import argparse

from typing import Tuple, List, Dict
from typing import Sequence, Mapping

from gffpal.gff import GFF
from gffpal.gff import GFF3Record
from gffpal.gff import Strand, Phase
from gffpal.attributes import GFF3Attributes
from gffpal.parsers.trnascan import TRNAScanRecord, TRNAScanSS

import logging
logger = logging.getLogger(__name__)

TYPE_MAP: Dict[str, str] = {
    "ala": "alanyl_tRNA",
    "gln": "glutaminyl_tRNA",
    "pro": "prolyl_tRNA",
    "glu": "glutamyl_tRNA",
    "met": "methionyl_tRNA",
    "asn": "asparaginyl_tRNA",
    "thr": "threonyl_tRNA",
    "gly": "glycyl_tRNA",
    "val": "valyl_tRNA",
    "tyr": "tyrosyl_tRNA",
    "cys": "cysteinyl_tRNA",
    "iso": "isoleucyl_tRNA",
    "ser": "seryl_tRNA",
    "leu": "leucyl_tRNA",
    "trp": "tryptophanyl_tRNA",
    "sec": "selenocysteinyl_tRNA",
    "pyl": "pyrrolysyl_tRNA",
    "lys": "lysyl_tRNA",
    "asp": "aspartyl_tRNA",
    "arg": "arginyl_tRNA",
    "his": "histidyl_tRNA",
    "phe": "phenylalanyl_tRNA",
}


def cli_trnascan2gff(parser):
    parser.add_argument(
        "txt",
        type=argparse.FileType('r'),
        help="Input trnascan result file. Use '-' for stdin.",
    )

    parser.add_argument(
        "ss",
        type=argparse.FileType('r'),
        help=(
            "Input trnascan secondary structure result file."
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
        default="tRNAScan-SE",
        help=f"What to put in the source gff column.",
    )

    return parser


def fix_strand(start: int, end: int) -> Tuple[int, int, Strand]:
    if start > end:
        strand = Strand.MINUS
        tmp = end
        end = start
        start = tmp

    elif start < end:
        strand = Strand.PLUS

    else:
        strand = Strand.UNSTRANDED

    start -= 1
    return start, end, strand


def match_to_gene(
    match: TRNAScanRecord,
    source: str,
    type: str
) -> GFF3Record:
    start, end, strand = fix_strand(match.start, match.end)

    gene = GFF3Record(
        seqid=match.seqid,
        source=source,
        type=type,
        start=start,
        end=end,
        score=match.infernal_score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=GFF3Attributes(
            id=f"{match.seqid}.{type}{match.num}",
        )
    )

    return gene


def match_to_trna(
    match: TRNAScanRecord,
    ss: TRNAScanSS,
    source: str,
    type_map: Mapping[str, str] = TYPE_MAP,
    parents: Sequence[GFF3Record] = []
) -> GFF3Record:
    start, end, strand = fix_strand(match.start, match.end)

    parent_ids = [
        p.attributes.id
        for p
        in parents
        if (p.attributes is not None
            and p.attributes.id is not None)
    ]

    if match.note is None or match.note == "":
        notes: List[str] = []
    else:
        notes = [match.note]

    trna = GFF3Record(
        seqid=match.seqid,
        source=source,
        type=type_map.get(match.trna_type.lower(), "tRNA"),
        start=start,
        end=end,
        score=match.infernal_score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=GFF3Attributes(
            id=f"{match.seqid}.tRNA{match.num}",
            parent=parent_ids,
            note=notes,
            custom={
                "secondary_structure": ss.ss,
                "anticodon": match.anticodon,
                "amino_acid": match.trna_type,
            }
        ),
        parents=parents
    )
    return trna


def match_to_introns(
    match: TRNAScanRecord,
    source: str,
    type: str = "tRNA_intron",
    parents: Sequence[GFF3Record] = [],
) -> List[GFF3Record]:
    introns = []

    parent_ids = [
        p.attributes.id
        for p
        in parents
        if (p.attributes is not None
            and p.attributes.id is not None)
    ]

    for istart, iend in zip(match.intron_starts, match.intron_ends):
        start, end, strand = fix_strand(istart, iend)
        intron = GFF3Record(
            seqid=match.seqid,
            source=source,
            type=type,
            start=start,
            end=end,
            score=match.infernal_score,
            strand=strand,
            phase=Phase.NOT_CDS,
            attributes=GFF3Attributes(
                id=f"{match.seqid}.{type}{match.num}",
                parent=parent_ids,
            ),
            parents=parents
        )
        introns.append(intron)
    return introns


def match_to_anticodon(
    match: TRNAScanRecord,
    ss: TRNAScanSS,
    source: str,
    type: str = "anticodon",
    parents: Sequence[GFF3Record] = []
) -> GFF3Record:
    start, end, strand = fix_strand(ss.anticodon_start, ss.anticodon_end)

    parent_ids = [
        p.attributes.id
        for p
        in parents
        if (p.attributes is not None
            and p.attributes.id is not None)
    ]

    anticodon = GFF3Record(
        seqid=match.seqid,
        source=source,
        type=type,
        start=start,
        end=end,
        score=match.infernal_score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=GFF3Attributes(
            id=f"{match.seqid}.{type}{match.num}",
            parent=parent_ids,
        ),
        parents=parents
    )
    return anticodon


def trnascan2gff(args: argparse.Namespace) -> None:
    genes: List[GFF3Record] = []

    matches = TRNAScanRecord.parse(args.txt)
    sses = TRNAScanSS.parse(args.ss)
    num_to_ss = {f"{r.seqid}.{r.num}": r for r in sses}

    for match in matches:
        ss = num_to_ss[f"{match.seqid}.{match.num}"]

        if match.note is not None and "pseudo" in match.note:
            type_ = "pseudogene"
        else:
            type_ = "tRNA_gene"

        gene = match_to_gene(match, args.source, type=type_)
        genes.append(gene)

        trna = match_to_trna(
            match,
            ss,
            args.source,
            type_map=TYPE_MAP,
            parents=[gene]
        )
        genes.append(trna)

        introns = match_to_introns(
            match,
            args.source,
            type="tRNA_intron",
            parents=[trna]
        )
        genes.extend(introns)

        anticodon = match_to_anticodon(
            match,
            ss,
            args.source,
            type="anticodon",
            parents=[trna]
        )
        genes.append(anticodon)

    for record in GFF(genes).traverse_children(sort=True):
        print(record, file=args.outfile)

    return
