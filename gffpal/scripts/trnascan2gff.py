import sys
import argparse

from gffpal.gff import GFF
from gffpal.gff import GFFRecord
from gffpal.gff import Strand, Phase
from gffpal.attributes import GFFAttributes
from gffpal.parsers.trnascan import TRNAScanRecord, TRNAScanSS

import logging
logger = logging.getLogger(__name__)

TYPE_MAP = {
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


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        """
    )

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
    return parser.parse_args(args)


def fix_strand(start, end):
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


def match_to_gene(match, source, type):
    start, end, strand = fix_strand(match.start, match.end)

    gene = GFFRecord(
        seqid=match.seqid,
        source=source,
        type=type,
        start=start,
        end=end,
        score=match.infernal_score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=GFFAttributes(
            id=f"{type}{match.num}",
        )
    )

    return gene


def match_to_trna(match, ss, source, type_map=TYPE_MAP, parents=[]):
    start, end, strand = fix_strand(match.start, match.end)
    trna = GFFRecord(
        seqid=match.seqid,
        source=source,
        type=type_map.get(match.trna_type.lower(), "tRNA"),
        start=start,
        end=end,
        score=match.infernal_score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=GFFAttributes(
            id=f"tRNA{match.num}",
            parent=[p.attributes.id for p in parents],
            custom={
                "secondary_structure": ss.ss,
                "anticodon": match.anticodon,
                "amino_acid": match.trna_type,
            }
        ),
        parents=parents
    )
    return trna


def match_to_introns(match, source, type="tRNA_intron", parents=[]):
    introns = []
    for istart, iend in zip(match.intron_starts, match.intron_ends):
        start, end, strand = fix_strand(istart, iend)
        intron = GFFRecord(
            seqid=match.seqid,
            source=source,
            type=type,
            start=start,
            end=end,
            score=match.infernal_score,
            strand=strand,
            phase=Phase.NOT_CDS,
            attributes=GFFAttributes(
                id=f"{type}{match.num}",
                parent=[p.attributes.id for p in parents],
            ),
            parents=parents
        )
        introns.append(intron)
    return introns


def match_to_anticodon(match, ss, source, type="anticodon", parents=[]):
    start, end, strand = fix_strand(ss.anticodon_start, ss.anticodon_end)
    anticodon = GFFRecord(
        seqid=match.seqid,
        source=source,
        type=type,
        start=start,
        end=end,
        score=match.infernal_score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=GFFAttributes(
            id=f"{type}{match.num}",
            parent=[p.attributes.id for p in parents],
        ),
        parents=parents
    )
    return anticodon


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])
    genes = []

    matches = TRNAScanRecord.parse(args.txt)
    sses = TRNAScanSS.parse(args.ss)
    num_to_ss = {r.num: r for r in sses}

    for match in matches:
        ss = num_to_ss[match.num]

        gene = match_to_gene(match, args.source, type="tRNA_gene")
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
