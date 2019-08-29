import sys
import argparse

from typing import List, Tuple, Dict
from typing import Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation
from Bio.Data import CodonTable

from gffpal.gff import GFF, GFFRecord, Strand, Phase
from gffpal.attributes import GFFAttributes


def cli_expandcds(parser):
    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The gff file fix the CDSs of.",
    )

    parser.add_argument(
        "-i", "--infasta",
        type=argparse.FileType('r'),
        default=None,
        help=(
            "The genome that the gff refers to. "
            "If this is provided, start and stops will be checked to "
            "see that they match the translation table, and will raise "
            "warnings if they don't."
        ),
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the tidied gff to."
    )

    parser.add_argument(
        "-g", "--gencode",
        type=int,
        default=1,
        help="The ncbi translation table number to use to determine start"
    )

    parser.add_argument(
        "-s", "--start",
        action="store_true",
        default=False,
        help="Expand the CDS to the start too. Default only shifts the stop."
    )

    parser.add_argument(
        "-n", "--no-stop",
        dest="stop",
        action="store_false",
        default=True,
        help="Don't expand the stop codon. Default only shifts the stop."
    )

    parser.add_argument(
        "-w", "--warnings",
        type=argparse.FileType("w"),
        default=sys.stderr,
        help="Where to write warnings to."
    )

    parser.add_argument(
        "--cds-type",
        type=str,
        default="CDS",
        help="The type in the GFF to treat as CDS (Case sensitive)."
    )

    return


def bump_start(cdss: List[GFFRecord], strand: Strand) -> None:
    assert len(cdss) > 0

    if strand == Strand.MINUS:
        cds = cdss[-1]
        cds.start += 3
    else:
        cds = cdss[0]
        cds.start -= 3

    return


def bump_end(cdss: List[GFFRecord], strand: Strand) -> None:
    assert len(cdss) > 0

    if strand == Strand.MINUS:
        cds = cdss[0]
        cds.start -= 3
    else:
        cds = cdss[-1]
        cds.end += 3

    return


def get_start(
    cdss: List[GFFRecord],
    seq: SeqRecord,
    strand: Strand,
) -> Tuple[str, int, int]:
    assert len(cdss) > 0

    if strand == Strand.MINUS:
        cds = cdss[-1]
        end = cds.end
        start = end - 3
        feat = FeatureLocation(start, end, -1)
    else:
        cds = cdss[0]
        start = cds.start
        end = start + 3
        feat = FeatureLocation(start, end, +1)

    return str(feat.extract(seq).seq), start, end


def get_stop(
    cdss: List[GFFRecord],
    seq: SeqRecord,
    strand: Strand,
) -> Tuple[str, int, int]:
    assert len(cdss) > 0

    if strand == Strand.MINUS:
        cds = cdss[0]
        start = cds.start
        end = start + 3
        feat = FeatureLocation(start, end, -1)
    else:
        cds = cdss[-1]
        end = cds.end
        start = end - 3
        feat = FeatureLocation(start, end, +1)

    return str(feat.extract(seq).seq), start, end


def check_start(
    cdss: List[GFFRecord],
    parent: GFFRecord,
    strand: Strand,
    seqs: Dict[str, SeqRecord],
    codon_table: CodonTable,
):
    codon, start, end = get_start(cdss, seqs[parent.seqid], strand)

    if codon in codon_table.start_codons:
        return None
    else:
        return get_non_canon_start_codon(
            parent.seqid,
            start,
            end,
            strand,
            codon,
            parent.attributes.id
        )


def check_stop(
    cdss: List[GFFRecord],
    parent: GFFRecord,
    strand: Strand,
    seqs: Dict[str, SeqRecord],
    codon_table: CodonTable,
) -> Optional[GFFRecord]:
    codon, start, end = get_stop(cdss, seqs[parent.seqid], strand)

    if codon in codon_table.stop_codons:
        return None
    else:
        return get_non_canon_stop_codon(
            parent.seqid,
            start,
            end,
            strand,
            codon,
            parent.attributes.id
        )


def find_strand(cdss: List[GFFRecord], parent: GFFRecord) -> Strand:
    strands = list({f.strand for f in cdss})
    if len(strands) > 1:
        return parent.strand
    else:
        return strands[0]


def get_non_canon_start_codon(
    seqid: str,
    start: int,
    end: int,
    strand: Strand,
    codon: str,
    parent_id: str,
) -> GFFRecord:

    return GFFRecord(
        seqid,
        "gffpal",
        "non_canonical_start_codon",
        start,
        end,
        None,
        strand,
        Phase.NOT_CDS,
        GFFAttributes(
            ontology_term=["SO:0000680"],
            custom={
                "codon": codon,
                "cds_parent": parent_id,
            }
        )
    )


def get_non_canon_stop_codon(
    seqid: str,
    start: int,
    end: int,
    strand: Strand,
    codon: str,
    parent_id: str,
) -> GFFRecord:

    return GFFRecord(
        seqid,
        "gffpal",
        "stop_codon",
        start,
        end,
        None,
        strand,
        Phase.NOT_CDS,
        GFFAttributes(
            ontology_term=["SO:0000319"],
            note=["Non-canonical stop codon"],
            custom={
                "codon": codon,
                "cds_parent": parent_id,
            }
        )
    )


def expandcds(args: argparse.Namespace) -> None:

    gff = GFF.parse(args.infile)

    if args.infasta is None:
        seqs = None
    else:
        seqs = SeqIO.to_dict(SeqIO.parse(args.infasta, format="fasta"))

    codon_table = CodonTable.unambiguous_dna_by_id[args.gencode]

    cds_parents = set()
    for record in gff.select_type(args.cds_type):
        cds_parents.update(record.parents)

    for parent in cds_parents:
        cdss = sorted(
            [f for f in parent.children if f.type == args.cds_type],
            key=lambda f: (f.start, f.end)
        )

        strand = find_strand(cdss, parent)

        if args.start:
            bump_start(cdss, strand)

        if args.stop:
            bump_end(cdss, strand)

        if seqs is not None and parent.seqid in seqs:

            start_codon = check_start(
                cdss,
                parent,
                strand,
                seqs,
                codon_table
            )

            if start_codon is not None:
                print(start_codon, file=args.warnings)

            stop_codon = check_stop(
                cdss,
                parent,
                strand,
                seqs,
                codon_table
            )

            if stop_codon is not None:
                print(stop_codon, file=args.warnings)

    for parent in gff.traverse_parents(gff.select_type(args.cds_type)):
        parent.expand_to_children()

    print("##gff-version 3", file=args.outfile)
    for feature in gff.traverse_children(sort=True):
        print(feature, file=args.outfile)
    return
