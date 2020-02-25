import sys
import argparse

from typing import List
from typing import Dict
from typing import Optional

from gffpal.gff import GFFRecord, GFF3Record
from gffpal.attributes import GFF3Attributes, GTFAttributes


def cli_exonerate2gff(parser):
    parser.add_argument(
        "infile",
        type=argparse.FileType('r'),
        help="The exonerate gff2 file to convert."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the tidied gff to."
    )

    return


def deal_with_block(block: List[str], gene_num: int) -> List[GFF3Record]:

    parsed: Dict[str, List[GFFRecord[GTFAttributes]]] = dict()
    for line in block:
        rec = GFFRecord.parse(line, attr=GTFAttributes)

        if rec.type in parsed:
            parsed[rec.type].append(rec)
        else:
            parsed[rec.type] = [rec]

    assert len(parsed["gene"]) == 1
    assert len(parsed["similarity"]) == 1

    gene = GFF3Record(
        parsed["gene"][0].seqid,
        "exonerate",
        type="gene",
        start=parsed["gene"][0].start,
        end=parsed["gene"][0].end,
        score=parsed["gene"][0].score,
        strand=parsed["gene"][0].strand,
        phase=parsed["gene"][0].phase,
        attributes=GFF3Attributes(
            id=f"gene{gene_num}",
            custom={
                "query": parsed["similarity"][0].attributes.custom["Query"],
                "identity": parsed["gene"][0].attributes.custom["identity"],
                "similarity": parsed["gene"][0].attributes.custom["similarity"]
            }
        )
    )

    cdss = [
        GFF3Record(
            e.seqid,
            "exonerate",
            "CDS",
            e.start,
            e.end,
            e.score,
            e.strand,
            e.phase,
            attributes=GFF3Attributes(
                id=f"CDS{gene_num}",
                parent=[f"mRNA{gene_num}"],
                custom=e.attributes.custom
            )
        )
        for e
        in parsed["exon"]
    ]

    mrna = GFF3Record.infer_from_children(
        cdss,
        id=f"mRNA{gene_num}",
        seqid=gene.seqid,
        source="exonerate",
        type="mRNA",
        strand=gene.strand,
        score=gene.score,
    )

    mrna.add_parent(gene)
    mrna.attributes.parent = [gene.attributes.id]
    out = [gene, mrna]
    out.extend(cdss)
    return out


def exonerate2gff(args: argparse.Namespace) -> None:

    gene_num = 1

    print("##gff-version 3", file=args.outfile)

    block: Optional[List[str]] = None
    for line in args.infile:
        line = line.strip()
        if (line.startswith("Command line: ") or
                line.startswith("Hostname: ") or
                line.startswith("--")):
            continue

        if line.startswith("# --- START OF GFF DUMP ---"):
            assert block is None
            block = []
        elif line.startswith("# --- END OF GFF DUMP ---"):
            dealt_with = deal_with_block(block, gene_num)
            print("\n".join(str(f) for f in dealt_with), file=args.outfile)
            print("###", file=args.outfile)

            block = None
            gene_num += 1

        elif not line.startswith("#"):
            assert block is not None, line
            block.append(line)

    return
