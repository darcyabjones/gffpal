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
    gene_parsed = parsed["gene"][0]
    similarity_parsed = parsed["similarity"][0]

    custom: Dict[str, str] = dict()
    if similarity_parsed.attributes is not None:
        custom["query"] = similarity_parsed.attributes.custom["Query"]

    if gene_parsed.attributes is not None:
        custom["identity"] = gene_parsed.attributes.custom["identity"]
        custom["similarity"] = gene_parsed.attributes.custom["similarity"]

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
            custom=custom,
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
                custom=(e.attributes.custom
                        if e.attributes is not None
                        else None)
            )
        )
        for e
        in parsed["exon"]
    ]

    for c in cdss:
        if gene.attributes is not None:
            # This is safe because we added attributes.
            assert c.attributes is not None
            c.attributes.custom["query"] = gene.attributes.custom["query"]

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

    if gene.attributes is not None:
        # This is safe because infer_from_children adds an ID to attributes.
        assert mrna.attributes is not None
        if gene.attributes.id is not None:
            mrna.attributes.parent = [gene.attributes.id]
        mrna.attributes.custom["query"] = gene.attributes.custom["query"]

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
            assert block is not None
            dealt_with = deal_with_block(block, gene_num)
            print("\n".join(str(f) for f in dealt_with), file=args.outfile)
            print("###", file=args.outfile)

            block = None
            gene_num += 1

        elif not line.startswith("#"):
            assert block is not None, line
            block.append(line)

    return
