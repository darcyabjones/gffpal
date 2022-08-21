import sys
import argparse

from copy import deepcopy

from bio import SeqIO
from bio import SeqFeature

from typing import Iterable
from typing import List, Sequence, Tuple
from typing import Dict
from typing import Any, Optional
from typing import TextIO

from gffpal.gff import GFF3Record, Phase, Strand
from gffpal.attributes import GFF3Attributes, Target


def cli_antismash2gff(parser):
    parser.add_argument(
        "infiles",
        type=argparse.FileType('r'),
        nargs="+",
        help="The antismash GBK files to convert."
    )

    parser.add_argument(
        "-o", "--outfile",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Where to write the gff to."
    )

    return


def int_to_strand(i: int) -> Strand:
    if i == 1:
        return Strand.PLUS
    elif i == -1:
        return Strand.MINUS
    elif i == 0:
        return Strand.UNSTRANDED
    else:
        return Strand.UNKNOWN


def int_to_phase(i: int) -> Phase:
    return Phase.parse(str(i))


def find_source(cdss: Iterable[GFF3Record]) -> str:
    sources = {c.source for c in cdss}
    if len(sources) > 1:
        raise ValueError(f"Got multiple strands {sources}")

    return sources.pop()


def find_strand(cdss: Iterable[GFF3Record]) -> Strand:
    strands = {c.strand for c in cdss}
    if len(strands) > 1:
        raise ValueError(f"Got multiple strands {strands}")

    return strands.pop()


def find_seqid(cdss: Iterable[GFF3Record]) -> str:
    seqids = {c.seqid for c in cdss}
    if len(seqids) > 1:
        raise ValueError(f"Got multiple seqids {seqids}")

    return seqids.pop()


def calculate_phase(
    parts: List[SeqFeature],
    strand: Strand
) -> List[Tuple[SeqFeature, Phase]]:

    if strand == Strand.PLUS:
        parts = sorted(parts, key=lambda p: p.start.position)
    elif strand == Strand.MINUS:
        parts = sorted(parts, reverse=True, key=lambda p: p.start.position)
    else:
        raise ValueError("Must be plus or minus")

    phases = []
    lag = 0
    for part in parts:
        phases.append(int_to_phase(lag))
        length = (part.end.position - part.start.position - 1)
        length -= lag

        if strand == Strand.PLUS:
            lag = 2 - length % 3
        else:
            lag = 2 - length % 3

    return list(zip(parts, phases))


def pop_list(
    key: str,
    qualifiers: Dict[Any, List[Any]]
) -> Any:
    tmp = qualifiers.pop(key, [])
    if len(tmp) == 0:
        return None
    else:
        return tmp[0]


def cds_to_gff(
    seqid: str,
    f: SeqFeature
) -> List[GFF3Record]:
    from copy import deepcopy
    cdss = []

    f = deepcopy(f)

    id_ = pop_list("ID", f.qualifiers)
    name = pop_list("Name", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)
    target = pop_list("Target", f.qualifiers)
    if target is not None:
        target = Target.parse(target)

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    attributes = GFF3Attributes(
        id=id_,
        name=name,
        alias=f.qualifiers.pop("Alias", None),
        parent=f.qualifiers.pop("Parent", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start: int, end: int, strand: Strand, phase: Phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type="CDS",
            start=start,
            end=end,
            strand=strand,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    if f.location_operator == "join":
        strands = {int_to_strand(fi.strand) for fi in f.location.parts}
        assert len(strands) == 1
        strand = strands.pop()
        del strands
        for part, phase in calculate_phase(f.location.parts, strand):
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, phase)
            cdss.append(cds_part)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.FIRST)
        cdss.append(cds_part)
    return cdss


def gene_to_gff(
    seqid: str,
    f: SeqFeature
) -> List[GFF3Record]:
    from copy import deepcopy
    cdss = []

    f = deepcopy(f)
    id_ = pop_list("ID", f.qualifiers)
    name = pop_list("Name", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)
    target = pop_list("Target", f.qualifiers)
    if target is not None:
        target = Target.parse(target)

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    attributes = GFF3Attributes(
        id=id_,
        name=name,
        alias=f.qualifiers.pop("Alias", None),
        parent=f.qualifiers.pop("Parent", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(
        start: int,
        end: int,
        strand: Strand,
        phase: Phase
    ) -> GFF3Record:
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type=f.type,
            start=start,
            end=end,
            strand=strand,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cdss.append(cds_part)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cdss.append(cds_part)
    return cdss


def candcluster_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str
) -> GFF3Record:
    from copy import deepcopy
    cdss = []

    f = deepcopy(f)

    id_ = pop_list("candidate_cluster_number", f.qualifiers)
    assert id_ is not None
    id_ = str(region) + "-" + str(f.type) + str(id_)

    name = pop_list("Name", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)
    target = pop_list("Target", f.qualifiers)
    if target is not None:
        target = Target.parse(target)

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    f.qualifiers.pop("tool", [])

    attributes = GFF3Attributes(
        id=id_,
        name=name,
        alias=f.qualifiers.pop("Alias", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    attributes.custom["type"] = f.type

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type="region",
            start=start,
            end=end,
            strand=strand,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cdss.append(cds_part)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cdss.append(cds_part)
    return cdss


def region_to_gff(
    seqid: str,
    f: SeqFeature,
    region: Optional[str] = None
) -> GFF3Record:
    from copy import deepcopy
    cdss = []

    f = deepcopy(f)

    if region is None:
        region = pop_list("region_number", f.qualifiers)
        assert region is not None
    id_ = str(region)

    name = pop_list("Name", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)
    target = pop_list("Target", f.qualifiers)
    if target is not None:
        target = Target.parse(target)

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    f.qualifiers.pop("tool", [])

    attributes = GFF3Attributes(
        id=id_,
        name=name,
        alias=f.qualifiers.pop("Alias", None),
        parent=f.qualifiers.pop("Parent", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    attributes.custom["type"] = f.type

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type="region",
            start=start,
            end=end,
            strand=Strand.UNSTRANDED,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cdss.append(cds_part)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cdss.append(cds_part)
    return cdss


def protocluster_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str
) -> GFF3Record:
    from copy import deepcopy
    cdss = []

    f = deepcopy(f)

    id_ = pop_list("protocluster_number", f.qualifiers)
    assert id_ is not None
    id_ = str(region) + "-" + str(f.type) + str(id_)

    name = pop_list("Name", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)
    target = pop_list("Target", f.qualifiers)
    if target is not None:
        target = Target.parse(target)

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    f.qualifiers.pop("tool", [])
    f.qualifiers.pop("core_location", [])

    attributes = GFF3Attributes(
        id=id_,
        name=name,
        alias=f.qualifiers.pop("Alias", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    attributes.custom["type"] = f.type

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type="region",
            start=start,
            end=end,
            strand=Strand.UNSTRANDED,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cdss.append(cds_part)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cdss.append(cds_part)
    return cdss


def proto_core_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str,
    i: int
) -> Tuple[GFF3Record, int]:
    from copy import deepcopy
    cdss = []

    f = deepcopy(f)

    id_ = pop_list("protocluster_number", f.qualifiers)
    assert id_ is not None
    id_ = str(region) + "-protocluster" + str(id_) + "-core" + str(i)
    i += 1

    name = pop_list("Name", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)
    target = pop_list("Target", f.qualifiers)
    if target is not None:
        target = Target.parse(target)

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    f.qualifiers.pop("tool", [])

    attributes = GFF3Attributes(
        id=id_,
        name=name,
        alias=f.qualifiers.pop("Alias", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    attributes.custom["type"] = f.type

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type="region",
            start=start,
            end=end,
            strand=Strand.UNSTRANDED,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cdss.append(cds_part)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cdss.append(cds_part)
    return cdss, i


def create_match(
    parts: Sequence[GFF3Record],
    type_: str,
    region: str,
) -> GFF3Record:
    from copy import deepcopy
    start = min(p.start for p in parts)
    end = max(p.end for p in parts)

    seqid = find_seqid(parts)
    strand = find_strand(parts)
    score = parts[0].score

    source = find_source(parts)

    attributes = deepcopy(parts[0].attributes)
    attributes.id = str(region)

    parent = GFF3Record(
        seqid=seqid,
        source=source,
        type=type_,
        start=start,
        end=end,
        score=score,
        strand=strand,
        phase=Phase.NOT_CDS,
        attributes=attributes
    )

    return parent


def cds_match_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str
) -> GFF3Record:
    from copy import deepcopy
    f = deepcopy(f)

    evalue = pop_list("evalue", f.qualifiers)
    if evalue is not None:
        evalue = float(evalue)

    label = pop_list("label", f.qualifiers)
    gap = pop_list("Gap", f.qualifiers)

    protein_start = pop_list("protein_start", f.qualifiers)
    protein_end = pop_list("protein_end", f.qualifiers)
    if (protein_start is not None) and (protein_end is not None):
        target = Target(label, int(protein_start), int(protein_end))
    else:
        target = None

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    domain_id = pop_list("domain_id", f.qualifiers)
    domain_id = str(region) + "-" + domain_id

    f.qualifiers.pop("tool", [])
    f.qualifiers.pop("locus_tag", [])

    attributes = GFF3Attributes(
        id=domain_id,
        alias=f.qualifiers.pop("Alias", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    if target is None:
        attributes.custom["label"] = label

    attributes.custom["type"] = f.type

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type=f.type,
            start=start,
            end=end,
            strand=strand,
            score=evalue,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    cdss = []
    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cds_part.attributes.id = None
            cds_part.type = "match_part"
            cdss.append(cds_part)

        parent = create_match(cdss, "protein_hmm_match", domain_id)
        parent.children = cdss
        for cds in cdss:
            cds.attributes.parent = [parent.attributes.id]
        cdss.append(parent)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cds_part.type = "protein_hmm_match"
        cdss.append(cds_part)
    return cdss


def pfam_domain_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str
) -> GFF3Record:
    from copy import deepcopy
    f = deepcopy(f)

    evalue = pop_list("evalue", f.qualifiers)
    if evalue is not None:
        evalue = float(evalue)

    f.qualifiers.pop("label", [])
    gap = pop_list("Gap", f.qualifiers)

    aSDomain = pop_list("aSDomain", f.qualifiers)
    protein_start = pop_list("protein_start", f.qualifiers)
    protein_end = pop_list("protein_end", f.qualifiers)
    if (protein_start is not None) and (protein_end is not None):
        target = Target(aSDomain, int(protein_start), int(protein_end))
    else:
        target = None

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    domain_id = pop_list("domain_id", f.qualifiers)
    domain_id = str(region) + "-" + domain_id

    f.qualifiers.pop("tool", [])
    f.qualifiers.pop("locus_tag", [])

    attributes = GFF3Attributes(
        id=domain_id,
        alias=f.qualifiers.pop("Alias", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    attributes.custom["type"] = "PFAM_domain"

    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    if target is None:
        attributes.custom["aSDomain"] = aSDomain

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type=f.type,
            start=start,
            end=end,
            strand=strand,
            score=evalue,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    cdss = []
    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cds_part.attributes.id = None
            cds_part.type = "match_part"
            cdss.append(cds_part)

        parent = create_match(cdss, "protein_hmm_match", domain_id)
        parent.children = cdss
        for cds in cdss:
            cds.attributes.parent = [parent.attributes.id]
        cdss.append(parent)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cds_part.type = "protein_hmm_match"
        cdss.append(cds_part)
    return cdss


def asdomain_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str
) -> GFF3Record:
    from copy import deepcopy
    f = deepcopy(f)

    evalue = pop_list("evalue", f.qualifiers)
    if evalue is not None:
        evalue = float(evalue)

    f.qualifiers.pop("label", [])
    gap = pop_list("Gap", f.qualifiers)

    aSDomain = pop_list("aSDomain", f.qualifiers)
    protein_start = pop_list("protein_start", f.qualifiers)
    protein_end = pop_list("protein_end", f.qualifiers)
    if (protein_start is not None) and (protein_end is not None):
        target = Target(aSDomain, int(protein_start), int(protein_end))
    else:
        target = None

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    domain_id = pop_list("domain_id", f.qualifiers)
    domain_id = str(region) + "-" + domain_id

    f.qualifiers.pop("tool", [])
    f.qualifiers.pop("locus_tag", [])

    attributes = GFF3Attributes(
        id=domain_id,
        alias=f.qualifiers.pop("Alias", None),
        target=target,
        gap=gap,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    attributes.custom["type"] = "aSDomain"
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    if target is None:
        attributes.custom["aSDomain"] = aSDomain

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type=f.type,
            start=start,
            end=end,
            strand=strand,
            score=evalue,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    cdss = []
    if f.location_operator == "join":
        for part in f.location.parts:
            start = part.start.position
            end = part.end.position
            strand = int_to_strand(part.strand)
            cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
            cds_part.attributes.id = None
            cds_part.type = "match_part"
            cdss.append(cds_part)

        parent = create_match(cdss, "protein_hmm_match", domain_id)
        parent.children = cdss
        for cds in cdss:
            cds.attributes.parent = [parent.attributes.id]
        cdss.append(parent)
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cds_part.type = "protein_hmm_match"
        cdss.append(cds_part)
    return cdss


def asmodule_to_gff(
    seqid: str,
    f: SeqFeature,
    region: str
) -> GFF3Record:
    from copy import deepcopy
    f = deepcopy(f)

    complete = pop_list("complete", f.qualifiers)
    complete = complete is not None
    iterative = pop_list("iterative", f.qualifiers)
    iterative = iterative is not None
    starter_module = pop_list("starter_module", f.qualifiers)
    starter_module = starter_module is not None
    final_module = pop_list("final_module", f.qualifiers)
    final_module = final_module is not None

    f.qualifiers.pop("label", [])
    f.qualifiers.pop("tool", [])

    is_circular = pop_list("Is_circular", f.qualifiers)
    is_circular = is_circular == "true"

    contig_edge = pop_list("contig_edge", f.qualifiers)
    contig_edge = contig_edge not in (None, "False")

    attributes = GFF3Attributes(
        alias=f.qualifiers.pop("Alias", None),
        target=None,
        gap=None,
        derives_from=f.qualifiers.pop("Derives_from", None),
        note=f.qualifiers.pop("Note", None),
        dbxref=f.qualifiers.pop("Dbxref", None),
        ontology_term=f.qualifiers.pop("Ontology_term", None),
        is_circular=is_circular,
    )
    attributes.custom = dict(f.qualifiers)
    attributes.custom["type"] = "aSModule"
    if contig_edge:
        attributes.custom["contig_edge"] = "true"

    if "translation" in attributes.custom:
        del attributes.custom["translation"]

    def get_feat(start, end, strand, phase):
        from copy import deepcopy
        return GFF3Record(
            seqid=seqid,
            source="antiSMASH",
            type=f.type,
            start=start,
            end=end,
            strand=strand,
            score=None,
            phase=phase,
            attributes=deepcopy(attributes),
        )

    cdss = []
    if (f.location_operator is not None) and (f.location_operator == "join"):
        raise ValueError(
            "ERROR: ASModule entries should not have compound locations"
        )
    else:
        start = f.location.start.position
        end = f.location.end.position
        strand = int_to_strand(f.location.strand)
        cds_part = get_feat(start, end, strand, Phase.NOT_CDS)
        cds_part.type = "region"
        cdss.append(cds_part)
    return cdss


def run(handle: TextIO):  # noqa: C901
    seq = SeqIO.read(handle, "genbank")
    region_id = handle.name[:-len(".gbk")]

    orig_start = int(
        seq
        .annotations["structured_comment"]
        ["antiSMASH-Data"]
        ["Orig. start"]
    )

    features = deepcopy(seq.features)
    seqid = seq.id
    print(seqid)
    for feature in features:
        feature.location += orig_start

    for f in features:
        if f.type not in (
            "gene",
            "mRNA",
            "CDS",
            "exon",
            "region",
            "cand_cluster",
            "protocluster",
            "proto_core",
            "CDS_motif",
            "PFAM_domain",
            "aSDomain",
            "aSModule",
        ):
            raise ValueError(str(f))

    out = []

    genes = [gene_to_gff(seqid, f)[0] for f in features if f.type == "gene"]
    genes_dict = {g.attributes.id: g for g in genes}
    out.extend(genes)

    mrnas = [gene_to_gff(seqid, f)[0] for f in features if f.type == "mRNA"]
    mrnas_dict = {m.attributes.id: m for m in mrnas}
    out.extend(mrnas)

    for m in mrnas:
        parent = genes_dict[m.attributes.parent[0]]
        parent.children.append(m)
        m.attributes.parent = [parent.attributes.id]

    cdss = [cds_to_gff(seqid, f) for f in features if f.type == "CDS"]
    for c in cdss:
        cd = []
        for ci in c:
            try:
                # TODO make this more general
                parent = mrnas_dict[ci.attributes.id.replace("SNOT", "SNOR")]
                parent.children.append(ci)
                ci.attributes.parent = [parent.attributes.id]
                cd.append(ci)
            except Exception:
                pass
        out.extend(cd)

    exons = [gene_to_gff(seqid, f) for f in features if f.type == "exon"]
    for c in exons:
        cd = []
        for ci in c:
            try:
                parent = mrnas_dict[ci.attributes.parent[0]]
                parent.children.append(ci)
                ci.attributes.parent = [parent.attributes.id]
                cd.append(ci)
            except Exception:
                pass
        out.extend(cd)

    regions = [
        region_to_gff(seqid, f, region_id)[0]
        for f
        in features
        if f.type == "region"
    ]
    assert len(regions) == 1

    region = regions[0]
    out.append(region)

    candclusters = [
        candcluster_to_gff(seqid, f, region_id)
        for f
        in features
        if f.type == "cand_cluster"
    ]

    for c in candclusters:
        out.extend(c)

    protoclusters = [
        protocluster_to_gff(seqid, f, region_id)
        for f
        in features
        if f.type == "protocluster"
    ]

    for c in protoclusters:
        out.extend(c)

    proto_core: List[GFF3Record] = []
    i = 1
    for f in features:
        if f.type != "proto_core":
            continue
        pc, i = proto_core_to_gff(seqid, f, region_id, i)
        proto_core.extend(pc)
        out.extend(pc)

    cds_matches = [
        cds_match_to_gff(seqid, f, region_id)
        for f
        in features
        if f.type == "CDS_motif"
    ]

    for c in cds_matches:
        out.extend(c)

    pfam_domains = [
        pfam_domain_to_gff(seqid, f, region_id)
        for f
        in features
        if f.type == "PFAM_domain"
    ]

    for c in pfam_domains:
        out.extend(c)

    asdomains = [
        asdomain_to_gff(seqid, f, region_id)
        for f
        in features
        if f.type == "aSDomain"
    ]

    for c in asdomains:
        out.extend(c)

    asmodules = [
        asmodule_to_gff(seqid, f, region_id)
        for f
        in features
        if f.type == "aSModule"
    ]

    for c in asmodules:
        out.extend(c)

    return out


def antismash2gff(args: argparse.Namespace) -> None:
    print("##gff-version 3", file=args.outfile)
    out = []
    for f in args.infiles:
        x = run(f)
        out.extend(x)

    seen = set()
    for xi in sorted(out, key=lambda xi: (xi.seqid, xi.start, xi.end)):
        if len(xi.parents) > 0:
            continue

        for child in xi.traverse_children():
            if child in seen:
                continue
            else:
                seen.add(child)

            print(xi, file=args.outfile)
        print("###", file=args.outfile)
    return
