#!/usr/bin/env python3

import logging
from typing import Optional
from typing import List
from typing import Set
from typing import Dict
from typing import Union
from typing import Sequence
from typing import Mapping
from typing import Iterator

from enum import Enum
from collections import defaultdict

from gffpal.attributes import GFFAttributes
from gffpal.higher import fmap, or_else

logger = logging.getLogger(__name__)


TYPE_PARENT_MAP = {
    "mRNA": "gene",
    "CDS": "mRNA",
    "exon": "mRNA",
    "intron": "mRNA",
    "UTR": "mRNA",
    "five_prime_UTR": "mRNA",
    "three_prime_UTR": "mRNA",
    "non_canonical_five_prime_splice_site": "mRNA",
    "non_canonical_three_prime_splice_site": "mRNA",
    "transcription_start_site": "mRNA",
    "transcription_end_site": "mRNA",
    "start_codon": "mRNA",
    "stop_codon": "mRNA",
}


class GFFFormats(Enum):
    GTF2 = 0
    GFF3 = 1


class Strand(Enum):
    PLUS = 0
    MINUS = 1
    UNSTRANDED = 2
    UNKNOWN = 3

    def __str__(self):
        into_str_map: List[str] = ["+", "-", ".", "?"]
        return into_str_map[self.value]

    def __repr__(self):
        return f"Strand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "Strand":
        from_str_map: Dict[str, Strand] = {
            "+": cls.PLUS,
            "-": cls.MINUS,
            ".": cls.UNSTRANDED,
            "?": cls.UNKNOWN,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Phase(Enum):
    FIRST = 0
    SECOND = 1
    THIRD = 2
    NOT_CDS = 3

    def __str__(self):
        into_str_map: List[str] = ["0", "1", "2", "."]
        return into_str_map[self.value]

    def __repr__(self):
        return f"Phase.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "Phase":
        from_str_map: Dict[str, Phase] = {
            "0": cls.FIRST,
            "1": cls.SECOND,
            "2": cls.THIRD,
            ".": cls.NOT_CDS,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class GFFRecord(object):

    columns: List[str] = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    ]

    def __init__(
        self,
        seqid: str,
        source: str,
        type: str,
        start: int,
        end: int,
        score: Optional[float] = None,
        strand: Strand = Strand.UNSTRANDED,
        phase: Phase = Phase.NOT_CDS,
        attributes: Union[GFFAttributes, None] = None,
        parents: List["GFFRecord"] = [],
        children: List["GFFRecord"] = [],
    ) -> None:
        self.seqid = seqid
        self.source = source
        self.type = type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attributes

        self.parents: List[GFFRecord] = []
        self.add_parents(parents)

        self.children: List[GFFRecord] = []
        self.add_children(children)
        return

    def __repr__(self) -> str:
        parameters = []
        for col in self.columns:
            val = repr(getattr(self, col))
            parameters.append(f"{col}={val}")

        joined_parameters = ", ".join(parameters)
        return f"GFFRecord({joined_parameters})"

    def add_child(self, child: "GFFRecord") -> None:
        if child not in self.children:
            self.children.append(child)
        if self not in child.parents:
            child.parents.append(self)
        return

    def add_parent(self, parent: "GFFRecord") -> None:
        if parent not in self.parents:
            self.parents.append(parent)

        if self not in parent.children:
            parent.children.append(self)
        return

    def add_children(self, children: List["GFFRecord"]) -> None:
        for child in children:
            self.add_child(child)
        return

    def add_parents(self, parents: List["GFFRecord"]) -> None:
        for parent in parents:
            self.add_parent(parent)
        return

    @classmethod
    def parse(
        cls,
        string: str,
        format: GFFFormats = GFFFormats.GFF3,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GFFRecord":
        sline = string.strip().split("\t")
        sline_len = len(sline)
        columns_len = len(cls.columns)

        if sline_len < columns_len:
            raise ValueError((
                "Line has too few columns. "
                f"Expected: {columns_len}, Encountered: {sline_len}"
            ))
        elif sline_len > columns_len:
            logger.warning(
                "Line has too many columns. Expected: %s, Encountered: %s",
                columns_len,
                sline_len
            )

        fields = dict(zip(cls.columns, sline))

        # 0-based indexing exclusive
        start = int(fields["start"]) - 1
        end = int(fields["end"])

        if start > end:
            tmp = start
            start = end
            end = tmp
            del tmp

        if fields["score"] in (".", ""):
            score = None
        else:
            score = float(fields["score"])

        strand = Strand.parse(fields["strand"])
        phase = Phase.parse(fields["phase"])

        if format == GFFFormats.GFF3:
            attributes = GFFAttributes.parse(
                fields["attributes"],
                strip_quote=strip_quote,
                unescape=unescape,
            )
        else:
            raise ValueError("Currently only support GFF3 formats.")

        return cls(
            fields["seqid"],
            fields["source"],
            fields["type"],
            start,
            end,
            score,
            strand,
            phase,
            attributes
        )

    def length(self) -> int:
        return self.end - self.start

    def trim_ends(self, length: int) -> None:
        from math import ceil

        if self.length() <= 2:
            length = 0
        elif self.length() < (2 * length):
            length = ceil(self.length() / 4)

        self.start += length
        self.end -= length
        return

    def expand_to_children(self) -> None:
        if len(self.children) == 0:
            return

        min_ = min(c.start for c in self.children)
        max_ = max(c.end for c in self.children)

        if min_ < self.start:
            self.start = min_

        if max_ > self.end:
            self.end = max_
        return

    def shrink_to_children(self) -> None:
        if len(self.children) == 0:
            return

        min_ = min(c.start for c in self.children)
        max_ = max(c.end for c in self.children)

        if min_ > self.start:
            self.start = min_

        if max_ < self.end:
            self.end = max_
        return

    def __str__(self):
        values = []
        for name in self.columns:
            value = getattr(self, name)

            if value is None:
                values.append(".")
            # Convert back to 1-based inclusive indexing.
            elif name == "start":
                values.append(str(value + 1))
            else:
                values.append(str(value))

        return "\t".join(values)


class GFF(object):

    def __init__(self, records: Sequence[GFFRecord]) -> None:
        self.inner: List[GFFRecord] = []
        self.index: Dict[str, List[GFFRecord]] = defaultdict(list)
        self.missing_parents: Dict[str, List[GFFRecord]] = defaultdict(list)
        self.add_records(records)
        return

    def add_record(self, record: GFFRecord) -> None:
        """ """

        # Attributes can be None
        id_ = fmap(lambda a: getattr(a, "id"), record.attributes)

        if id_ is not None:
            self.index[id_].append(record)

        children = self.missing_parents.pop(id_, [])
        record.add_children(children)

        parent_ids: List[str] = or_else(
            [],
            fmap(lambda a: getattr(a, "parent"), record.attributes)
        )

        for parent_id in parent_ids:
            parents = self.index.get(parent_id, [])

            if len(parents) == 0:
                self.missing_parents[parent_id].append(record)
            else:
                record.add_parents(parents)

        self.inner.append(record)
        return

    def add_records(self, records: Sequence[GFFRecord]) -> None:
        for record in records:
            self.add_record(record)
        return

    def get(self, key: str) -> List[GFFRecord]:
        return self.index.get(key, [])

    @classmethod
    def infer_parent(
        cls,
        id: str,
        type: str,
        records: List[GFFRecord]
    ) -> GFFRecord:
        """ """

        if len(records) == 0:
            raise ValueError("Cannot get the parent of an empty set.")

        seqids = list(set(r.seqid for r in records))
        if len(seqids) > 1:
            raise ValueError("Multiple seqids encountered.")
        else:
            seqid = seqids[0]

        start = min(r.start for r in records)
        end = max(r.end for r in records)

        strands = list(set(r.strand for r in records))
        if len(strands) > 1:
            logger.warning("Ambiguous strand found during parent inference.")
            strand = Strand.UNKNOWN
        else:
            strand = strands[0]

        attributes = GFFAttributes(id=id)
        record = GFFRecord(
            seqid=seqid,
            source=".",
            type=type,
            start=start,
            end=end,
            score=None,
            strand=strand,
            phase=Phase.NOT_CDS,
            attributes=attributes,
            children=records,
        )
        return record

    def infer_missing_parents(
        self,
        type_map: Mapping[str, str] = TYPE_PARENT_MAP
    ) -> None:
        parent_ids = list(self.missing_parents.keys())

        for parent_id in parent_ids:
            children = self.missing_parents.pop(parent_id, [])
            if len(children) == 0:
                continue

            type_candidates = list(set(
                type_map[c.type]
                for c
                in children
                if c.type in type_map
            ))
            if len(type_candidates) > 1 or len(type_candidates) == 0:
                logger.warning(
                    "Could not infer parent type for parent %s.",
                    parent_id
                )
                type = "."
            else:
                type = type_candidates[0]

            record = self.infer_parent(parent_id, type, children)
            self.add_record(record)
        return

    def traverse_children(
        self,
        records: Optional[Sequence[GFFRecord]] = None,
        sort: bool = False,
    ) -> Iterator[GFFRecord]:
        """ A depth first traversal of the GFF from parents to children.

        Optionally sort the children by position.
        """

        seen: Set[GFFRecord] = set()
        if records is None:
            to_visit = [f for f in reversed(self.inner) if len(f.parents) == 0]
        else:
            to_visit = [f for f in reversed(records)]

        if sort:
            to_visit.sort(
                key=lambda f: (f.seqid, f.start, f.end, f.type),
                reverse=True,
            )

        while len(to_visit) > 0:
            node = to_visit.pop()

            # NB uses id() for hashing
            if node in seen:
                continue
            else:
                seen.add(node)

            children = list(node.children)
            if sort:
                children.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=True
                )

            to_visit.extend(children)
            yield node

        return None

    def traverse_parents(
        self,
        records: Optional[Sequence[GFFRecord]] = None,
        sort: bool = False,
    ) -> Iterator[GFFRecord]:
        """ A depth first traversal of the GFF from children to parents.

        Optionally sort the parents by position.
        """

        seen: Set[GFFRecord] = set()
        if records is None:
            to_visit = [f for f in self.inner if len(f.children) == 0]
        else:
            to_visit = [f for f in records]  # Quick and dirty copy

        if sort:
            to_visit.sort(key=lambda f: (f.seqid, f.start, f.end, f.type))

        while len(to_visit) > 0:
            node = to_visit.pop()

            if node in seen:
                continue
            else:
                seen.add(node)

            parents = list(node.parents)
            if sort:
                parents.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                )

            to_visit.extend(parents)
            yield node

        return None

    @classmethod
    def parse(
        cls,
        handle: Sequence[str],
        format: GFFFormats = GFFFormats.GFF3,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GFF":
        self = cls([])
        for i, line in enumerate(handle):
            if line.startswith("#"):
                continue
            try:
                record = GFFRecord.parse(
                    line.strip(),
                    strip_quote=strip_quote,
                    unescape=unescape,
                )
            except ValueError as e:
                print(line)
                raise e
            self.add_record(record)

        return self

    def __str__(self) -> str:
        return "\n".join(str(r) for r in self.inner)

    def __getitem__(self, key: int) -> GFFRecord:
        return self.inner[key]

    def __iter__(self) -> Iterator[GFFRecord]:
        return iter(self.inner)

    def select_type(self, type: str) -> Iterator[GFFRecord]:
        for f in self:
            if f.type == type:
                yield f
        return
