#!/usr/bin/env python3

import logging
from typing import cast
from typing import Optional, Union
from typing import Set, List, Dict, Tuple
from typing import Sequence, Mapping
from typing import Iterator
from typing import Hashable, Any
from typing import TypeVar, Type
from typing import Generic

from enum import Enum
from collections import defaultdict
from collections import deque

from gffpal.parsers.parsers import (
    parse_int,
    parse_float,
    parse_or_none,
    is_one_of,
    parse_string_not_empty
)

from gffpal.attributes import Attributes, GFF3Attributes
from gffpal.higher import fmap, or_else

logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)


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


T = TypeVar('T')


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

    @classmethod
    def infer_from_many(cls, strands: Sequence["Strand"]) -> "Strand":
        """ Infer the strand from a sequence of strands. """

        if len(set(strands)) > 1:
            logger.warning("Ambiguous strand found during parent inference.")
            return Strand.UNKNOWN
        else:
            return list(strands)[0]


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


HashT = TypeVar("HashT", Hashable, Any)
AttrT = TypeVar("AttrT", bound=Attributes)


class GFFRecord(Generic[AttrT]):

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
        attributes: Optional[AttrT] = None,
        parents: Optional[Sequence["GFFRecord"]] = None,
        children: Optional[Sequence["GFFRecord"]] = None,
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

        self.parents: List[GFFRecord[AttrT]] = []
        if parents is not None:
            self.add_parents(parents)

        self.children: List[GFFRecord[AttrT]] = []
        if children is not None:
            self.add_children(children)
        return

    def __str__(self) -> str:
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

    def __repr__(self) -> str:
        parameters = []
        for col in self.columns:
            val = repr(getattr(self, col))
            parameters.append(f"{val}")

        joined_parameters = ", ".join(parameters)
        return f"GFFRecord({joined_parameters})"

    def length(self) -> int:
        """ Returns the distance between start and end. """
        return self.end - self.start

    def __len__(self) -> int:
        return self.length()

    def add_child(self, child: "GFFRecord[AttrT]") -> None:
        if child not in self.children:
            self.children.append(child)
        if self not in child.parents:
            child.parents.append(self)
        return

    def add_parent(self, parent: "GFFRecord[AttrT]") -> None:
        if parent not in self.parents:
            self.parents.append(parent)

        if self not in parent.children:
            parent.children.append(self)
        return

    def add_children(self, children: Sequence["GFFRecord[AttrT]"]) -> None:
        for child in children:
            self.add_child(child)
        return

    def add_parents(self, parents: Sequence["GFFRecord[AttrT]"]) -> None:
        for parent in parents:
            self.add_parent(parent)
        return

    def traverse_children(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFFRecord[AttrT]"]:
        """ A graph traversal of this `GFFRecord`s children.

        Keyword arguments:
        sort -- Sort the children so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFFRecord[AttrT]] = set()

        to_visit = deque([self])

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
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
                    reverse=should_reverse
                )

            to_visit.extend(children)
            yield node

        return None

    def traverse_parents(
        self,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator["GFFRecord[AttrT]"]:
        """ A graph traversal of this `GFFRecord`s parents.

        Keyword arguments:
        sort -- Sort the parents so that they appear in increasing order.
        breadth -- Perform a breadth first search rather than the default
            depth first search.

        Returns:
        A generator yielding `GFFRecord`s.
        """

        should_reverse = not breadth

        seen: Set[GFFRecord[AttrT]] = set()
        to_visit = deque([self])

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            if node in seen:
                continue
            else:
                seen.add(node)

            parents = list(node.parents)
            if sort:
                parents.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(parents)
            yield node

        return None

    @classmethod
    def parse(
        cls,
        string: str,
        attr: Type[AttrT] = cast(Type[AttrT], GFF3Attributes),
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GFFRecord[AttrT]":
        """ Parse a gff line string as a `GFFRecord`.

        Keyword arguments:
        string -- The gff line to parse.
        format -- What format the gff file is in.
            Currently only GFF3 is supported.
        strip_quote -- Strip quotes from attributes values. The specification
            says that they should not be stripped, so we don't by default.
        unescape -- Unescape reserved characters in the attributes to their
            original values. I.E. some commas, semicolons, newlines etc.

        Returns:
        A `GFFRecord`
        """

        sline = string.strip().split("\t")
        sline_len = len(sline)
        columns_len = len(cls.columns)

        if sline_len == columns_len - 1:
            logger.warning(
                ("Line has has too few columns columns. "
                 "Probably it is missing the attributes"),
            )
        elif sline_len < columns_len:
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

        fields: Dict[str, str] = dict(zip(cls.columns, sline))
        if sline_len == columns_len - 1:
            fields["attributes"] = ""

        # 0-based indexing exclusive
        start = parse_int(fields["start"], "start") - 1
        end = parse_int(fields["end"], "end")

        if start > end:
            tmp = start
            start = end
            end = tmp
            del tmp

        score = parse_or_none(fields["score"], "score", ".", parse_float)

        strand = Strand.parse(is_one_of(
            fields["strand"],
            ["-", "+", ".", "?"],
            "strand"
        ))

        phase = Phase.parse(is_one_of(
            fields["phase"],
            ["0", "1", "2", "."],
            "phase"
        ))

        attributes = cast(
            AttrT,
            attr.parse(
                fields["attributes"],
                strip_quote=strip_quote,
                unescape=unescape,
            )
        )

        return cls(
            parse_string_not_empty(fields["seqid"], "seqid"),
            parse_string_not_empty(fields["source"], "source"),
            parse_string_not_empty(fields["type"], "type"),
            start,
            end,
            score,
            strand,
            phase,
            attributes
        )

    def trim_ends(self, length: int) -> None:
        """ Trim a specified number of bases from each end of the feature. """

        from math import ceil

        if self.length() <= 2:
            length = 0
        elif self.length() < (2 * length):
            length = ceil(self.length() / 4)

        self.start += length
        self.end -= length
        return

    def pad_ends(self, length: int, max_value: Optional[int] = None) -> None:

        if (self.start - length) < 0:
            self.start = 0
        else:
            self.start -= length

        if max_value is not None and (self.end + length) > max_value:
            self.end = max_value
        else:
            self.end += length
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


class GFF3Record(GFFRecord[GFF3Attributes]):

    @classmethod
    def infer_from_children(
        cls,
        children: Sequence["GFF3Record"],
        id: Optional[str] = None,
        seqid: Optional[str] = None,
        source: str = ".",
        type: str = ".",
        strand: Optional[Strand] = None,
        score: Optional[float] = None,
        phase: Phase = Phase.NOT_CDS,
    ) -> "GFF3Record":
        """ """

        if len(children) == 0:
            raise ValueError("Cannot get the parent of an empty set.")

        if seqid is None:
            seqid = cls._infer_seqid_from_children(children)

        if id is None:
            id = cls._infer_id_from_children(children)

        if strand is None:
            strand = Strand.infer_from_many([f.strand for f in children])

        start = min(r.start for r in children)
        end = max(r.end for r in children)

        attributes = GFF3Attributes(id=id)
        record = cast(
            GFF3Record,
            GFFRecord(
                seqid,
                source,
                type,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
                children=children,
            )
        )
        return record

    @staticmethod
    def _infer_id_from_children(children: Sequence["GFF3Record"]) -> str:
        """ Try to predict the id to use based on the children's parents. """

        ids: Dict[str, int] = defaultdict(int)
        for child in children:
            if child.attributes is None:
                continue

            for parent_id in set(child.attributes.parent):
                ids[parent_id] += 1

        correct_num = [k for k, v in ids.items() if v == len(children)]
        if len(correct_num) == 1:
            id = correct_num[0]
        else:
            raise ValueError("Cannot infer parent id.")

        return id

    @staticmethod
    def _infer_seqid_from_children(children: Sequence["GFF3Record"]) -> str:
        seqids = set(r.seqid for r in children)
        if len(seqids) > 1:
            raise ValueError("Multiple seqids encountered.")
        else:
            return list(seqids)[0]

    def copy(self) -> "GFF3Record":
        """ You'll still need to update the ID """
        from copy import copy, deepcopy

        node_copy = copy(self)
        node_copy.attributes = deepcopy(self.attributes)
        node_copy.attributes.parent = []
        node_copy.children = []
        node_copy.parents = []
        return node_copy

def make_root_id(
    template: Optional[str],
    id_: Optional[str],
    type_: str,
    global_index: int,
) -> str:
    if type_ == ".":
        type_ = "notype"

    if template is None:
        return None

    if "{id}" not in tem
    return

    def break_bubbles(
        self,
        id_template: str = "{id}_{index}",
        singleton_id_template: Optional[str] = None,
        root_id_template: Optional[str] = None,
        global_index_start: int = 1
    ) -> "GFF3Record":
        from collections import deque

        feature_map: Dict["GFF3Record", List["GFF3Record"]] = defaultdict(list)
        to_visit = deque([self])
        global_indices: Dict[str, int] = defaultdict(
            lambda: global_index_start
        )

        while len(to_visit) > 0:

            # Doing breadth first traversal.
            node = to_visit.popleft()
            if node in feature_map:
                continue

            to_visit.extend(node.children)

            if (len(node.parents) == 0) or (node == self):
                node_copy = node.copy()
                feature_map[node].append(node_copy)
                if ((node.attributes is not None) and
                        (node.attributes.id is not None) and
                        (root_id_template )):
                    new_id = root_id_template.format(
                        id=node.attributes.id,
                        type=node.type,
                        global_index=global_indices[node.type],
                    )
                    node_copy.attributes.id = new_id

                global_indices[node.type] += 1
                continue

            elif len(node.parents) == 1:
                template_id = "{id}"

            elif len(node.parents) > 1:
                template_id = "{id}_{index}"

            flattened_parents = flatten_list_of_lists(
                feature_map[p]
                for p
                in node.parents
            )

            for local_index, parent in enumerate(flattened_parents, 1):
                node_copy = node.copy()
                node_copy.attributes.id = template_id.format(
                    id=node.attributes.id,
                    index=local_index,
                    parent_id=parent.attributes.id
                )

                node_copy.add_parent(parent)
                node_copy.attributes.parent.append(parent.attributes.id)
                feature_map[node].append(node_copy)

        assert len(feature_map[self]) == 1
        return feature_map[self][0]


class GFF(object):

    def __init__(self, records: Sequence[GFF3Record]) -> None:
        self.inner: List[GFF3Record] = []
        self.index: Dict[str, List[GFF3Record]] = defaultdict(list)
        self.missing_parents: Dict[str, List[GFF3Record]] = defaultdict(list)
        self.add_records(records)
        return

    def __str__(self) -> str:
        return "\n".join(str(r) for r in self.inner)

    def __getitem__(
        self,
        key: Union[int, slice, List[int]]
    ) -> Union[GFF3Record, "GFF"]:
        if isinstance(key, int):
            return self.inner[key]
        elif isinstance(key, slice):
            return self.__class__(self.inner[key])
        elif isinstance(key, list):
            return self.__class__([self.inner[k] for k in key])
        else:
            raise ValueError(
                "GFF __getitem__ can only take an int, slice or list of ints."
            )

    def __iter__(self) -> Iterator[GFF3Record]:
        return iter(self.inner)

    def __len__(self) -> int:
        return len(self.inner)

    def get(self, key: Union[str, List[str]]) -> List[GFF3Record]:
        if isinstance(key, str):
            return self.index.get(key, [])
        elif isinstance(key, list):
            out = list()
            for k in key:
                out.extend(self.index.get(k, []))

            return out
        else:
            raise ValueError("GFF.get can only take a str or list of str.")

    def add_record(self, record: GFF3Record) -> None:
        """ """

        # Attributes can be None
        id_: Optional[str] = fmap(
            lambda a: getattr(a, "id"),
            record.attributes
        )

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

    def add_records(self, records: Sequence[GFF3Record]) -> None:
        for record in records:
            self.add_record(record)
        return

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

            if len(type_candidates) != 1:
                logger.warning(
                    "Could not infer parent type for parent %s.",
                    parent_id
                )
                type = "."
            else:
                type = type_candidates[0]

            record = GFF3Record.infer_from_children(
                children,
                id=parent_id,
                type=type,
            )
            self.add_record(record)
        return

    def select_type(self, type: str) -> Iterator[GFF3Record]:
        for f in self:
            if f.type == type:
                yield f
        return

    def groupby(
        self,
        keys: List[HashT]
    ) -> Iterator[Tuple[HashT, "GFF"]]:
        """ """

        from collections import defaultdict
        out: Dict[HashT, List[GFF3Record]] = defaultdict(list)

        for key, row in zip(keys, self.inner):
            out[key].append(row)

        for key, features in out.items():
            yield key, GFF(features)

        return

    def traverse_children(
        self,
        records: Optional[Sequence[GFF3Record]] = None,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator[GFF3Record]:
        """ A graph traversal of the GFF from parents to children.

        Optionally sort the children by position.
        """

        should_reverse = not breadth

        seen: Set[GFF3Record] = set()

        if records is None:
            init_nodes = [
                f for
                f in
                self.inner
                if len(f.parents) == 0
            ]
        else:
            init_nodes = [f for f in records]

        if sort:
            init_nodes.sort(
                key=lambda f: (f.seqid, f.start, f.end, f.type),
                reverse=should_reverse,
            )
        elif should_reverse:
            init_nodes = init_nodes[::-1]

        to_visit = deque(init_nodes)

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            # NB uses id() for hashing
            if node in seen:
                continue
            else:
                seen.add(node)

            children = cast(List[GFF3Record], list(node.children))
            if sort:
                children.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(children)
            yield node

        return None

    def traverse_parents(
        self,
        records: Optional[Sequence[GFF3Record]] = None,
        sort: bool = False,
        breadth: bool = False,
    ) -> Iterator[GFF3Record]:
        """ A graph traversal of the GFF from children to parents.

        Optionally sort the parents by position.
        """

        should_reverse = not breadth

        seen: Set[GFF3Record] = set()
        if records is None:
            init_nodes = [f for f in self.inner if len(f.children) == 0]
        else:
            init_nodes = [f for f in records]  # Quick and dirty copy

        if sort:
            init_nodes.sort(
                key=lambda f: (f.seqid, f.start, f.end, f.type),
                reverse=should_reverse
            )
        elif should_reverse:
            init_nodes = init_nodes[::-1]

        to_visit = deque(init_nodes)

        while len(to_visit) > 0:

            if breadth:
                node = to_visit.popleft()
            else:
                node = to_visit.pop()

            if node in seen:
                continue
            else:
                seen.add(node)

            parents = cast(List[GFF3Record], list(node.parents))
            if sort:
                parents.sort(
                    key=lambda f: (f.seqid, f.start, f.end, f.type),
                    reverse=should_reverse
                )

            to_visit.extend(parents)
            yield node

        return None

    @classmethod
    def parse(
        cls,
        handle: Sequence[str],
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GFF":
        self = cls([])
        for i, line in enumerate(handle):
            if line.startswith("#"):
                continue
            try:
                record = GFF3Record.parse(
                    line.strip(),
                    strip_quote=strip_quote,
                    unescape=unescape,
                )
            except ValueError as e:
                print(line)
                raise e
            self.add_record(cast(GFF3Record, record))

        return self


def flatten_list_of_lists(li: Sequence[Sequence[T]]) -> Iterator[T]:
    for i in li:
        for j in i:
            yield j
    return
