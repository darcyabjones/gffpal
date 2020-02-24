from copy import copy

from enum import Enum
from typing import cast
from typing import Optional
from typing import List
from typing import Dict
from typing import Union
from typing import Mapping
from typing import Sequence

from gffpal.higher import fmap

GFF3_KEY_TO_ATTR: Dict[str, str] = {
    "ID": "id",
    "Name": "name",
    "Alias": "alias",
    "Parent": "parent",
    "Target": "target",
    "Gap": "gap",
    "Derives_from": "derives_from",
    "Note": "note",
    "Dbxref": "dbxref",
    "Ontology_term": "ontology_term",
    "Is_circular": "is_circular",
}

GFF3_ATTR_TO_KEY: Dict[str, str] = {v: k for k, v in GFF3_KEY_TO_ATTR.items()}

GFF3_WRITE_ORDER: List[str] = [
    "ID",
    "Name",
    "Alias",
    "Parent",
    "Target",
    "Gap",
    "Derives_from",
    "Note",
    "Dbxref",
    "Ontology_term",
    "Is_circular",
]


GTF_WRITE_ORDER: List[str] = [
    "gene_id",
    "transcript_id",
]


class TargetStrand(Enum):
    PLUS = 0
    MINUS = 1

    def __str__(self) -> str:
        into_str_map = ["+", "-"]
        return into_str_map[self.value]

    def __repr__(self) -> str:
        return "TargetStrand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "TargetStrand":
        from_str_map = {
            "+": cls.PLUS,
            "-": cls.MINUS,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")


class Target(object):

    """
    Indicates the target of a nucleotide-to-nucleotide or
    protein-to-nucleotide alignment.
    The format of the value is "target_id start end [strand]",
    where strand is optional and may be "+" or "-".
    If the target_id contains spaces, they must be escaped as hex escape %20.
    """

    def __init__(
        self,
        target_id: str,
        start: int,
        end: int,
        strand: Optional[TargetStrand] = None,
    ) -> None:
        self.target_id = target_id
        self.start = start
        self.end = end
        self.strand = strand
        return

    def __repr__(self) -> str:
        if self.strand is None:
            return f"Target('{self.target_id}', {self.start}, {self.end})"
        else:
            strand = repr(self.strand)
            return (
                f"Target('{self.target_id}', {self.start}, "
                f"{self.end}, {strand})"
            )

    def __str__(self) -> str:
        target_id = self.target_id.replace(" ", "%20")

        # Recode back to 1 based (inclusive)
        start = self.start + 1

        if self.strand is None:
            return "{} {} {}".format(target_id, start, self.end)
        else:
            return "{} {} {} {}".format(
                target_id,
                start,
                self.end,
                self.strand
            )

    @classmethod
    def parse(cls, string: str) -> "Target":
        split_string = string.strip().split(" ")

        if len(split_string) < 3:
            raise ValueError("Too few fields in Target string.")

        elif len(split_string) == 3:
            target_id = split_string[0]
            start = int(split_string[1]) - 1  # We want 0 based exclusive
            end = int(split_string[2])
            return cls(target_id, start, end)

        elif len(split_string) == 4:
            target_id = split_string[0]
            start = int(split_string[1]) - 1  # We want 0 based exclusive
            end = int(split_string[2])
            strand = TargetStrand.parse(split_string[3])
            return cls(target_id, start, end, strand)

        else:
            raise ValueError("Too many fields in Target string.")


class GapCode(Enum):
    MATCH = 0
    INSERT = 1
    DELETE = 2
    FORWARD = 3
    REVERSE = 4

    def __str__(self) -> str:
        into_str_map = ["M", "I", "D", "F", "R"]
        return into_str_map[self.value]

    def __repr__(self) -> str:
        return f"GapCode.{self.name}"

    @classmethod
    def parse(cls, string: str) -> "GapCode":
        from_str_map = {
            "M": cls.MATCH,
            "I": cls.INSERT,
            "D": cls.DELETE,
            "F": cls.FORWARD,
            "R": cls.REVERSE,
        }

        try:
            return from_str_map[string]
        except KeyError:
            valid = list(from_str_map.keys())
            raise ValueError(f"Invalid option. Must be one of {valid}")
        return


class GapElement(object):

    def __init__(self, code: GapCode, length: int) -> None:
        """ """
        assert length >= 0

        self.code = code
        self.length = length
        return

    def __str__(self) -> str:
        return "{}{}".format(self.code, self.length)

    def __repr__(self) -> str:
        code = repr(self.code)
        return f"GapElement({code}, {self.length})"

    @classmethod
    def parse(cls, string: str) -> "GapElement":
        string = string.strip()
        code = GapCode.parse(string[:1])
        length = int(string[1:])
        return cls(code, length)


class Gap(object):

    def __init__(self, elements: Sequence[GapElement]) -> None:
        self.elements = list(elements)
        return

    def __str__(self) -> str:
        return " ".join(map(str, self.elements))

    def __repr__(self) -> str:
        elems = repr(list(self.elements))
        return f"Gap({elems})"

    @classmethod
    def parse(cls, string: str) -> "Gap":
        split_string = string.strip().split(" ")
        elements = [GapElement.parse(s) for s in split_string if s != '']
        return cls(elements)


def parse_attr_list(string: str) -> List[str]:
    return list(f.strip() for f in string.strip(", ").split(","))


def parse_bool(string: str) -> bool:
    if string in ("true", "True", "TRUE"):
        return True
    elif string in ("false", "False", "FALSE"):
        return False
    else:
        raise ValueError(f"Invalid boolean encountered: {string}")


class Attributes(object):

    def __init__(
        self,
        custom: Mapping[str, str] = {},
    ) -> None:
        self.custom = dict(custom)
        return

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "Attributes":
        raise NotImplementedError("Attributes is a baseclass.")
        return

    def get(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        raise NotImplementedError("Attributes is a baseclass.")
        return

    def pop(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        raise NotImplementedError("Attributes is a baseclass.")
        return


class GTFAttributes(Attributes):

    def __init__(
        self,
        gene_id: Optional[str] = None,
        transcript_id: Optional[str] = None,
        custom: Mapping[str, str] = {},
    ) -> None:
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        super().__init__(custom)
        return

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GTFAttributes":
        if string.strip() in (".", ""):
            return cls()

        seen: Dict[str, int] = dict()

        fields = (
            f.strip().split("=", maxsplit=1)
            for f in
            string.strip(" ;").split(";")
        )

        kvpairs: Dict[str, str] = dict()

        for k, v in fields:
            k = k.strip()
            v = v.strip()

            if k in seen:
                seen[k] += 1
                # If it's a duplicate entry, add a number to it.
                k = f"{k}_{seen[k]}"
            else:
                seen[k] = 1

            if strip_quote:
                v = v.strip(" '\"")

            kvpairs[k] = v

        for k, count in seen.items():
            if count > 1:
                v = kvpairs.pop(k)
                kvpairs[f"{k}_1"] = v

        gene_id = kvpairs.pop("gene_id", None)
        transcript_id = kvpairs.pop("transcript_id", None)
        return cls(gene_id, transcript_id, kvpairs)

    def __str__(self) -> str:
        keys = []
        keys.extend(GTF_WRITE_ORDER)
        keys.extend(self.custom.keys())

        kvpairs = []
        for key in keys:
            value = self[key]
            if value is None:
                continue

            # Deal with newlines etc.
            value = value.encode("unicode_escape").decode("utf-8")
            kvpairs.append((key, value))

        return " ; ".join(f"{k} {v}" for k, v in kvpairs)

    def __repr__(self) -> str:
        param_names = copy(GTF_WRITE_ORDER)
        param_names.append("custom")

        parameters = []
        for param in param_names:
            value = getattr(self, param)

            if value is None or value == []:
                continue

            if isinstance(value, list):
                value = repr(list(value))
            else:
                value = repr(value)

            parameters.append(f"{param}={value}")

        joined_parameters = ", ".join(parameters)
        return f"GTFAttributes({joined_parameters})"


    def __getitem__(
        self,
        key: str
    ) -> Optional[str]:
        if key in GTF_WRITE_ORDER:
            return getattr(self, key)
        else:
            return self.custom[key]

    def __setitem__(
        self,
        key: str,
        value: Optional[str],
    ) -> None:
        if key in GTF_WRITE_ORDER:
            setattr(self, key, value)
        elif value is not None:
            self.custom[key] = value
        return

    def __delitem__(self, key: str) -> None:
        """ Removes an item from the custom dictionary or resets
        attribute to default """

        if key in GTF_WRITE_ORDER:
            setattr(self, key, None)

        else:
            del self.custom[key]

        return

    def get(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Gets an atrribute or element from the custom dict. """

        if key in GTF_WRITE_ORDER:
            return getattr(self, key)
        else:
            return self.custom.get(key, default)

    def pop(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Removes an item from the attributes and returns its value.

        If the item is one of the reserved GFF3 terms, the
        value is reset to the default.
        """

        if key in GTF_WRITE_ORDER:
            value = getattr(self, key)
            del self[key]
            return value

        else:
            return self.custom.pop(key, default)


class GFF3Attributes(Attributes):

    def __init__(
        self,
        id: Optional[str] = None,
        name: Optional[str] = None,
        alias: Sequence[str] = [],
        parent: Sequence[str] = [],
        target: Optional[Target] = None,
        gap: Optional[Gap] = None,
        derives_from: Sequence[str] = [],
        note: Sequence[str] = [],
        dbxref: Sequence[str] = [],
        ontology_term: Sequence[str] = [],
        is_circular: Optional[bool] = None,
        custom: Mapping[str, str] = {},
    ) -> None:
        self.id = id
        self.name = name
        self.alias = alias
        self.parent = parent
        self.target = target
        self.gap = gap
        self.derives_from = derives_from
        self.note = note
        self.dbxref = dbxref
        self.ontology_term = ontology_term
        self.is_circular = is_circular

        super().__init__(custom=dict(custom))
        return

    @classmethod
    def _parse_list_of_strings(
        cls,
        value: str,
        strip_quote: bool = False,
        unescape: bool = False
    ) -> List[str]:
        """ Parses a gff attribute list of strings. """
        if value == "":
            return []

        if strip_quote:
            value = value.strip("\"' ")

        if unescape:
            return [cls._attr_unescape(v) for v in parse_attr_list(value)]
        else:
            return parse_attr_list(value)

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> "GFF3Attributes":
        if string.strip() in (".", ""):
            return cls()

        fields = (
            f.split("=", maxsplit=1)
            for f
            in string.strip(" ;").split(";")
        )

        if strip_quote:
            kvpairs: Dict[str, str] = {
                k.strip(): v.strip(" '\"")
                for k, v
                in fields
            }
        else:
            kvpairs = {
                k.strip(): v.strip()
                for k, v
                in fields
            }

        id = kvpairs.pop("ID", None)
        name = kvpairs.pop("Name", None)

        alias = cls._parse_list_of_strings(
            kvpairs.pop("Alias", ""),
            strip_quote,
            unescape
        )

        parent = cls._parse_list_of_strings(
            kvpairs.pop("Parent", ""),
            strip_quote,
            unescape
        )

        target: Optional[Target] = fmap(
            Target.parse,
            fmap(
                lambda x: cls._attr_unescape(x) if unescape else x,
                kvpairs.pop("Target", None)
            )
        )

        gap = fmap(
            Gap.parse,
            fmap(
                lambda x: cls._attr_unescape(x) if unescape else x,
                kvpairs.pop("Gap", None)
            )
        )

        derives_from = cls._parse_list_of_strings(
            kvpairs.pop("Derives_from", ""),
            strip_quote,
            unescape
        )

        note = cls._parse_list_of_strings(
            kvpairs.pop("Note", ""),
            strip_quote,
            unescape
        )

        dbxref = cls._parse_list_of_strings(
            kvpairs.pop("Dbxref", ""),
            strip_quote,
            unescape
        )

        ontology_term = cls._parse_list_of_strings(
            kvpairs.pop("Ontology_term", ""),
            strip_quote,
            unescape
        )

        is_circular = fmap(parse_bool, kvpairs.pop("Is_circular", None))

        return cls(
            id,
            name,
            alias,
            parent,
            target,
            gap,
            derives_from,
            note,
            dbxref,
            ontology_term,
            is_circular,
            kvpairs
        )

    def __str__(self) -> str:
        keys = []
        keys.extend(GFF3_WRITE_ORDER)
        keys.extend(self.custom.keys())

        kvpairs = []
        for key in keys:
            value = self[key]
            if value is None or value == []:
                continue

            key = self._attr_escape(key)

            if isinstance(value, list):
                value = ", ".join(self._attr_escape(str(v)) for v in value)

            elif isinstance(value, bool):
                value = "true" if value else "false"

            else:
                value = self._attr_escape(str(value))

            kvpairs.append((key, value))

        return ";".join(f"{k}={v}" for k, v in kvpairs)

    def __repr__(self) -> str:
        param_names = [GFF3_KEY_TO_ATTR[k] for k in GFF3_WRITE_ORDER]
        param_names.append("custom")

        parameters = []
        for param in param_names:
            value = getattr(self, param)

            if value is None or value == []:
                continue

            if isinstance(value, list):
                value = repr(list(value))
            else:
                value = repr(value)

            parameters.append(f"{param}={value}")

        joined_parameters = ", ".join(parameters)
        return f"GFF3Attributes({joined_parameters})"

    @staticmethod
    def _attr_escape(string: str) -> str:
        return (
            string
            .replace(";", "%3B")
            .replace(",", "%2C")
            .replace("=", "%3D")
            .replace("\t", "%09")
            .replace("\n", "%0A")
            .replace("\r", "%0D")
        )

    @staticmethod
    def _attr_unescape(string: str) -> str:
        return (
            string
            .replace("%3B", ";")
            .replace("%2C", ",")
            .replace("%3D", "=")
            .replace("%09", "\t")
            .replace("%0A", "\n")
            .replace("%0D", "\r")
        )

    def __getitem__(
        self,
        key: str,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        if key in GFF3_KEY_TO_ATTR:
            return getattr(self, GFF3_KEY_TO_ATTR[key])
        else:
            return self.custom[key]

    def __setitem__(
        self,
        key: str,
        value: Union[str, Sequence[str], Target, Gap, bool, None],
    ) -> None:
        """ If the key is an attr set it, otherwise update the custom dict."""

        if key in GFF3_KEY_TO_ATTR:
            setattr(self, GFF3_KEY_TO_ATTR[key], value)
        else:
            self.custom[key] = cast(str, value)  # Cast is for mypy
        return

    def __delitem__(self, key: str) -> None:
        """ Removes an item from the custom dictionary or resets
        attribute to default """

        if key in ("ID", "Name", "Target", "Gap"):
            setattr(self, GFF3_KEY_TO_ATTR[key], None)

        elif key in ("Alias", "Parent", "Derives_from", "Note",
                     "Dbxref", "Ontology_term"):
            setattr(self, GFF3_KEY_TO_ATTR[key], [])

        elif key == "Is_circular":
            setattr(self, GFF3_KEY_TO_ATTR[key], False)

        else:
            del self.custom[key]

        return

    def get(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Gets an atrribute or element from the custom dict. """

        if key in GFF3_KEY_TO_ATTR:
            return getattr(self, GFF3_KEY_TO_ATTR[key])
        else:
            return self.custom.get(key, default)

    def pop(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None]:
        """ Removes an item from the attributes and returns its value.

        If the item is one of the reserved GFF3 terms, the
        value is reset to the default.
        """

        if key in GFF3_KEY_TO_ATTR:
            value = getattr(self, GFF3_KEY_TO_ATTR[key])
            del self[GFF3_KEY_TO_ATTR[key]]
            return value

        else:
            return self.custom.pop(key, default)
