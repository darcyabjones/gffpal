from enum import Enum
from typing import Optional
from typing import TypeVar
from typing import List
from typing import Dict
from typing import Union
from typing import Mapping
from typing import Sequence

TargetT = TypeVar("TargetT", bound="Target")
TargetStrandT = TypeVar("TargetStrandT", bound="TargetStrand")
T = TypeVar("T")

_KEY_TO_ATTR = {
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

_ATTR_TO_KEY = {v: k for k, v in _KEY_TO_ATTR.items()}

_WRITE_ORDER = [
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


class TargetStrand(Enum):
    PLUS = 0
    MINUS = 1

    def __str__(self):
        into_str_map = ["+", "-"]
        return into_str_map[self.value]

    def __repr__(self):
        return "TargetStrand.{self.name}"

    @classmethod
    def parse(cls, string: str) -> TargetStrandT:
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
        strand: Optional[TargetStrandT] = None,
    ):
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
    def parse(cls, string: str) -> TargetT:
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
    def parse(cls, string: str):
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

    def __str__(self):
        return "{}{}".format(self.code, self.length)

    def __repr__(self):
        code = repr(self.code)
        return f"GapElement({code}, {self.length})"

    @classmethod
    def parse(cls, string: str):
        string = string.strip()
        code = GapCode.parse(string[:1])
        length = int(string[1:])
        return cls(code, length)


class Gap(object):

    def __init__(self, elements: Sequence[GapElement]):
        self.elements = list(elements)
        return

    def __str__(self):
        return " ".join(map(str, self.elements))

    def __repr__(self):
        elems = repr(list(self.elements))
        return f"Gap({elems})"

    @classmethod
    def parse(cls, string: str):
        split_string = string.strip().split(" ")
        elements = (GapElement.parse(s) for s in split_string if s != '')
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


class GFFAttributes(object):

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
        self.custom = dict(custom)
        return


    @classmethod
    def _parse_elem(
        cls,
        key: str,
        value: str,
        strip_quote: bool = False,
        unescape: bool = False,
    ) -> str:
        if strip_quote:
            value = value.strip("\"' ")

        if key in ("Alias", "Parent", "Derives_from",
                   "Note", "Dbxref", "Ontology_term"):
            if unescape:
                return [cls._attr_unescape(v) for v in parse_attr_list(value)]
            else:
                return parse_attr_list(value)

        if unescape:
            value = cls._attr_unescape(value)

        if key == "Target":
            return Target.parse(value)
        elif key == "Gap":
            return Gap.parse(value)
        elif key == "Is_circular":
            return parse_bool(value)
        else:
            return str(value)

    @classmethod
    def parse(
        cls,
        string: str,
        strip_quote: bool = False,
        unescape: bool = False,
    ):
        fields = (
            f.split("=", maxsplit=1)
            for f
            in string.strip(" ;").split(";")
        )

        kvpairs = {
            k.strip(): v.strip()
            for k, v
            in fields
        }

        kwargs = {}

        for key, attr in _KEY_TO_ATTR.items():
            value = kvpairs.pop(key, None)
            if value is not None:
                kwargs[attr] = cls._parse_elem(
                    key,
                    value,
                    strip_quote=strip_quote,
                    unescape=unescape
                )

        return cls(**kwargs, custom=kvpairs)

    def strip_quotes(self):

        return

    def __str__(self):
        keys = []
        keys.extend(_WRITE_ORDER)
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
        param_names = [_KEY_TO_ATTR[k] for k in _WRITE_ORDER]
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

        parameters = ", ".join(parameters)
        return f"GFFAttributes({parameters})"

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
        if key in _KEY_TO_ATTR:
            return getattr(self, _KEY_TO_ATTR[key])
        else:
            return self.custom[key]

    def __setitem__(
        self,
        key: str,
        value: Union[str, Sequence[str], Target, Gap, bool, None],
    ) -> None:
        """ If the key is an attr set it, otherwise update the custom dict."""

        if key in _KEY_TO_ATTR:
            setattr(self, _KEY_TO_ATTR[key], value)
        else:
            self.custom[key] = value
        return

    def __delitem__(self, key: str) -> None:
        """ Removes an item from the custom dictionary or resets
        attribute to default """

        if key in ("ID", "Name", "Target", "Gap"):
            setattr(self, _KEY_TO_ATTR[key], None)

        elif key in ("Alias", "Parent", "Derives_from", "Note",
                     "Dbxref", "Ontology_term"):
            setattr(self, _KEY_TO_ATTR[key], [])

        elif key == "Is_circular":
            setattr(self, _KEY_TO_ATTR[key], False)

        else:
            del self.custom[key]

        return

    def get(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None, T]:
        """ Gets an atrribute or element from the custom dict. """

        if key in _KEY_TO_ATTR:
            return getattr(self, _KEY_TO_ATTR[key])
        else:
            return self.custom.get(key, default)

    def pop(
        self,
        key: str,
        default: Union[str, Sequence[str], Target, Gap, bool, None] = None,
    ) -> Union[str, Sequence[str], Target, Gap, bool, None, T]:
        """ Removes an item from the attributes and returns its value.

        If the item is one of the reserved GFF3 terms, the
        value is reset to the default.
        """

        if key in _KEY_TO_ATTR:
            value = getattr(self, _KEY_TO_ATTR[key])
            del self[_KEY_TO_ATTR[key]]
            return value

        else:
            return self.custom.pop(key, default)
