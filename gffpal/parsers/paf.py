from typing import NamedTuple
from typing import Iterator, List
from typing import TextIO
from typing import Optional

from gffpal.parsers.parsers import (
    parse_int,
    is_one_of,
    parse_string_not_empty,
    LineParseError,
    ParseError
)


class PAF(NamedTuple):

    query: str
    qlen: int
    qstart: int
    qend: int
    strand: str
    target: str
    tlen: int
    tstart: int
    tend: int
    nmatch: int
    alilen: int
    mq: int
    attrs: List[str]

    @staticmethod
    def columns() -> List[str]:
        return [
            "query",
            "qlen",
            "qstart",
            "qend",
            "strand",
            "target",
            "tlen",
            "tstart",
            "tend",
            "nmatch",
            "alilen",
            "mq",
        ]

    @classmethod
    def from_line(cls, line: str) -> 'PAF':
        sline = line.strip().split("\t")
        if len(sline) < len(cls.columns()):
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected at least {len(cls.columns())} but got {len(sline)}."
            )

        dline = dict(zip(cls.columns(), sline))
        attrs = sline[len(cls.columns()):]

        return cls(
            parse_string_not_empty(dline["query"], "query"),
            parse_int(dline["qlen"], "qlen"),
            parse_int(dline["qstart"], "qstart"),
            parse_int(dline["qend"], "qend"),
            is_one_of(dline["strand"], ["+", "-"], "strand"),
            parse_string_not_empty(dline["target"], "target"),
            parse_int(dline["tlen"], "tlen"),
            parse_int(dline["tstart"], "tstart"),
            parse_int(dline["tend"], "tend"),
            parse_int(dline["nmatch"], "nmatch"),
            parse_int(dline["alilen"], "alilen"),
            parse_int(dline["mq"], "mq"),
            attrs
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator['PAF']:
        for i, line in enumerate(handle):
            sline = line.strip()

            if sline.startswith("#"):
                continue
            elif sline == "":
                continue

            try:
                yield cls.from_line(sline)

            except LineParseError as e:
                if hasattr(handle, "name"):
                    filename: Optional[str] = handle.name
                else:
                    filename = None

                raise ParseError(
                    filename,
                    i,
                    e.message
                )
        return

    def __str__(self) -> str:
        line = [str(getattr(self, c)) for c in self.columns()]
        line.extend(self.attrs)
        return "\t".join(line)
