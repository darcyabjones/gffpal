from typing import NamedTuple
from typing import TextIO
from typing import Iterator
from typing import Optional

from intervaltree import Interval

from gffpal.gff import GFF3Record, Strand
from gffpal.attributes import GFF3Attributes, Target
from gffpal.parsers.parsers import (
    parse_int,
    parse_float,
    parse_string_not_empty,
    LineParseError,
    ParseError
)


class Coords(NamedTuple):

    rstart: int
    rend: int
    qstart: int
    qend: int
    strand: Strand
    ralnlen: int
    qalnlen: int
    pid: float
    rlen: int
    qlen: int
    rcov: float
    qcov: float
    ref: str
    query: str

    @classmethod
    def from_line(cls, line: str) -> "Coords":
        if line.strip() == "":
            raise LineParseError("The line was empty")

        sline = line.strip().split("\t")

        if len(sline) < 13:
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected at least 13 but got {len(sline)}."
            )

        rstart = parse_int(sline[0], "reference_start")
        rend = parse_int(sline[1], "reference_end")

        # This shouldn't ever happen AFAIK
        assert rstart <= rend, line
        rstart -= 1

        qstart = parse_int(sline[2], "query_start")
        qend = parse_int(sline[3], "query_end")

        if qstart > qend:
            strand = Strand.MINUS
            temp = qstart
            qstart = qend - 1
            qend = temp

            del temp
        else:
            strand = Strand.PLUS
            qstart -= 1

        # This isn't as type safe as I'd like
        return cls(
            rstart,
            rend,
            qstart,
            qend,
            strand,
            parse_int(sline[4], "reference_alnlen"),
            parse_int(sline[5], "query_alnlen"),
            parse_float(sline[6], "pid"),
            parse_int(sline[7], "reference_len"),
            parse_int(sline[8], "query_len"),
            parse_float(sline[9], "reference_cov"),
            parse_float(sline[10], "query_cov"),
            parse_string_not_empty(sline[11], "reference"),
            parse_string_not_empty(sline[12], "query"),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["Coords"]:
        for i, line in enumerate(handle):

            # Skip the headers.
            if i < 4:
                continue

            sline = line.strip()
            if sline == "":
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

        if self.strand == Strand.MINUS:
            qstart = self.qend
            qend = self.qstart + 1
        else:
            qstart = self.qstart + 1
            qend = self.qend

        cols = [str(self.rstart + 1), str(self.rend), str(qstart), str(qend)]

        for attr in [
            "ralnlen",
            "qalnlen",
            "pid",
            "rlen",
            "qlen",
            "rcov",
            "qcov",
            "ref",
            "query"
        ]:
            cols.append(str(getattr(self, attr)))

        return "\t".join(cols)

    def as_interval(self) -> Interval:
        return Interval(self.rstart, self.rend, self)

    def as_gffrecord(
        self,
        source="MUMmer",
        type="nucleotide_match"
    ) -> GFF3Record:
        return GFF3Record(
            self.ref,
            source,
            type,
            self.rstart,
            self.rend,
            score=self.qcov,
            strand=self.strand,
            attributes=GFF3Attributes(
                target=Target(self.query, self.qstart, self.qend),
                custom={
                    "pid": str(self.pid),
                    "contig_id": str(self.query),
                    "contig_coverage": str(self.qcov),
                    "contig_length": str(self.qlen),
                    "contig_alignment_length": str(self.qalnlen),
                    "scaffold_alignment_length": str(self.ralnlen),
                }
            )
        )
