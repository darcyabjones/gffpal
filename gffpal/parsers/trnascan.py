import re
import logging
from typing import TextIO
from typing import Optional
from typing import Sequence, Iterator, List, Dict

from gffpal.parsers.parsers import ParseError, LineParseError
from gffpal.parsers.parsers import MULTISPACE_REGEX
from gffpal.parsers.parsers import (
    parse_int,
    parse_float,
    parse_string_not_empty
)

logger = logging.getLogger(__name__)


class TRNAScanRecord(object):

    column_order = [
        "seqid",
        "num",
        "start",
        "end",
        "trna_type",
        "anticodon",
        "intron_starts",
        "intron_ends",
        "infernal_score",
        "note",
    ]

    def __init__(
        self,
        seqid: str,
        start: int,
        end: int,
        trna_type: str,
        anticodon: str,
        num: Optional[int],
        intron_starts: Sequence[int],
        intron_ends: Sequence[int],
        infernal_score: Optional[float],
        note: Optional[str],
    ) -> None:
        self.seqid = seqid
        self.start = start
        self.end = end
        self.trna_type = trna_type
        self.anticodon = anticodon
        self.num = num
        self.intron_starts = intron_starts
        self.intron_ends = intron_ends
        self.infernal_score = infernal_score
        self.note = note
        return

    @classmethod
    def from_line(cls, line: str) -> "TRNAScanRecord":
        sline = MULTISPACE_REGEX.split(line.rstrip())

        if ((len(sline) != len(cls.column_order)) and
                (len(sline) != len(cls.column_order) - 1)):
            raise LineParseError(
                "Line had the wrong number of columns. "
                f"Expected {len(cls.column_order)} or "
                f"{len(cls.column_order) - 1} but got {len(sline)}."
            )

        record: Dict[str, str] = {
            k.strip(): v.strip()
            for k, v
            in zip(cls.column_order, sline)
        }

        start = parse_int(record["start"], "start")
        end = parse_int(record["end"], "end")
        num = parse_int(record["num"], "num")

        infernal_score = parse_float(
            record["infernal_score"],
            "infernal_score"
        )

        if record["intron_starts"] == "0" and record["intron_ends"] == "0":
            intron_starts: List[int] = []
            intron_ends: List[int] = []
        else:
            intron_starts = [
                parse_int(i.strip(), "intron_starts")
                for i
                in record["intron_starts"].split(",")
            ]

            intron_ends = [
                parse_int(i.strip(), "intron_ends")
                for i
                in record["intron_ends"].split(",")
            ]

        return cls(
            seqid=parse_string_not_empty(record["seqid"], "seqid"),
            start=start,
            end=end,
            trna_type=parse_string_not_empty(record["trna_type"], "trna_type"),
            anticodon=parse_string_not_empty(record["anticodon"], "anticodon"),
            num=num,
            intron_starts=intron_starts,
            intron_ends=intron_ends,
            infernal_score=infernal_score,
            note=record.get("note", None),
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["TRNAScanRecord"]:

        started = False
        for i, line in enumerate(handle, 1):
            if line.startswith("-"):
                started = True
                continue

            elif not started:
                continue

            try:
                yield cls.from_line(line)
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


class TRNAScanSS(object):

    def __init__(
        self,
        seqid: str,
        start: int,
        end: int,
        trna_type: str,
        anticodon: str,
        anticodon_start: int,
        anticodon_end: int,
        score: float,
        seq: str,
        ss: str,
        num: Optional[int],
    ) -> None:
        self.seqid = seqid
        self.start = int(start)
        self.end = int(end)
        self.trna_type = trna_type
        self.anticodon = anticodon
        self.anticodon_start = anticodon_start
        self.anticodon_end = anticodon_end
        self.score = score
        self.seq = seq
        self.ss = ss
        self.num = num
        return

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["TRNAScanSS"]:
        record = {}
        id_regex = re.compile(r"(.+)\.trna(\d+)\s+\((\d+)-(\d+)\)")
        type_regex = re.compile((
            r"Type: (\S+)\s+Anticodon: "
            r"(\S+).*\((\d+)-(\d+).*Score: (\d+\.?\d*)"
        ))
        for line in handle:
            line = line.strip()

            if line == "" or line.startswith("*"):
                continue
            elif line.startswith("Type"):
                match = type_regex.search(line)

                if match is None:
                    raise ValueError(f"Could not parse type line: {line}")

                (trna_type, anticodon, anticodon_start,
                 anticodon_end, score) = match.groups()
                record["trna_type"] = trna_type
                record["anticodon"] = anticodon
                record["anticodon_start"] = anticodon_start
                record["anticodon_end"] = anticodon_end
                record["score"] = score

            elif line.startswith("Possible"):
                continue

            elif line.startswith("Seq"):
                seq = line[5:].strip()
                record["seq"] = seq
            elif line.startswith("Str"):
                ss = line[5:].strip()
                record["ss"] = ss
            else:
                match = id_regex.search(line)
                if match is None:
                    logger.warning("Encountered unexpected line: %s", line)
                    continue
                elif "seqid" in record:
                    # Its not the first record
                    yield cls(
                        record["seqid"],
                        int(record["start"]),
                        int(record["end"]),
                        record["trna_type"],
                        record["anticodon"],
                        int(record["anticodon_start"]),
                        int(record["anticodon_end"]),
                        float(record["score"]),
                        record["seq"],
                        record["ss"],
                        int(record["num"])
                    )

                seqid, num, start, end = match.groups()
                record = {
                    "seqid": seqid,
                    "num": num,
                    "start": start,
                    "end": end,
                }

        if "seqid" in record:
            yield cls(
                record["seqid"],
                int(record["start"]),
                int(record["end"]),
                record["trna_type"],
                record["anticodon"],
                int(record["anticodon_start"]),
                int(record["anticodon_end"]),
                float(record["score"]),
                record["seq"],
                record["ss"],
                int(record["num"])
            )
        return
