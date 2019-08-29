import re
import logging
from typing import Optional
from typing import Sequence, Iterator

logger = logging.getLogger(__name__)

COLUMN_ORDER = [
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


class TRNAScanRecord(object):

    def __init__(
        self,
        seqid: str,
        start: int,
        end: int,
        trna_type: str,
        anticodon: str,
        num: Optional[int] = None,
        intron_starts: Sequence[int] = [],
        intron_ends: Sequence[int] = [],
        infernal_score: Optional[float] = None,
        note: Optional[str] = None,
    ):
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
    def parse(cls, handle: Sequence[str]) -> Iterator["TRNAScanRecord"]:
        regex = re.compile(r"\s+")

        started = False
        for line in handle:
            if line.startswith("-"):
                started = True
                continue

            elif not started:
                continue

            sline = regex.split(line.rstrip())
            if len(COLUMN_ORDER) < len(sline):
                logger.warning("Line has unexpected number of columns.")
                logger.warning("offending line is: %s", line)

            record = {
                k.strip(): v.strip()
                for k, v
                in zip(COLUMN_ORDER, sline)
            }

            record["start"] = int(record["start"])
            record["end"] = int(record["end"])
            record["num"] = int(record["num"])
            record["infernal_score"] = float(record["infernal_score"])

            if record["intron_starts"] == "0" and record["intron_ends"] == "0":
                del record["intron_starts"]
                del record["intron_ends"]
            else:
                record["intron_starts"] = [
                    int(i.strip())
                    for i
                    in record["intron_starts"].split(",")
                ]

                record["intron_ends"] = [
                    int(i.strip())
                    for i
                    in record["intron_ends"].split(",")
                ]

            yield cls(**record)

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
        num: Optional[int] = None,
    ):
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
    def parse(cls, handle: Sequence[str]) -> Iterator["TRNAScanSS"]:
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
                record["anticodon_start"] = int(anticodon_start)
                record["anticodon_end"] = int(anticodon_end)
                record["score"] = float(score)

            elif line.startswith("Possible intron"):
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
                    yield cls(**record)

                seqid, num, start, end = match.groups()
                record = {
                    "seqid": seqid,
                    "num": int(num),
                    "start": int(start),
                    "end": int(end),
                }

        if "seqid" in record:
            yield cls(**record)
        return
