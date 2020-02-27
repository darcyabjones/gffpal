from typing import Optional
from typing import TextIO
from typing import Iterator

from gffpal.parsers.parsers import (
    parse_int,
    parse_float,
    parse_string_not_empty,
    MULTISPACE_REGEX,
    LineParseError,
    ParseError
)


class DomTbl(object):

    """ """

    def __init__(
        self,
        target_name: str,
        target_acc: Optional[str],
        target_len: int,
        query_name: str,
        query_acc: Optional[str],
        query_len: int,
        full_evalue: float,
        full_score: float,
        full_bias: float,
        match_num: int,
        nmatches: int,
        domain_c_evalue: float,
        domain_i_evalue: float,
        domain_score: float,
        domain_bias: float,
        hmm_from: int,
        hmm_to: int,
        ali_from: int,
        ali_to: int,
        env_from: int,
        env_to: int,
        acc: float,
        description: Optional[str]
    ) -> None:
        self.target_name = target_name
        self.target_acc = target_acc
        self.target_len = target_len
        self.query_name = query_name
        self.query_acc = query_acc
        self.query_len = query_len
        self.full_evalue = full_evalue
        self.full_score = full_score
        self.full_bias = full_bias
        self.match_num = match_num
        self.nmatches = nmatches
        self.domain_c_evalue = domain_c_evalue
        self.domain_i_evalue = domain_i_evalue
        self.domain_score = domain_score
        self.domain_bias = domain_bias
        self.hmm_from = hmm_from
        self.hmm_to = hmm_to
        self.ali_from = ali_from
        self.ali_to = ali_to
        self.env_from = env_from
        self.env_to = env_to
        self.acc = acc
        self.description = description
        return

    @classmethod
    def from_line(cls, line: str) -> "DomTbl":
        if line == "":
            raise LineParseError("The line was empty.")

        sline = MULTISPACE_REGEX.split(line.strip(), maxsplit=22)
        if len(sline) != 22 and len(sline) != 23:
            # Technically because of the max_split this should be impossible.
            # the description line is allowed to have spaces.
            print(sline)
            raise LineParseError(
                "The line had the wrong number of columns. "
                f"Expected 22 or 23 but got {len(sline)}"
            )

        if sline[22] == "-" or sline[22] == "":
            description: Optional[str] = None
        else:
            description = sline[22]

        return cls(
            parse_string_not_empty(sline[0], "target_name"),
            parse_string_not_empty(sline[1], "target_acc"),
            parse_int(sline[2], "target_len"),
            parse_string_not_empty(sline[3], "query_name"),
            parse_string_not_empty(sline[4], "query_acc"),
            parse_int(sline[5], "query_len"),
            parse_float(sline[6], "full_evalue"),
            parse_float(sline[7], "full_score"),
            parse_float(sline[8], "full_bias"),
            parse_int(sline[9], "match_num"),
            parse_int(sline[10], "nmatches"),
            parse_float(sline[11], "domain_c_evalue"),
            parse_float(sline[12], "domain_i_evalue"),
            parse_float(sline[13], "domain_score"),
            parse_float(sline[14], "domain_bias"),
            parse_int(sline[15], "hmm_from"),
            parse_int(sline[16], "hmm_to"),
            parse_int(sline[17], "ali_from"),
            parse_int(sline[18], "ali_to"),
            parse_int(sline[19], "env_from"),
            parse_int(sline[20], "env_to"),
            parse_float(sline[21], "acc"),
            description
        )

    @classmethod
    def from_file(cls, handle: TextIO) -> Iterator["DomTbl"]:
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
