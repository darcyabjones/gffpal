#!/usr/bin/env python3

import re
from typing import Optional
from typing import Sequence
from typing import Callable
from typing import TypeVar


MULTISPACE_REGEX = re.compile(r"\s+")
T = TypeVar("T")


class ParseError(Exception):
    """ Some aspect of parsing failed. """

    def __init__(
        self,
        filename: Optional[str],
        line: Optional[int],
        message: str
    ):
        self.filename = filename
        self.line = line
        self.message = message
        return


class LineParseError(Exception):

    def __init__(self, message: str):
        self.message = message
        return


def parse_or_none(
    field: str,
    field_name: str,
    none_value: str,
    fn: Callable[[str, str], T],
) -> Optional[T]:
    """ If the value is the same as the none value, will return None.
    Otherwise will attempt to run the fn with field and field name as the
    first and 2nd arguments.
    """

    if field == none_value:
        return None

    try:
        val = fn(field, field_name)
    except LineParseError as e:
        msg = e.message + (
            f"\nThe value may also be '{none_value}', which will be"
            "interpreted as None."
        )
        raise LineParseError(msg)

    return val


def parse_int(field: str, field_name: str) -> int:
    """ Parse a string as an integer, raising a custom error if fails. """

    try:
        return int(field)
    except ValueError:
        raise LineParseError(
            f"Could not parse value in {field_name} column as an integer. "
            f"The offending value was: '{field}'."
        )


def parse_float(field: str, field_name: str) -> float:
    """ Parse a string as a float, raising a custom error if fails. """

    try:
        return float(field)
    except ValueError:
        raise LineParseError(
            f"Could not parse value in {field_name} column as a float. "
            f"The offending value was: '{field}'."
        )


def parse_bool(
    field: str,
    field_name: str,
    true_value: str,
    false_value: str
) -> bool:
    """ """

    if field == true_value:
        return True
    elif field == false_value:
        return False
    else:
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            f"Must be either {true_value} or {false_value}."
        )


def is_value(field: str, field_name: str, value: str) -> bool:
    """ """
    return field == value


def parse_string_not_empty(field: str, field_name: str) -> str:
    """ """

    if field.strip() == "":
        raise LineParseError(f"The value in column: '{field_name}' was empty.")
    else:
        return field


def is_one_of(field: str, options: Sequence[str], field_name: str) -> str:
    """ """

    if field in options:
        return field
    else:
        raise LineParseError(
            f"Invalid value: '{field}' in the column: '{field_name}'. "
            f"Must be one of {options}."
        )
