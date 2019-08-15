""" Some higher order functions for dealing with static analysis. """

from typing import TypeVar
from typing import Optional
from typing import Callable


T = TypeVar("T")
U = TypeVar("U")


def fmap(function: Callable[[T], U], option: Optional[T]) -> Optional[U]:
    if option is None:
        return None
    else:
        return function(option)


def applicative(
    function: Callable[[T], Optional[U]],
    option: Optional[T],
) -> Optional[U]:
    """ Same as fmap except for the type signature. """
    if option is None:
        return None
    else:
        return function(option)
    return


def or_else(default: T, option: Optional[T]) -> T:
    if option is None:
        return default
    else:
        return option
