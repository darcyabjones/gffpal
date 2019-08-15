""" Some higher order functions for dealing with static analysis. """

from typing import TypeVar
from typing import Optional
from typing import Callable
from typing import Union


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


def or_else(default: U, option: Optional[T]) -> Union[T, U]:
    """ Replaces None with some default value. """
    if option is None:
        return default
    else:
        return option
