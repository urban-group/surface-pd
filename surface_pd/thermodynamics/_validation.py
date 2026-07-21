"""Internal validation helpers for thermodynamic value objects."""

from collections.abc import Mapping
from numbers import Integral, Real
from types import MappingProxyType

import numpy as np


def validate_variable_name(value: object) -> str:
    """Return a valid state-variable name or raise ``ValueError``."""
    if not isinstance(value, str) or not value or not value.isidentifier():
        raise ValueError(
            "state variable name must be a nonempty identifier-like string"
        )
    return value


def validate_finite_real(value: object, field_name: str) -> float:
    """Return *value* as a finite float or raise a descriptive exception."""
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{field_name} must be a finite real number")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{field_name} must be a finite real number")
    return result


def validate_positive_integer(value: object, field_name: str) -> int:
    """Return *value* as a positive integer or raise an exception."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{field_name} must be a positive integer")
    result = int(value)
    if result <= 0:
        raise ValueError(f"{field_name} must be a positive integer")
    return result


def validate_identifier(value: object, field_name: str) -> str:
    """Return a nonempty, single-line identifier without ``:``."""
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a string")
    result = value.strip()
    if not result or ":" in result or "\n" in result or "\r" in result:
        raise ValueError(
            f"{field_name} must be a nonempty single-line string without ':'"
        )
    return result


def validate_provenance(value: object, field_name: str) -> str:
    """Return nonempty, single-line provenance text."""
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a string")
    result = value.strip()
    if not result or "\n" in result or "\r" in result:
        raise ValueError(f"{field_name} must be nonempty single-line text")
    return result


def validate_composition(value: object) -> MappingProxyType:
    """Return an immutable copy of an absolute elemental composition."""
    if not isinstance(value, Mapping):
        raise TypeError("composition must be a mapping of component counts")
    if not value:
        raise ValueError("composition must not be empty")

    result = {}
    for component, count in value.items():
        try:
            valid_component = validate_variable_name(component)
        except ValueError as error:
            raise ValueError(
                "composition component names must be identifier-like strings"
            ) from error
        if isinstance(count, bool) or not isinstance(count, Real):
            raise TypeError("composition counts must be nonnegative integers")
        numeric_count = float(count)
        if not np.isfinite(numeric_count) or not numeric_count.is_integer():
            raise ValueError("composition counts must be nonnegative integers")
        integer_count = int(numeric_count)
        if integer_count < 0:
            raise ValueError("composition counts must be nonnegative integers")
        result[valid_component] = integer_count

    if not any(result.values()):
        raise ValueError(
            "composition must contain at least one positive count"
        )
    return MappingProxyType(result)
