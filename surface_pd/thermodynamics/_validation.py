"""Internal validation helpers for thermodynamic value objects."""

from numbers import Integral, Real

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
