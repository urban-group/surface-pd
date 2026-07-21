"""Immutable named thermodynamic state."""

from collections.abc import Iterator, Mapping
from types import MappingProxyType
from typing import Any

import numpy as np

from ._validation import validate_variable_name


class ThermodynamicState(Mapping[str, np.ndarray]):
    """Store finite, broadcast-compatible thermodynamic variables.

    Parameters
    ----------
    values : mapping of str to array-like
        Named scalar or array-valued state variables. Names must be nonempty
        Python-style identifiers. Values are converted to floating NumPy
        arrays, defensively copied, broadcast to one common shape, and made
        read-only. Every value must contain only finite real numbers.

    Raises
    ------
    TypeError
        If ``values`` is not a mapping or a value is not real numeric data.
    ValueError
        If a name is invalid, a value is non-finite, or value shapes cannot be
        broadcast together.

    Notes
    -----
    An empty state is valid and has scalar shape ``()``. The class is an
    immutable mapping; arrays returned through indexing are also read-only.
    """

    __slots__ = ("_shape", "_values")

    def __init__(self, values: Mapping[str, Any]):
        if not isinstance(values, Mapping):
            raise TypeError("values must be a mapping of state variables")

        converted = {}
        for name, value in values.items():
            valid_name = validate_variable_name(name)
            array = np.asarray(value)
            if (
                np.issubdtype(array.dtype, np.bool_)
                or np.issubdtype(array.dtype, np.complexfloating)
                or not np.issubdtype(array.dtype, np.number)
            ):
                raise TypeError(
                    f"{valid_name} must contain only finite real numbers"
                )
            array = np.array(array, dtype=float, copy=True)
            if not np.isfinite(array).all():
                raise ValueError(
                    f"{valid_name} must contain only finite real numbers"
                )
            converted[valid_name] = array

        try:
            shape = np.broadcast_shapes(
                *(array.shape for array in converted.values())
            )
        except ValueError as error:
            raise ValueError(
                "thermodynamic state values must be broadcast-compatible"
            ) from error

        broadcast_values = {}
        for name, array in converted.items():
            broadcast = np.broadcast_to(array, shape)
            broadcast.setflags(write=False)
            broadcast_values[name] = broadcast

        self._shape = shape
        self._values = MappingProxyType(broadcast_values)

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the common NumPy broadcast shape of all state variables."""
        return self._shape

    def __getitem__(self, key: str) -> np.ndarray:
        """Return one read-only, broadcast state array."""
        return self._values[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate over state-variable names in input order."""
        return iter(self._values)

    def __len__(self) -> int:
        """Return the number of named state variables."""
        return len(self._values)
