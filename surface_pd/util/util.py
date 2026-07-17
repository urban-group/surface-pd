"""
Utility functions for surface-pd operations.

This module provides various helper functions for scaling matrices, data
conversions, validation, and dictionary manipulations used throughout the
surface-pd package.
"""

import ast
import copy
import re

from surface_pd.error import NonIntegerError, TooLargeSlabError

_INT_PATTERN = re.compile(r"^[+-]?\d+$")
_FLOAT_PATTERN = re.compile(
    r"^[+-]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+[eE][+-]?\d+)|"
    r"(?:(?:\d+\.\d*)|(?:\.\d+))[eE][+-]?\d+)$"
)


def define_scaling_matrix(a: float, b: float, multiple: int) -> list[int]:
    """
    Generate the scaling matrix for creating supercells.

    This function determines the appropriate scaling matrix based on the ratio
    of lattice parameters a and b, and the desired supercell size. The scaling
    is constrained to a maximum size of 4 to prevent excessive computation.

    Args:
        a: Lattice parameter a in Angstroms.
        b: Lattice parameter b in Angstroms.
        multiple: Maximum target cell size (supercell multiplier). Must be ≤ 4.

    Returns
    -------
        Scaling matrix [scale_a, scale_b, scale_c] where scale_c is always 1
        (no scaling in the slab direction).

    Raises
    ------
        TooLargeSlabError: If multiple > 4, as this would generate too many
            structures.
    """
    if multiple > 4:
        raise TooLargeSlabError
    if multiple == 4:
        scaling_matrix = [2, 2, 1]
    else:
        if round(a / b, 0) == 2:
            scaling_matrix = [1, multiple, 1]
        elif round(b / a, 0) == 2:
            scaling_matrix = [multiple, 1, 1]
        else:
            scaling_matrix = [multiple, 1, 1]
    return scaling_matrix


def _parse_csv_scalar(value: str):
    """
    Parse a scalar value from ``csv2dict`` input.

    Accepted values are integers, floats, quoted strings, and unquoted string
    tokens. Other Python expressions are rejected instead of evaluated.
    """
    value = value.strip()
    if not value:
        raise ValueError("Expected non-empty value")

    if _INT_PATTERN.fullmatch(value):
        return int(value)
    if _FLOAT_PATTERN.fullmatch(value):
        return float(value)
    if value[0] in {"'", '"'}:
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError) as exc:
            raise ValueError(f"Invalid quoted string: {value}") from exc
        if not isinstance(parsed, str):
            raise ValueError(f"Unsupported literal value: {value}")
        return parsed

    try:
        expression = ast.parse(value, mode="eval")
    except SyntaxError:
        if value[0].isdigit() or value[0] in {"+", "-", "."}:
            raise ValueError(
                f"Unsupported expression in csv2dict value: {value}"
            ) from None
        return value
    if isinstance(expression.body, ast.Name):
        return value
    raise ValueError(f"Unsupported expression in csv2dict value: {value}")


def csv2dict(csvlist: list[str]) -> dict:
    """
    Convert list of comma-separated values to a structured dictionary.

    The hierarchy of the dictionary is expressed by equal signs ("=") and
    colons (":"). This is useful for parsing user input for species
    substitutions and compositions. Scalar values are parsed as integers,
    floats, quoted strings, or unquoted strings. Expressions such as
    ``1+1`` are rejected instead of evaluated.

    Example:
        Input: ["Co=Co:0.5&Ni:0.5", "Li=Na"]
        Output: {"Co": {"Co": 0.5, "Ni": 0.5}, "Li": "Na"}

    Args:
        csvlist: List of comma-separated values with hierarchy defined by
            "=" (key-value) and ":" (nested dictionary) or "&" (separator
            for nested entries).

    Returns
    -------
        Structured dictionary with deterministic scalar conversion where
        applicable.

    Raises
    ------
        ValueError: If an item is malformed or contains an unsupported
            expression.
    """
    outdict = {}
    for item in csvlist:
        if "=" not in item:
            raise ValueError(f"Expected KEY=VALUE entry: {item}")
        key, value = (part.strip() for part in item.split("=", 1))
        if not key or not value:
            raise ValueError(f"Expected KEY=VALUE entry: {item}")
        if ":" in value:
            nested_value = {}
            for nested_item in value.split("&"):
                if ":" not in nested_item:
                    raise ValueError(
                        f"Expected nested KEY:VALUE entry: {nested_item}"
                    )
                nested_key, nested_scalar = (
                    part.strip() for part in nested_item.split(":", 1)
                )
                if not nested_key or not nested_scalar:
                    raise ValueError(
                        f"Expected nested KEY:VALUE entry: {nested_item}"
                    )
                nested_value[nested_key] = _parse_csv_scalar(nested_scalar)
            value = nested_value
        else:
            value = _parse_csv_scalar(value)
        outdict[key] = value
    return outdict


def check_int(num: float) -> float:
    """
    Validate that a number is effectively an integer within tolerance.

    This function checks if a floating-point number is close enough to an
    integer (within 1e-2) and returns the integer value if valid.

    Args:
        num: Number to check (typically a float from calculations).

    Returns
    -------
        The rounded integer value as a float if within tolerance.

    Raises
    ------
        NonIntegerError: If the number differs from the nearest integer by
            more than 0.01.
    """
    eps = 1e-2
    int_f = float(round(num))
    if abs(num - int_f) < eps:
        return int_f
    else:
        raise NonIntegerError


def get_values_nested_dict(m: dict):
    """
    Recursively extract all values from a nested dictionary.

    This generator function traverses through nested dictionaries and yields
    all leaf values (non-dictionary values).

    Args:
        m: Dictionary that may contain nested dictionaries.

    Yields
    ------
        All non-dictionary values found in the nested structure.
    """
    for val in m.values():
        if isinstance(val, dict):
            yield from get_values_nested_dict(val)
        else:
            yield val


def all_int(n: list) -> bool:
    """
    Check if every item in the list is an integer.

    Args:
        n: List of numbers to check (must have is_integer() method).

    Returns
    -------
        True if all items are integers, False otherwise.
    """
    m = [x.is_integer() for x in n]
    if all(m):
        return True
    else:
        return False


def have_zero(n: list) -> bool:
    """
    Check if any item in the list is zero.

    Args:
        n: List of numbers to check.

    Returns
    -------
        True if any item equals zero, False otherwise.
    """
    return any(x == 0 for x in n)


def replace_dummy(subs_dict: dict, dummy_species: list) -> dict:
    """
    Replace real species keys with dummy species in substitution dictionary.

    This function replaces actual chemical species with dummy/placeholder
    species in the substitution dictionary. This is useful for enumeration
    algorithms that need to track specific sites without chemical meaning.
    Entries with zero occupancy are removed from the output.

    Args:
        subs_dict: Species and occupancy dictionaries containing the species
            mapping in string-string pairs. Format:
            {species: {species: occupancy}}.
        dummy_species: List of dummy species names (e.g., ["X", "Y"]) to
            replace real species as dictionary keys.

    Returns
    -------
        Updated dictionary with dummy species as keys and real species as
        nested values. Entries with zero occupancy are removed.
    """
    for i, (key, value) in enumerate(list(subs_dict.items())):
        subs_dict[dummy_species[i]] = subs_dict.pop(key)
        value[dummy_species[i]] = value.pop(key)

    updated_dict = copy.deepcopy(subs_dict)
    for key in subs_dict.keys():
        if subs_dict[key][key] == 0:
            updated_dict.pop(key)

    return updated_dict
