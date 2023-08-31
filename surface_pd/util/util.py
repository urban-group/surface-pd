import copy

import pandas as pd

from surface_pd.error import TooLargeSlabError, NonIntegerError


def define_scaling_matrix(a, b, multiple):
    """
    Generate the scaling matrix based on the different values of the
    lattice parameters a and b.

    Args:
        a: Lattice parameter a
        b: Lattice parameter a
        multiple: Maximum target cell size

    Returns:
        Scaling matrix that will be used to create the supercell.
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


def csv2dict(csvlist: list):
    """
     Convert list of comma separated values to a dictionary.

    The hierarchy of the dictionary is expressed by equal signs ("=") and
    colons (":").  For example

       ["Co=Co:0.5,Ni:0.5", "Li=Na"]

    will be converted to

       {"Co": {"Co": 0.5, "Ni": 0.5}, "Li": "Na"}

    Args:
        csvlist: list of comma separated values.

    Returns:
        dictionary
    """
    def trynumeric(v):
        try:
            out = int(v)
        except ValueError:
            try:
                out = pd.eval(v)
            except ValueError:
                out = v
        return out

    outdict = {}
    for item in csvlist:
        key, value = item.split("=")
        if ":" in value:
            value = {s.strip(): trynumeric(v.strip()) for s, v in
                     [el.split(":") for el in value.split("&")]}
        else:
            value = trynumeric(value.strip())
        outdict[key.strip()] = value
    return outdict


def check_int(num):
    """
    Check if the number is an integer.

    Args:
        num:

    Returns:

    """
    eps = 1e-2
    int_f = float(round(num))
    if abs(num - int_f) < eps:
        return int_f
    else:
        raise NonIntegerError


def get_values_nested_dict(m: dict):
    """
    Get the values inside the nested dictionary.

    Args:
        m: Dict

    Returns:
        values

    """
    for val in m.values():
        if isinstance(val, dict):
            yield from get_values_nested_dict(val)
        else:
            yield val


def all_int(n: list):
    """
    Check if every item in the list is an integer.

    Args:
        n: Input list

    """
    m = [x.is_integer() for x in n]
    if all(m):
        return True
    else:
        return False


def have_zero(n: list):
    """
    Check if any item in the list is zero.

    Args:
        n: Input list

    """
    return any([x == 0 for x in n])


def replace_dummy(subs_dict: dict,
                  dummy_species: list):
    """
    Replace the keys (species) in the "replace" dict with "dummy" species.

    Args:
        subs_dict: Species and occupancy dictionaries containing the species
            mapping in string-string pairs.
        dummy_species: A special specie for representing non-traditional
            elements or species.

    Returns:
        Updated dict with "dummy" species as the keys

    """
    for i, key, value in zip(range(len(dummy_species)),
                             subs_dict.keys(),
                             subs_dict.values()):
        subs_dict[dummy_species[i]] = subs_dict.pop(key)
        value[dummy_species[i]] = value.pop(key)

    updated_dict = copy.deepcopy(subs_dict)
    for key in subs_dict.keys():
        if subs_dict[key][key] == 0:
            updated_dict.pop(key)

    return updated_dict
