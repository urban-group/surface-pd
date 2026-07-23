# Surface-enumeration command-line walkthrough

Run every command in this guide from the `examples/` directory. The paths in
the committed JSON files are relative to that directory.

```bash
cd examples
```

The command-line enumerator requires the separately installed enumlib
executables `enum.x` and `makestr.x` (or a supported `makeStr` variant) on
`PATH`. Confirm the installed surface-pd command and its options with:

```bash
surface-enumeration --help
```

## Inspect the input

The LiCoO2 example uses
`enumeration-examples/input/input-LCO.json`. It identifies the VASP slab,
selects two Li-bearing layers and one O-bearing layer on each surface, requests
symmetric enumeration, and permits surface cells up to twice the parent area.
Its replacement fractions describe Li and O vacancy concentrations.

The referenced slab already contains selective-dynamics flags. Fixed sites
remain fixed in generated structures, while sites in the enumerated surface
regions remain relaxed.

## Enumerate structures

Run the complete LiCoO2 example without writing structure files:

```bash
surface-enumeration enumeration-examples/input/input-LCO.json
```

The command reports the target surface-cell size and the number of accepted
structures for every requested Li/O occupancy pair, followed by the total.
Candidate generation can take some time because enumlib evaluates the full
composition cross product.

For a smaller single-species example, use:

```bash
surface-enumeration enumeration-examples/input/input-Li.json
```

## Write VASP structures

Add `--generate-poscar` (or `-g`) to write every accepted structure:

```bash
surface-enumeration \
    enumeration-examples/input/input-LCO.json \
    --generate-poscar
```

The command organizes POSCAR files in composition-specific directories below
the current directory. Run it in a disposable working copy when you do not
want generated files beside the examples.

Before using a generated model in a calculation, inspect its lattice,
composition, surface ordering, and selective-dynamics flags with pymatgen or a
structure viewer.

## Related material

- [`enumeration-python-api.ipynb`](enumeration-python-api.ipynb) demonstrates
  the high-level Python API, asymmetric and symmetric enumeration, result
  metadata, and controlled POSCAR export.
- The [Sphinx command reference](../docs/source/tutorials-surface-enum.rst)
  documents the input fields and durable command semantics.
