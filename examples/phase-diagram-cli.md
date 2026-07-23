# Phase-diagram command-line walkthrough

Run every command in this guide from the `examples/` directory so that all
configuration-relative table paths resolve consistently.

```bash
cd examples
```

Inspect the installed command and current options with:

```bash
surface-pd-plot --help
```

## Render the maintained charge example

The SCAN+rVV10+U LiNiO2 (001) charge configuration is
`plotting-examples/lno-001-scan/charge.json`. It declares the thermodynamic
model, axes, explicit reference energies, phase table, and column mappings.

Render it over an explicit voltage-temperature domain:

```bash
surface-pd-plot \
    plotting-examples/lno-001-scan/charge.json \
    --x-range 0 5 \
    --y-range 1 1500 \
    --output lno-001-charge.pdf
```

The filename extension selects any Matplotlib-supported format, including
PDF, PNG, and SVG. The command does not invent scientific axis ranges or an
output destination.

Use `--show` for an interactive window. It may be combined with `--output`:

```bash
surface-pd-plot \
    plotting-examples/lno-001-scan/charge.json \
    --x-range 0 5 \
    --y-range 1 1500 \
    --output lno-001-charge.png \
    --show
```

## Adjust command-line presentation

`--mesh-points N` changes the shared linear sampling density. The default is
201 points per axis. Repeat `--condition NAME=VALUE` for configured state
variables that are not assigned to an axis.

The default coloring shows the atomic fraction of the first independent
component. Select another component with `--color-component`, or request
discrete stable-phase colors with `--phase-identity-colors`. Those two options
are mutually exclusive.

For publication styling, unequal coordinate arrays, or programmatic access to
the thermodynamic result, use the Python API rather than expanding the command
line with plotting details.

## Charge and discharge datasets

The migrated discharge examples contain a scientifically selected subset of
the corresponding charge candidates while preserving their original DFT
energies. Surface-pd does not infer an oxygen-loss history or treat zero energy
as an exclusion marker; dataset accessibility is explicit user policy.

## Related material

- [`phase-diagram-python-api.ipynb`](phase-diagram-python-api.ipynb) evaluates
  and renders the same scientific workflow step by step.
- The [Sphinx command reference](../docs/source/tutorials-surface-plot.rst)
  documents configuration semantics, workflow guarantees, and coloring
  behavior.
