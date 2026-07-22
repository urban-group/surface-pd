# surface-pd Examples

This directory contains input structures, configuration files, and phase data
for the installed surface-pd commands.

## Directory Structure

```
examples/
├── README.md                         # This file
├── phase-diagram-python-api.ipynb    # Reviewed Python API tutorial
├── enumeration-examples/             # Input files and structures for enumeration
│   ├── input/                        # JSON configuration files
│   └── structure/                    # VASP structure files
│       ├── catalyst_surface/         # Catalyst structures (Pt)
│       ├── electrode/                # Electrode structures (Li, Li2O, LCO, LNO)
│       ├── metal/                    # Semiconductor structures (GaAs)
│       └── other/                    # Other structures (NaTiS2)
└── plotting-examples/                # Generalized phase-diagram inputs
    ├── lno-001-pbe/                  # PBE+U charge configuration
    ├── lno-001-scan/                 # SCAN charge and discharge examples
    ├── lno-104-pbe/                  # PBE+U charge configuration
    └── lno-104-scan/                 # SCAN charge and discharge examples
```

## Quick Start

The maintained
[phase-diagram Python API tutorial](phase-diagram-python-api.ipynb) evaluates
the SCAN+rVV10+U LiNiO2 (001) example step by step. It is intentionally more
explicit than the command-line workflow so that the model, datasets,
thermodynamic result, and rendering API can be inspected independently.

See the docs for reference details.

### 1. Surface Enumeration

The `surface-enumeration` command enumerates vacancies or substitutions in
surface slab structures. It reads a JSON configuration containing the slab
path, replacement fractions, selected layers, maximum cell size, and symmetry
choice. Example configurations are in `enumeration-examples/input`.

**Example:**

```bash
surface-enumeration examples/enumeration-examples/input/input-Li.json
```

### 2. Phase Diagram Plotting

The `surface-pd-plot` command evaluates one versioned generalized configuration:

**Example:**

```bash
surface-pd-plot \
    examples/plotting-examples/lno-001-scan/charge.json \
    --x-range 0 5 --y-range 1 1500 \
    --output lno-001-charge.pdf
```

The plotting examples use version-1 JSON and explicitly mapped phase tables.
Distinct PBE+U and SCAN+rVV10+U calculations are retained. The discharge
examples contain only accessible candidate phases and preserve their original
charge-state DFT energies; they do not use zero-energy exclusion sentinels.

## Getting Help

- Full documentation: https://surface-pd.readthedocs.io
- Report issues: https://github.com/urban-group/surface-pd/issues
- Contact: a.urban@columbia.edu

## Citation

If you use these examples or the surface-pd package in your research, please cite:

```bibtex
@software{surface_pd,
  author = {Li, Xinhao and Urban, Alexander},
  title = {surface-pd: Surface Phase Diagram Generator},
  url = {https://github.com/urban-group/surface-pd},
  version = {0.1.0},
  year = {2024}
}
```
