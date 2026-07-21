# surface-pd Examples

This directory contains input structures, configuration files, and phase data
for the installed surface-pd commands.

## Directory Structure

```
examples/
├── README.md                         # This file
├── enumeration-examples/             # Input files and structures for enumeration
│   ├── input/                        # JSON configuration files
│   └── structure/                    # VASP structure files
│       ├── catalyst_surface/         # Catalyst structures (Pt)
│       ├── electrode/                # Electrode structures (Li, Li2O, LCO, LNO)
│       ├── metal/                    # Semiconductor structures (GaAs)
│       └── other/                    # Other structures (NaTiS2)
└── plotting-examples/                # DFT data for phase diagram construction
    ├── *.dat                         # Energy data files
    └── LNO-001/, LNO-104/            # Specific system data
```

## Quick Start

See the docs for further details.

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
surface-pd-plot CONFIG.json --output diagram.pdf
```

The existing files under `plotting-examples/` use the former Li/O-specific
format. They are retained temporarily as numerical regression inputs and are
not accepted directly by the generalized command. Issue 42 tracks reproduction
of their validated results, migration to version-1 JSON plus explicitly mapped
tables, and removal of the legacy files and implementation before the initial
release.

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
