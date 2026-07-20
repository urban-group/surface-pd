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

The `surface-pd-plot` command generates surface phase diagrams for electrode
surfaces with the approach detailed in [Li et al., *ACS Appl. Energy Mater.*
**5**, 2022, 5730–5741](https://doi.org/10.1021/acsaem.2c00012).

**Example:**

```bash
surface-pd-plot examples/plotting-examples/SCAN-Li-surface.dat \
    --lithium-like-species Li \
    --oxygen-like-species O
```

Each plotting data file starts with the Li, O2, and bulk LiTMO2 reference
energies required for that dataset. Method details are retained as provenance;
unknown settings are marked explicitly, for example ``U_Ni=? eV``.

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
