# surface-pd Examples

This directory contains example data for use with the scripts in the main `scripts` directory.

## Directory Structure

```
examples/
├── README.md                         # This file
├── enumeration-examples/             # Input files and structures for enumeration
│   ├── input/                        # JSON configuration files
│   └── structure/                    # VASP structure files
│       ├── catalyst_surface/         # Catalyst examples (Pt)
│       ├── electrode/                # Electrode examples (Li, Li2O, LCO, LNO)
│       ├── metal/                    # Metal surface examples (GaAs)
│       └── other/                    # Other examples (NaTiS2)
└── plotting-examples/                # DFT data for phase diagram construction
    ├── *.dat                         # Energy data files
    └── LNO-001/, LNO-104/            # Specific system data
```

## Quick Start

See the docs for further details.

### 1. Surface Enumeration

The script `surface-enumeration.py` enumerates species replacements in surface slab structures.  It takes as input a JSON file with the species the replacement instructions (species to be replaced, fractions to be sampled, substitute species), the path to the slab model structure, and other parameters (e.g., whether or not to create symmetric slab models).  Example input files are in `enumeration-examples/input`.  

**Example:**

```bash
surface-enumeration.py examples/enumeration-examples/input/input-Li.json
```

### 2. Phase Diagram Plotting

The script `surface-pd-plot.py` generates surface phase diagrams for electrode surfaces with the approach detailed in [Li et al., *ACS Appl. Energy Mater.* **5**, 2022, 5730–5741](https://doi.org/10.1021/acsaem.2c00012).  

**Example:**

```bash
surface-pd-plot.py examples/plotting-examples/Li-surface_pd_data-symm-all.dat \
    --lithium-like-species Li \
    --oxygen-like-species O \
    --functional PBE
```

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
  version = {1.0.0},
  year = {2025}
}
```
