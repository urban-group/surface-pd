# surface-pd Examples

This directory contains input structures, configuration files, and phase data
for the installed surface-pd commands.

## Directory Structure

```
examples/
├── README.md                         # This file
├── enumeration-python-api.ipynb      # Reviewed enumeration API tutorial
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

Run all commands and notebooks in this directory from the `examples/`
directory. This keeps paths identical across the maintained examples:

```bash
cd examples
```

The maintained
[surface-enumeration Python API tutorial](enumeration-python-api.ipynb)
loads a committed Li(100) slab, identifies its layers and selected surface
sites, and prepares a vacancy enumeration. Its deterministic setup runs
without external tools; the final enumlib-dependent operation is clearly
marked and opt-in.

The maintained
[phase-diagram Python API tutorial](phase-diagram-python-api.ipynb) evaluates
the SCAN+rVV10+U LiNiO2 (001) example step by step. It is intentionally more
explicit than the command-line workflow so that the model, datasets,
thermodynamic result, and rendering API can be inspected independently.

The authoritative runnable command-line guides are:

- [Surface-enumeration command-line walkthrough](surface-enumeration-cli.md)
- [Phase-diagram command-line walkthrough](phase-diagram-cli.md)

Both guides use the committed inputs in this directory and keep every command
relative to `examples/`. Sphinx retains the corresponding input-format,
command-option, and scientific reference material without duplicating these
walkthroughs.

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
