# Table of Contents
- [Table of Contents](#table-of-contents)
- [Introduction](#introduction)
- [surface-pd Framework](#surface-pd-framework)
- [How to cite](#how-to-cite)
- [Documentation](#documentation)
- [Installation](#installation)
- [Using surface-pd](#using-surface-pd)


<a name="introduction"></a>
# Introduction
The surface degradation of layered transition metal (TM) oxides and other
related cathode compositions has been characterized extensively in
experiments. However, how the reactions take place on the cathode surface
and how surface reconstructions form on the atomic scale remains unclear.
The surface phase diagram as a function of state of charge and temperature
can offer insights into atomic-scale processes that are challenging to probe
experimentally.

<p align="center">
    <img src="docs/source/images/2D-3D-surface-pd-demo.png" width="500">
</p>

To construct the surface phase diagram, slab models
with different surface compositions need to be created. The surface phase
diagram shown on the left above only has four phases that are predicted to be
stable in the range of 0 to 5 V and 0 to 1500 K. However, behind the most
stable phases shown on the surface phase diagram, hundreds/thousands
meta-stable phases should also be calculated. The figure shown above on the
right illustrates how the three-dimensional surface phase diagram looks like
when "all" enumerated slab models are considered.

<a name="surface-pd-framework"></a>
# surface-pd Framework
This package uses the method developed in the [enumlib code](https://github.com/msg-byu/enumlib).
In contrast to the conventional enumlib code, we apply systematic enumeration
to only part of the slab model, creating vacancies on the surface. By taking
advantage of the inherent inversion symmetry center of the slab models, the
top surface is symmetrized to the bottom surface, i.e., a symmetrically
equivalent modification is introduced at the bottom of the slab models. The
flowchart below shows the framework of the surface enumeration.

<p align="center">
    <img src="docs/source/images/flowchart.png" width="500">
</p>


<a name="how-to-cite"></a>
# How to cite
If you have used ***surface-pd***, please cite the GitHub repository:

```bibtex
@software{surface_pd,
  author = {Li, Xinhao and Urban, Alexander},
  title = {surface-pd: Surface Phase Diagram Generator},
  url = {https://github.com/urban-group/surface-pd},
  version = {1.0.0},
  year = {2024}
}
```

Or in text format:
> Li, X., & Urban, A. (2024). surface-pd: Surface Phase Diagram Generator (Version 1.0.0) [Computer software]. https://github.com/urban-group/surface-pd

Thank you for your interest.

<a name="documentation"></a>
# Documentation
If you want to learn more about ***surface-pd***, please find the [user manual](https://surface-pd.readthedocs.io/en/latest/index.html).


Or please contact: Xinhao Li (<xinhao.li@columbia.edu>) and Alex Urban (<a.urban@columbia.edu>)

<a name="installation"></a>
# Installation
surface-pd is a Python 3 package and the stable version can be installed via pip.
```
$ pip install surface-pd
```

<a name="using-surface-pd"></a>
# Using surface-pd
Please refer to the official [surface-pd documentation](https://surface-pd.readthedocs.io/en/latest/index.html)
for tutorials and examples.
