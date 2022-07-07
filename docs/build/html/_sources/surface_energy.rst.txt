The brief overview of how the surface energy is calculated for the both
non-polar and polar surface facets is described here.
For a detailed description of how all the equations listed below are
derived, please refer to the paper.

`Understanding the Onset of Surface Degradation in LiNiO2 Cathodes
<https://pubs.acs.org/doi/10.1021/acsaem.2c00012>`__

==============
Surface energy
==============

..
    Substitutions

.. |LiTMO2| replace:: LiTMO\ :sub:`2`
.. |O2| replace:: O\ :sub:`2`

..
    Content

In thermodynamic equilibrium, the |LiTMO2| (TM represents transition metal
species) layered surface
forms the reconstruction
with the lowest Gibbs free energy for the given conditions.
The grand-canonical surface phase diagram is thus determined by the surface
free energies of all possible surface structures with different arrangements of
lithium and oxygen vacancies, subject to the chemical potentials of lithium and
oxygen.

To derive an expression of the surface free energy, we consider the formal
truncation of the extended |LiTMO2| crystal structure along a lattice plane.
In thermodynamic equilibrium with oxygen and lithium reservoirs (e.g., the
reactants during the synthesis of the material), the surfaces may absorb or
release Li and O atoms, so that the surface stoichiometry can differ from the
stoichiometry of the |LiTMO2| bulk.
A surface slab model of any surface, whether ideal stoichiometric or
reconstructed, can thus be thought of as the result of the formal formation
reaction:

.. math::

 \frac{n_\text{TM}^\text{slab}}{n_\text{TM}^\text{bulk}}\text{LiTMO}_2^{\text
    {bulk}}
    + \Bigg(n_\text{Li}^\text{slab} -
            \frac{n_{\text{TM}}^{\text{slab}}}{n_{\text{TM}}^{\text{bulk}}}n_
    {\text{Li}}^{\text{bulk}}
    \Bigg)\text{Li}
    + \\
    \frac{1}{2} \Bigg(n_\text{O}^\text{slab} -
            \frac{n_\text{TM}^\text{slab}}{n_\text{TM}^\text{bulk}}n_{\text{O
    }}^{\text{bulk}}
    \Bigg)\text{O}_2
    \rightarrow
    \text{Li}_{n_{\text{Li}}^{\text{slab}}}\text{TM}_{n_{\text{TM}}^{\text{slab
    }}}O_{n_{\text{O}}^{\text{slab}}}

The `surface free energy` is the reaction free energy normalized by the
surface area `A` and is given by:

.. math::
  \gamma = \frac{1}{2\text{A}}\Bigg[
  G_{\text{slab}} - \frac{n_{\text{TM}}^{\text{slab}}}{n_{\text{TM}}^{\text
    {bulk}}} G_{\text{bulk}}
	- \sum_{{i}}^{\text{Li},
  \text{O}}(n_{{i}}^{\text{slab}} -
  \frac{n_{\text{TM}}^{\text{slab}}}{n_{\text{TM}}^{\text{bulk}}}
  n_{{i}}^{\text{bulk}}) \mu_{\text{i}}
  \Bigg]

where :math:`n_{i}^{\text{slab}}` and :math:`n_{i}^{\text{bulk}}` are the
number of atoms of species :math:`i` (O and Li) in the slab and bulk
models, respectively. :math:`G_{\text{slab}} = G
(\text{Li}_{n_{\text{Li}}^{\text{slab}}}\text{TM}_{n_{\text{TM}}^{\text{slab
}}}O_{n_{\text{O}}^{\text{slab}}})`
and :math:`G_{\text{bulk}}=G(\text{LiTMO}_2^{\text{bulk}})`
are the Gibbs free energy of the slab and
bulk models, respectively, and the TM content is
assumed to be constant.
Neglecting the temperature dependence of the solids, :math:`G_\text{slab}`
and :math:`G_\text{bulk}` can be obtained from DFT calculations,
representing the energies of the slab
model and the bulk structure (one |LiTMO2| formula unit).

..
    For the fully lithiated |LiTMO2| bulk composition, the number of Li
    and TM atoms is identical, and we identify :math:`n_{\text{Li}}^{\text{bulk}}
    = n_{\text{TM}}^{\text{bulk}}`,
    and :math:`n_{\text{O}}^{\text{bulk}} = 2n_{\text{TM}}^{\text{bulk}}`.

In equilibrium, the chemical potential of Li can be calculated as:

.. math::
    \mu_{\text{Li}}=E(\text{Li}_{\text{bcc}})- F\, V

where :math:`E(\text{Li}_{\text{bcc}})` is the free energy of Li metal
in the body-centered cubic structure, which is also approximated with the
zero-Kelvin DFT energy.
`F` and `V` are Faraday's constant and cell potential, respectively.

The chemical potential of oxygen depends, in principle, on the temperature and
pressure, though the pressure dependence is negligible compared to the
temperature dependence. Thus, we ignore the impact of pressure.

.. math::
    \mu_{\text{O}}(T) = \frac{1}{2}\Bigg\{
    \mu_{\text{O}_2}^{0\text{K}} + \Bigg[ \Delta H^{\circ} + \Delta H(T) \Bigg]
    - T\Bigg[ S^{\circ} + \Delta S(T) \Bigg]
    \Bigg\}

where :math:`\frac{1}{2}` is a normalization factor accounting for the two
oxygen atoms in each |O2| molecule, and :math:`\mu_{\text{O}_2}^{0\text{K}}`
is the oxygen chemical potential at 0 Kelvin (equal to the enthalpy),
which was obtained from DFT calculations as described above.
The remaining terms are the enthalpy and entropy contributions to the
relative oxygen chemical potential.

..
    :math:`\mu_{\text{Li}}`, is
    equal to the sum of the Li ion and electron chemical potentials,
    :math:`{\mu_{\text{Li}} = \mu_{\text{Li}^+} + \mu_{\text{e}^-}}`.
    Using the Li metal electrode :math:`\text{Li} \leftrightarrow \text{Li}^+ +
    e^-` as reference, we introduce the reference chemical potentials
    :math:`\mu_{\text{Li}^+}^{\circ}` and
    :math:`\mu_{\text{e}^-}^{\circ}` with
    :math:`\mu_{\text{Li}^+}^{\circ} + \mu_{\text{e}^-}^{\circ} = G
    (\text{Li}_{\text{bcc}})`,
    The difference of the :math:`\text{Li}^+` and :math:`\text{e}^-` chemical
    potentials from the reference Li metal electrode determines the cell
    potential `V`: :math:`\Delta\mu_{\text{Li}^+} + \Delta\mu_{\text{e}^-}= - F\,
    V`,
    where `F` is Faraday's constant.

..
    The difference of the enthalpy at standard conditions from its
    0~\ Kelvin value, :math:`\Delta{}H^{\circ}`, and the standard entropy
    :math:`S^{\circ}` were taken from the NIST-JANAF Thermochemical
    Tables.
    :math:`\Delta H(T) = C_p(T - T_{\text{0}})` and :math:`\Delta S(T) = C_p\ln{(T -T_{\textup{0}})}`,
    where the heat capacity :math:`C_{p} = 3.5k_B` was taken to be the value for
    an ideal gas of diatomic molecules at :math:`T \geq298`.

Finally, taken together, the surface energy of any |LiTMO2| surface
(ideal stoichiometric or defected/reconstructed) can be approximated as:

.. math::
    \gamma(T, V)
  	\approx \frac{1}{2A} \Bigg\{
    E_{\text{slab}}
    + (n_{\text{TM}} - n_{\text{Li}}) \Bigg[
    E(\text{Li}_{\text{bcc}})
    - F\, V
    \Bigg] \\
    + \frac{1}{2} (2n_{\text{TM}} - n_{\text{O}}) \mu_{\text{O}_2}(T)
   	- n_{\text{TM}} E(\text{LiTMO}_2)
  \Bigg\}