"""Chemical-potential models for surface thermodynamics."""

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from ._validation import (
    validate_finite_real,
    validate_positive_integer,
    validate_variable_name,
)
from .state import ThermodynamicState

_TEMPERATURE_VARIABLE = "temperature"
_VOLTAGE_VARIABLE = "voltage"
_STANDARD_TEMPERATURE_K = 298.15
_GAS_CONSTANT_KJ_PER_MOL_K = 0.008314463
_O2_ENTHALPY_INCREMENT_KJ_PER_MOL = 8.683
_O2_STANDARD_ENTROPY_KJ_PER_MOL_K = 205.147e-3
_O2_HEAT_CAPACITY_KJ_PER_MOL_K = 3.5 * _GAS_CONSTANT_KJ_PER_MOL_K
_KJ_PER_MOL_PER_EV = 96.487


@runtime_checkable
class ChemicalPotentialModel(Protocol):
    """Structural interface implemented by chemical-potential models.

    This protocol is satisfied by any object that declares its required state
    variables and evaluates a :class:`ThermodynamicState` to an array of
    chemical potentials in eV per associated composition component. Explicit
    inheritance is not required.
    """

    @property
    def required_state_variables(self) -> frozenset[str]:
        """Return state-variable names required for evaluation."""
        ...  # pragma: no cover - typing protocol declaration

    def evaluate(self, state: ThermodynamicState) -> np.ndarray:
        """Evaluate chemical potential over ``state.shape``.

        Parameters
        ----------
        state : ThermodynamicState
            Named thermodynamic variables for evaluation.

        Returns
        -------
        numpy.ndarray
            Chemical potential in eV per component with ``state.shape``.
        """
        ...  # pragma: no cover - typing protocol declaration


def _require_state(state: object) -> ThermodynamicState:
    """Return a validated state object or raise ``TypeError``."""
    if not isinstance(state, ThermodynamicState):
        raise TypeError("state must be a ThermodynamicState")
    return state


def _require_variable(
    state: ThermodynamicState,
    variable: str,
) -> np.ndarray:
    """Return one required state variable with a descriptive missing error."""
    if variable not in state:
        raise ValueError(
            f"Thermodynamic state is missing required variable {variable!r}"
        )
    return state[variable]


@dataclass(frozen=True)
class ConstantChemicalPotential:
    r"""Represent a state-independent chemical potential.

    Parameters
    ----------
    value_ev_per_component : float
        Finite chemical potential in eV per associated composition component.

    Notes
    -----
    The thermodynamic equation is

    .. math::

        \mu(\mathbf{x}) = \mu_0,

    where :math:`\mu_0` is ``value_ev_per_component``.
    """

    value_ev_per_component: float

    def __post_init__(self) -> None:
        """Validate and normalize the configured energy."""
        value = validate_finite_real(
            self.value_ev_per_component,
            "value_ev_per_component",
        )
        object.__setattr__(self, "value_ev_per_component", value)

    @property
    def required_state_variables(self) -> frozenset[str]:
        """Return an empty set because the model is state-independent."""
        return frozenset()

    def evaluate(self, state: ThermodynamicState) -> np.ndarray:
        """Return the constant potential broadcast over ``state.shape``.

        Parameters
        ----------
        state : ThermodynamicState
            State whose common shape determines the result shape.

        Returns
        -------
        numpy.ndarray
            Constant chemical potential in eV per component.

        Raises
        ------
        TypeError
            If ``state`` is not a :class:`ThermodynamicState`.
        """
        valid_state = _require_state(state)
        return np.full(
            valid_state.shape,
            self.value_ev_per_component,
            dtype=float,
        )


@dataclass(frozen=True)
class DirectChemicalPotential:
    r"""Vary a chemical potential directly relative to a reference.

    Parameters
    ----------
    state_variable : str
        Name of the state variable containing the offset in eV per component.
    reference_energy_ev_per_component : float
        Finite reference chemical potential in eV per component. Use zero for
        an absolute chemical-potential coordinate.

    Notes
    -----
    The thermodynamic equation is

    .. math::

        \mu = \mu_\mathrm{ref} + \Delta\mu.
    """

    state_variable: str
    reference_energy_ev_per_component: float

    def __post_init__(self) -> None:
        """Validate the state name and reference energy."""
        object.__setattr__(
            self,
            "state_variable",
            validate_variable_name(self.state_variable),
        )
        reference = validate_finite_real(
            self.reference_energy_ev_per_component,
            "reference_energy_ev_per_component",
        )
        object.__setattr__(
            self,
            "reference_energy_ev_per_component",
            reference,
        )

    @property
    def required_state_variables(self) -> frozenset[str]:
        """Return the configured chemical-potential offset variable."""
        return frozenset({self.state_variable})

    def evaluate(self, state: ThermodynamicState) -> np.ndarray:
        """Return reference plus direct offset over ``state.shape``.

        Parameters
        ----------
        state : ThermodynamicState
            State containing the configured offset variable.

        Returns
        -------
        numpy.ndarray
            Chemical potential in eV per component.

        Raises
        ------
        TypeError
            If ``state`` is not a :class:`ThermodynamicState`.
        ValueError
            If the configured offset variable is missing.
        """
        valid_state = _require_state(state)
        offset = _require_variable(valid_state, self.state_variable)
        return np.asarray(
            self.reference_energy_ev_per_component + offset,
            dtype=float,
        )


@dataclass(frozen=True)
class IntercalationChemicalPotential:
    r"""Represent an intercalant chemical potential controlled by voltage.

    Parameters
    ----------
    reference_energy_ev_per_component : float
        Finite reference energy in eV per intercalant component.
    electrons_per_component : int
        Positive number of electrons transferred per intercalant component.
    voltage_reference : str
        Nonempty, single-line provenance label for the voltage reference, such
        as ``"Li/Li+"``. The label is retained but not interpreted.

    Notes
    -----
    With voltage :math:`V` in volts and electron count :math:`z`, the model is

    .. math::

        \mu(V) = \mu_\mathrm{ref} - zV.
    """

    reference_energy_ev_per_component: float
    electrons_per_component: int
    voltage_reference: str

    def __post_init__(self) -> None:
        """Validate the reference energy, electron count, and provenance."""
        reference = validate_finite_real(
            self.reference_energy_ev_per_component,
            "reference_energy_ev_per_component",
        )
        electrons = validate_positive_integer(
            self.electrons_per_component,
            "electrons_per_component",
        )
        if (
            not isinstance(self.voltage_reference, str)
            or not self.voltage_reference.strip()
            or "\n" in self.voltage_reference
            or "\r" in self.voltage_reference
        ):
            raise ValueError(
                "voltage_reference must be a nonempty single-line string"
            )
        object.__setattr__(
            self,
            "reference_energy_ev_per_component",
            reference,
        )
        object.__setattr__(self, "electrons_per_component", electrons)
        object.__setattr__(
            self,
            "voltage_reference",
            self.voltage_reference.strip(),
        )

    @property
    def required_state_variables(self) -> frozenset[str]:
        """Return the standard voltage state variable."""
        return frozenset({_VOLTAGE_VARIABLE})

    def evaluate(self, state: ThermodynamicState) -> np.ndarray:
        """Evaluate the intercalant chemical potential over voltage.

        Parameters
        ----------
        state : ThermodynamicState
            State containing voltage in V.

        Returns
        -------
        numpy.ndarray
            Intercalant chemical potential in eV per component.

        Raises
        ------
        TypeError
            If ``state`` is not a :class:`ThermodynamicState`.
        ValueError
            If the voltage variable is missing.
        """
        valid_state = _require_state(state)
        voltage = _require_variable(valid_state, _VOLTAGE_VARIABLE)
        return np.asarray(
            self.reference_energy_ev_per_component
            - self.electrons_per_component * voltage,
            dtype=float,
        )


@dataclass(frozen=True)
class FixedPressureOxygenChemicalPotential:
    r"""Evaluate the preserved fixed-pressure O2 temperature approximation.

    Parameters
    ----------
    raw_o2_energy_ev_per_molecule : float
        Finite raw DFT reference energy in eV per O2 molecule.
    correction_ev_per_molecule : float
        Finite correction added to the raw reference in eV per O2 molecule.

    Raises
    ------
    ValueError
        If absolute temperature is not positive.

    Notes
    -----
    This model returns eV per O atom and preserves the legacy analytical
    approximation

    .. math::

        \mu_\mathrm{O}(T) = \frac{1}{2}\left[
        E_\mathrm{O_2}^\mathrm{DFT} + E_\mathrm{corr}
        + \frac{\Delta\mu_\mathrm{O_2}^\mathrm{thermal}(T)}
        {96.487}\right],

    where the thermal term is in kJ/mol O2,

    .. math::

        \Delta\mu_\mathrm{O_2}^\mathrm{thermal}(T) =
        \Delta H^\circ + C_p(T-T_0)
        - T\left[S^\circ + C_p\ln(T/T_0)\right].

    The constant-heat-capacity expression is anchored by NIST-JANAF values at
    298.15 K. It is not table interpolation, and pressure dependence is
    neglected.
    """

    raw_o2_energy_ev_per_molecule: float
    correction_ev_per_molecule: float

    def __post_init__(self) -> None:
        """Validate and normalize the O2 reference inputs."""
        raw = validate_finite_real(
            self.raw_o2_energy_ev_per_molecule,
            "raw_o2_energy_ev_per_molecule",
        )
        correction = validate_finite_real(
            self.correction_ev_per_molecule,
            "correction_ev_per_molecule",
        )
        object.__setattr__(self, "raw_o2_energy_ev_per_molecule", raw)
        object.__setattr__(self, "correction_ev_per_molecule", correction)

    @property
    def required_state_variables(self) -> frozenset[str]:
        """Return the standard absolute-temperature state variable."""
        return frozenset({_TEMPERATURE_VARIABLE})

    def evaluate(self, state: ThermodynamicState) -> np.ndarray:
        """Return oxygen chemical potential in eV per O atom.

        Parameters
        ----------
        state : ThermodynamicState
            State containing positive absolute temperature in K.

        Returns
        -------
        numpy.ndarray
            Oxygen chemical potential in eV per O atom.

        Raises
        ------
        TypeError
            If ``state`` is not a :class:`ThermodynamicState`.
        ValueError
            If temperature is missing or not positive.
        """
        valid_state = _require_state(state)
        temperature = _require_variable(valid_state, _TEMPERATURE_VARIABLE)
        if np.any(temperature <= 0):
            raise ValueError("Temperature must be positive in Kelvin.")

        delta_h = _O2_HEAT_CAPACITY_KJ_PER_MOL_K * (
            temperature - _STANDARD_TEMPERATURE_K
        )
        delta_s = _O2_HEAT_CAPACITY_KJ_PER_MOL_K * np.log(
            temperature / _STANDARD_TEMPERATURE_K
        )
        thermal_kj_per_mol = (
            _O2_ENTHALPY_INCREMENT_KJ_PER_MOL + delta_h
        ) - temperature * (_O2_STANDARD_ENTROPY_KJ_PER_MOL_K + delta_s)
        return np.asarray(
            0.5
            * (
                self.raw_o2_energy_ev_per_molecule
                + self.correction_ev_per_molecule
                + thermal_kj_per_mol / _KJ_PER_MOL_PER_EV
            ),
            dtype=float,
        )
