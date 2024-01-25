# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from __future__ import annotations

import numpy as np
from openff.units import unit
from typing import Optional, Tuple, Literal

from .component import Component

_CATIONS = {'Cs', 'K', 'Li', 'Na', 'Rb'}
_ANIONS = {'Cl', 'Br', 'F', 'I'}
_ALLOWED_BOX_TYPES = {'cube', 'dodecahedron', 'octahedron'}


def _box_vectors_are_in_reduced_form(box_vectors: unit.Quantity) -> bool:
    """
    Return ``True`` if the box is in OpenMM reduced form; ``False`` otherwise.

    These conditions are shared by OpenMM and GROMACS and greatly simplify
    working with triclinic boxes. Any periodic system can be represented in this
    form by rotating the system and lattice reduction.
    See http://docs.openmm.org/latest/userguide/theory/05_other_features.html#periodic-boundary-conditions

    Acknowledgement
    ---------------
    Taken from openff.interchange's openff.interchange.components._packmol
    """
    if box_vectors.shape != (3, 3):
        return False

    a, b, c = box_vectors.m
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    return (
        [ay, az] == [0, 0]
        and bz == 0
        and ax > 0
        and by > 0
        and cz > 0
        and ax >= 2 * np.abs(bx)
        and ax >= 2 * np.abs(cx)
        and by >= 2 * np.abs(cy)
    )


# really wanted to make this a dataclass but then can't sort & strip ion input
class SolventComponent(Component):
    """
    :class:`Component` representing solvent molecules in a chemical system.

    This component represents the abstract idea of the solvent and ions present
    around the other components, rather than a list of specific water molecules
    and their coordinates.  This abstract representation is later made concrete
    by specific MD engine methods.
    """
    _smiles: str
    _positive_ion: str
    _negative_ion: str
    _neutralize: bool
    _ion_concentration: unit.Quantity

    def __init__(self, *,  # force kwarg usage
                 smiles: str = 'O',
                 num_solvent: Optional[int] = None,
                 positive_ion: str = 'Na+',
                 negative_ion: str = 'Cl-',
                 neutralize: bool = True,
                 ion_concentration: Optional[unit.Quantity] = 0.15 * unit.molar,
                 solvent_padding: Optional[unit.Quantity] = 1.2 * unit.nanometer,
                 solvent_density: Optional[unit.Quantity] = None,
                 box_shape: Optional[Literal['cube', 'dodecahedron', 'octahedron']] = 'cube',
                 box_vectors: Optional[unit.Quantity] = None,):
        """
        Parameters
        ----------
        smiles : str
          smiles of the solvent, default 'O' (water)
        num_solvent : int, optional
          number of solvent molecules present, if ``None``, will be determined
          by either `solvent_density` if defined, or the utilized solvation
          backend (Protocol specific). Cannot be defined as the same time
          as both `solvent_density` and either one of `solvent_padding` or
          `box_vectors`. Default `None`
        positive_ion, negative_ion : str
          the pair of ions which is used to neutralize (if neutralize=True) and
          bring the solvent to the required ionic concentration.  Must be a
          positive and negative monoatomic ions, defaults "Na+", "Cl-"
        neutralize : bool
          if the net charge on the chemical state is neutralized by the ions in
          this solvent component.  Default `True`
        ion_concentration : openff.units.unit.Quantity, optional
          ionic concentration required, default 0.15 * unit.molar
          this must be supplied with units, e.g. "1.5 * unit.molar"
        solvent_padding : openff.units.unit.Quantity, optional
          padding distance away from solute to use in constructing the
          solvent box. Must be supplied with units, e.g.
          "1.2 * unit.nanometer". Cannot be defined at the same time as
          `box_vectors`. Cannot be defined alongside both `num_solvent`
          and `solvent_density`. Default 1.2 * unit.nanometer
        solvent_density : openff.units.unit.Quantity, optional
          Target density of solvated systems. Must be defined with units
          compatible with g / mL, e.g. "850 * unit.kilogram / unit.meter**3".
          Cannot be defined alongside `num_solvent` if either `box_vector`
          or `solvent_padding` are defined. Default `None`
        box_shape : str, optional
          Defined the shape of the solvent box being built. Can be one of
          'cube', 'dodecahedron', and 'octahedron'. Cannot be defined alongside
          `box_vectors`. Default 'cube'
        box_vectors : openff.units.unit.Quantity, optional
          Vectors defining the solvent box. Box vectors must be provided in
          `OpenMM reduced form <http://docs.openmm.org/latest/userguide/theory/05_other_features.html#periodic-boundary-conditions>`_
          as a unit.Quantity array of shape (3, 3).
          Cannot be defined alongside `box_shape`.
          Cannot be defined alongside both `num_solvent` and
          `solvent_density`. Default `None`

        Raises
        ------
        ValueError
          * If an unknown positive or negative ion type is passed.
          * If `ion_concentration` is given in units incompatible with
            openff.units.unit.molar.
          * If `ion_concentration` is not positive.
          * If `num_solvent` is not `None` and not positive.
          * If `num_solvent` and `density` are defined alongside either
            `solvent_padding` or `box_vectors`.
          * If `solvent_padding` is defined alongside `box_vectors`.
          * If `solvent_padding` is given in units incompatible with
            openff.units.unit.nanometer.
          * If `solvent_padding` is not `None` and not positive.
          * If `solvent_density` is defined and not compatible with units of
            "g / mL".
          * If `solvent_desnity` is defined and not positive.
          * If an unknown box shape is passed.

        Examples
        --------
        To create a sodium chloride solution at 0.2 molar concentration with
        a cubic box with edges up to 1.2 nm from the solute::

          >>> s = SolventComponent(position_ion='Na', negative_ion='Cl',
          ...                      ion_concentration=0.2 * unit.molar,
          ...                      solvent_padding=1.2 * unit.nanometer,
          ...                      box_shape='cube')

        """
        self._smiles = smiles
        norm = positive_ion.strip('-+').capitalize()
        if norm not in _CATIONS:
            raise ValueError(f"Invalid positive ion, got {positive_ion}")
        positive_ion = norm + '+'
        self._positive_ion = positive_ion
        norm = negative_ion.strip('-+').capitalize()
        if norm not in _ANIONS:
            raise ValueError(f"Invalid negative ion, got {negative_ion}")
        negative_ion = norm + '-'
        self._negative_ion = negative_ion

        self._neutralize = neutralize

        if (not isinstance(ion_concentration, unit.Quantity)
              or not ion_concentration.is_compatible_with(unit.molar)):
            raise ValueError(f"ion_concentration must be given in units of"
                             f" concentration, got: {ion_concentration}")
        if ion_concentration.m < 0:
            raise ValueError("ion_concentration must be positive, "
                             f"got: {ion_concentration}")

        self._ion_concentration = ion_concentration

        # Validate num_solvent
        if num_solvent is not None and num_solvent < 1:
            errmsg = "num_solvent must be greater than zero, got {num_solvent}"
            raise ValueError(errmsg)

        if (num_solvent is not None and solvent_density is not None and
            (box_vectors is not None or solvent_padding is not None)):
            errmsg = ("Cannot define the number of solvent molecules "
                      f"({num_solvent}) alongside the solvent density "
                      f"({solvent_density}) and either of the box vectors "
                      f"({box_vectors}) or the solvent padding "
                      f"({solvent_padding})")
            raise ValueError(errmsg)

        self._num_solvent = num_solvent

        # Validate the solvent padding
        if solvent_padding is not None:

            if box_vectors is not None:
                errmsg = (f"solvent_padding ({solvent_padding}) cannot be "
                          f"defined alongside box_vectors ({box_vectors})")
                raise ValueError(errmsg)

            if (not isinstance(solvent_padding, unit.Quantity) or
                not solvent_padding.is_compatible_with(unit.nanometer)):
                errmsg = ("solvent_padding must be given in units of "
                          f"distance, got: {solvent_padding}")
                raise ValueError(errmsg)

            if solvent_padding.m < 0:
                errmsg = ("solvent_padding must be positive, "
                          f"got: {solvent_padding}")
                raise ValueError(errmsg)

        self._solvent_padding = solvent_padding

        # Validate solvent density
        if solvent_density is not None:
            if (not isinstance(solvent_density, unit.Quantity) or
                not solvent_density.is_compatible_with(unit.gram / unit.milliliter)):
                errmsg = ("solvent_density must be given in units compatible "
                          f"with g/mL, got: {solvent_density}")
                raise ValueError(errmsg)

            if solvent_density.m < 0:
                errmsg = ("solvent_density cannot be negative, "
                          f"got: {solvent_density}")
                raise ValueError(errmsg)

        self._solvent_density = solvent_density

        # Validate box shape
        if box_shape is not None:
            if box_shape.lower() not in _ALLOWED_BOX_TYPES:
                errmsg = (f"Unknown box_shape passed, got: {box_shape}, "
                          f"must be one of {', '.join(_ALLOWED_BOX_TYPES)}")
                raise ValueError(errmsg)

        self._box_shape = box_shape

        # Validate box vectors
        if box_vectors is not None:
            if not isinstance(box_vectors, unit.Quantity):
                errmsg = ("box_vectors must be defined as a unit.Quantity, "
                          f"got: {box_vectors}")
                raise ValueError(errmsg)

            if not _box_vectors_are_in_reduced_form(box_vectors):
                errmsg = ("box_vectors are not in reduced form, "
                          f"got: {box_vectors}")
                raise ValueError(errmsg)

        self._box_vectors = box_vectors

    @property
    def name(self) -> str:
        return (f"{self.smiles}, {self.positive_ion}, {self.negative_ion}, "
                f"{self.num_solvent}, {self.ion_concentration}, "
                f"{self.solvent_padding}, {self.solvent_density}, "
                f"{self.box_shape}, {self.box_vectors}")

    @property
    def smiles(self) -> str:
        """SMILES representation of the solvent molecules"""
        return self._smiles

    @property
    def positive_ion(self) -> Optional[str]:
        """The cation in the solvent state"""
        return self._positive_ion

    @property
    def negative_ion(self) -> Optional[str]:
        """The anion in the solvent state"""
        return self._negative_ion

    @property
    def neutralize(self) -> bool:
        """If the solvent neutralizes the system overall"""
        return self._neutralize

    @property
    def ion_concentration(self) -> unit.Quantity:
        """Concentration of ions in the solvent state"""
        return self._ion_concentration

    @property
    def total_charge(self):
        """Solvents don't have a formal charge defined so this returns None"""
        return None

    @property
    def num_solvent(self) -> Optional[int]:
        """Number of solvent molecules, if defined"""
        return self._num_solvent

    @property
    def solvent_padding(self) -> Optional[unit.Quantity]:
        """The solvent padding distance, if defined"""
        return self._solvent_padding

    @property
    def solvent_density(self) -> Optional[unit.Quantity]:
        """The solvent density, if defined"""
        return self._solvent_padding

    @property
    def box_shape(self) -> Optional[str]:
        """The solvated system box shape, if defined"""
        return self._box_shape

    @property
    def box_vectors(self) -> Optional[unit.Quantity]:
        """The solvated system box vectors, if defined"""
        return self._box_vectors

    @classmethod
    def _from_dict(cls, d):
        """Deserialize from dict representation"""
        ion_conc = d['ion_concentration']
        d['ion_concentration'] = unit.parse_expression(ion_conc)

        return cls(**d)

    def _to_dict(self):
        """For serialization"""
        ion_conc = str(self.ion_concentration)

        return {'smiles': self.smiles, 'positive_ion': self.positive_ion,
                'negative_ion': self.negative_ion,
                'ion_concentration': ion_conc,
                'neutralize': self._neutralize,
                'num_solvent': self._num_solvent,
                'solvent_padding': self._solvent_padding,
                'solvent_density': self._solvent_density,
                'box_shape': self._box_shape,
                'box_vectors': self._box_vectors}

    @classmethod
    def _defaults(cls):
        return super()._defaults()
