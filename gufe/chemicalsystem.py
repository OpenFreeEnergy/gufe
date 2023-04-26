# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from collections import abc
from typing import Dict, Optional

import numpy as np

from .tokenization import GufeTokenizable
from .components import Component


class ChemicalSystem(GufeTokenizable, abc.Mapping):
    """A node of an alchemical network.

    Attributes
    ----------
    components
        The molecular representation of the chemical state, including
        connectivity and coordinates. This is a frozendict with user-defined
        labels as keys, :class:`.Component`\ s as values.
    name
        Optional identifier for the chemical state; used as part of the
        (hashable) graph node itself when the chemical state is added to an
        :class:`.AlchemicalNetwork`.

    """

    def __init__(
        self,
        components: Dict[str, Component],
        name: Optional[str] = "",
    ):
        """Create a node for an alchemical network.

        Parameters
        ----------
        components
            The molecular representation of the chemical state, including
            connectivity and coordinates. Given as a dict with user-defined
            labels as keys, :class:`.Component`\ s as values.
        name
            Optional identifier for the chemical state; included with the other
            attributes as part of the (hashable) graph node itself when the
            chemical state is added to an :class:`.AlchemicalNetwork`.

        """
        super().__init__()

        self._components = components
        self._name = name

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name={self.name}, components={self.components})"
        )

    def _to_dict(self):
        return {
            "components": {
                key: value for key, value in sorted(self.components.items())
            },
            "name": self.name,
        }

    @classmethod
    def _from_dict(cls, d):
        return cls(
            components={
                key: value for key, value in d["components"].items()
            },
            name=d["name"],
        )

    @property
    def components(self):
        return dict(self._components)

    @property
    def name(self):
        return self._name

    @property
    def total_charge(self):
        """Formal charge for the ChemicalSystem."""
        # This might evaluate the property twice?
        #return sum(component.total_charge
        #           for component in self._components.values()
        #           if component.total_charge is not None)
        total_charge = 0
        for c in self._components.values():
            fc = c.total_charge
            if fc is not None:
                total_charge += fc
        return total_charge

    def __getitem__(self, item):
        return self.components[item]

    def __iter__(self):
        return iter(self.components)

    def __len__(self):
        return len(self.components)

    @classmethod
    def as_protein_smallmolecule_solvent(cls):
        """ """
        # alternate initializer for typical protein+ligand+solvent system
        ...

    @classmethod
    def as_smallmolecule_solvent(cls):
        """ """
        # alternate initializer for typical ligand+solvent system
        ...

    @classmethod
    def _defaults(cls):
        return super()._defaults()
