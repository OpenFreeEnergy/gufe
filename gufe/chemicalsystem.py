# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from collections import abc
from typing import Dict, Optional

import numpy as np

from .tokenization import GufeTokenizable
from .components import Component


class ChemicalSystem(abc.Mapping, GufeTokenizable):
    """A node of an alchemical network.

    Attributes
    ----------
    components
        The molecular representation of the chemical state, including
        connectivity and coordinates. This is a frozendict with user-defined
        labels as keys, `Component`s as values.
    box_vectors
        Numpy array indicating shape and size of unit cell for the system. May
        be a partial definition to allow for variability on certain dimensions.
    name
        Optional identifier for the chemical state; used as part of the
        (hashable) graph node itself when the chemical state is added to an
        `AlchemicalNetwork`.

    """

    def __init__(
        self,
        components: Dict[str, Component],
        box_vectors: Optional[np.ndarray] = None,
        name: Optional[str] = "",
    ):
        """Create a node for an alchemical network.

        Attributes
        ----------
        components
            The molecular representation of the chemical state, including
            connectivity and coordinates. Given as a dict with user-defined
            labels as keys, `Component`s as values.
        box_vectors
            Optional `numpy` array indicating shape and size of unit cell for
            the system. May be a partial definition to allow for variability on
            certain dimensions.
        name
            Optional identifier for the chemical state; included with the other
            attributes as part of the (hashable) graph node itself when the
            chemical state is added to an `AlchemicalNetwork`.

        """
        super().__init__()

        self._components = components
        self._name = name

        if box_vectors is None:
            self._box_vectors = np.array([np.nan] * 9)
        else:
            self._box_vectors = box_vectors

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(name={self.name}, components={self.components})"
        )

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self._name != other.name:
            return False
        if not np.array_equal(self._box_vectors, other.box_vectors,
                              equal_nan=True):  # nan usually compares to false
            return False
        if self._components != other.components:
            return False

        return True

    def __hash__(self):
        return hash(
            (
                tuple(sorted(self._components.items())),
                self._box_vectors.tobytes(),
                self._name,
            )
        )

    def _to_dict(self):
        return {
            "components": {
                key: value for key, value in sorted(self.components.items())
            },
            "box_vectors": self.box_vectors.tolist(),
            "name": self.name,
        }

    @classmethod
    def _from_dict(cls, d):
        return cls(
            components={
                key: value for key, value in d["components"].items()
            },
            box_vectors=np.array(d["box_vectors"]),
            name=d["name"],
        )

    @property
    def components(self):
        return dict(self._components)

    @property
    def box_vectors(self):
        return np.array(self._box_vectors)

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

    def _defaults(self):
        return super()._defaults()
