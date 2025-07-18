# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from collections import abc
from typing import Optional

from .components import Component
from .tokenization import GufeTokenizable


class ChemicalSystem(GufeTokenizable, abc.Mapping):
    def __init__(
        self,
        components: dict[str, Component],
        name: str | None = "",
    ):
        r"""A combination of Components that form a system

        Containing a combination of :class:`.SmallMoleculeComponent`,
        :class:`.SolventComponent` and :class:`.ProteinComponent`, this object
        typically represents all the molecules in a simulation box.

        Used as a node for an :class:`.AlchemicalNetwork`.

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
            components={key: value for key, value in d["components"].items()},
            name=d["name"],
        )

    @property
    def components(self) -> dict[str, Component]:
        """
        The individual components of the chemical system.

        Components include atomic connectivity and coordinates. This is a
        dict with user-defined labels as keys and :class:`.Component`
        instances as values.
        """
        return dict(self._components)

    def component_diff(
        self, other
    ) -> tuple[tuple[Component, ...], tuple[Component, ...]]:
        """Compare the Components of this ChemicalSystem with the Components of another ChemicalSystem.

        Parameters
        ----------
        other : ChemicalSystem
            The ChemicalSystem to compare to.

        Returns
        -------
        tuple[tuple[Component, ...], tuple[Component, ...]]
            A tuple containing two tuples. The first tuple contains
            the components that are unique to this ChemicalSystem, and
            the second tuple contains the components that are unique
            to the other ChemicalSystem.

        Raises
        ------
        TypeError
            If `other` is not an instance of `ChemicalSystem`.

        """

        if not isinstance(other, ChemicalSystem):
            raise TypeError(
                f"`other` must be an instance of `{ChemicalSystem.__qualname__}`, not `{other.__class__.__qualname__}`"
            )

        self_comps = set(self._components.values())
        other_comps = set(other._components.values())

        self_uniques = tuple(self_comps.difference(other_comps))
        other_uniques = tuple(other_comps.difference(self_comps))

        return (self_uniques, other_uniques)

    @property
    def name(self):
        """
        Optional identifier for the chemical system.

        Used as part of the (hashable) graph node itself when the chemical state
        is added to an :class:`.AlchemicalNetwork`.
        """
        return self._name

    @property
    def total_charge(self):
        """Formal charge for the ChemicalSystem."""
        # This might evaluate the property twice?
        # return sum(component.total_charge
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
    def _defaults(cls):
        return super()._defaults()
