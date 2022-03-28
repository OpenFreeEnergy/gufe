# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Dict, Optional

import numpy as np
from frozendict import frozendict

from .custom_typing import RDKitMol


class ChemicalState:
    """A node of an alchemical network.

    Attributes
    ----------
    components
        The molecular representation of the chemical state, including
        connectivity and coordinates. This is a frozendict with user-defined
        labels as keys, RDKit molecules as values.
    box_vectors
        Numpy array indicating shape and size of unit cell for the system. May
        be a partial definition to allow for variability on certain dimensions.
    identifier
        Optional identifier for the chemical state; used as part of the
        (hashable) graph node itself when the chemical state is added to an
        `AlchemicalNetwork`

    """

    def __init__(
            self,
            components: Dict[str, RDKitMol],
            box_vectors: Optional[np.ndarray] = None,
            identifier: Optional[str] = None,
            ):
        """Create a node for an alchemical network.

        Attributes
        ----------
        components
            The molecular representation of the chemical state, including
            connectivity and coordinates. Given as a dict with user-defined
            labels as keys, RDKit molecules as values.
        box_vectors
            Numpy array indicating shape and size of unit cell for the system.
            May be a partial definition to allow for variability on certain
            dimensions.
        identifier
            Optional identifier for the chemical state; used as part of the
            (hashable) graph node itself when the chemical state is added to an
            `AlchemicalNetwork`

        """
        self._components = frozendict(components)
        self._identifier = identifier

        if box_vectors is None:
            self._box_vectors = np.array([np.nan]*9)
        else:
            self._box_vectors = box_vectors

    def __hash__(self):
        return hash((self._components, self._identifier, self._box_vectors))

    @property
    def components(self):
        return self._components

    @property
    def box_vectors(self):
        return np.array(self._box_vectors)

    @property
    def identifier(self):
        return self._identifier

    @classmethod
    def as_protein_ligand_solvent(cls):
        """

        """
        ...


    @classmethod
    def as_ligand_solvent(cls):
        """
    
        """
        ...
