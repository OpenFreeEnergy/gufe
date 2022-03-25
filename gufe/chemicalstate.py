# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Dict

import numpy as np
from frozendict import frozendict
from openff.toolkit.topology import Molecule, Topology
from openff.interchange.components.interchange import Interchange

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
            box_vectors: np.ndarray,
            identifier: str = None,
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
        self.components = frozendict(components)
        self.box_vectors = box_vectors
        self.identifier = identifier

    def __hash__(self):
        return hash((self.components, self.identifier))


class ProteinLigandSolventMicrostate(ChemicalState):
    """

    """
    ...


class LigandSolvent(ChemicalState):
    """

    """
    ...

    pass
