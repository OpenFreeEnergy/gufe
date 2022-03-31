"""

"""

from typing import FrozenSet, Iterable

import networkx as nx

from .chemicalstate import ChemicalState
from .transformation import Transformation


class AlchemicalNetwork:
    """A network of microstates as nodes, alchemical transformations as edges.

    Attributes
    ----------

    """

    def __init__(
            self,
            chemicalstates: Iterable[ChemicalState] = None,
            transformations: Iterable[Transformation] = None
            ):

        self._graph = None
