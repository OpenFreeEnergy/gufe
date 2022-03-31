# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import FrozenSet, Iterable

import networkx as nx

from .chemicalstate import ChemicalState
from .transformation import Transformation


class AlchemicalNetwork:
    """A network of `ChemicalState`s as nodes, `Transformation`s as edges.

    Attributes
    ----------
    chemicalstates :

    """

    def __init__(
        self,
        transformations: Iterable[Transformation] = None,
        chemicalstates: Iterable[ChemicalState] = None,
    ):

        self._transformations = tuple(transformations) if transformations else tuple()
        self._chemicalstates = tuple(chemicalstates) if chemicalstates else tuple()

        self._graph = self._generate_graph(chemicalstates, transformations)

    def _generate_graph(self, nodes, edges):
        g = nx.MultiDiGraph()
        g.add_nodes_from(nodes)

        # THIS WILL NOT WORK
        # for each transformation, will need to add an edge to the graph using
        # start and end chemicalstates
        # not sure if really necessary to go further, so long as we have an easy way to
        # go from an edge in the graph itself to the corresponding Transformation object
        # with its Protocol
        g.add_edges_from(edges)

        return g

    @property
    def graph(self):
        return self._graph.copy(as_view=True)
