# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import FrozenSet, Iterable, Union

import networkx as nx
from openff.toolkit.utils.serialization import Serializable

from .chemicalsystem import ChemicalSystem
from .transformation import Transformation


class AlchemicalNetwork(Serializable):
    """A network of `ChemicalSystem`s as nodes, `Transformation`s as edges.

    Attributes
    ----------
    chemicalstates :

    """

    def __init__(
        self,
        edges: Iterable[Transformation] = None,
        nodes: Iterable[ChemicalSystem] = None,
    ):

        self._edges = frozenset(edges) if edges else frozenset()

        # possible to get more nodes via edges above,
        # so we merge these together
        if nodes is None:
            self._nodes = frozenset()
        else:
            self._nodes = frozenset(nodes)

        self._nodes = (
                       self._nodes | 
                       frozenset(e.start for e in self._edges) |
                       frozenset(e.end for e in self._edges)
                      )

        self._graph = None

    def __eq__(self, other):
        return self.nodes == other.nodes and self.edges == other.edges

    @staticmethod
    def _generate_graph(edges, nodes):
        g = nx.MultiDiGraph()

        for transformation in edges:
            g.add_edge(transformation.start, transformation.end, object=transformation)

        g.add_nodes_from(nodes)

        return nx.freeze(g)

    @property
    def graph(self):
        if self._graph is None:
            self._graph = self._generate_graph(self._edges, self._nodes)
        return self._graph.copy(as_view=True)

    @property
    def edges(self) -> FrozenSet[Transformation]:
        return self._edges

    @property
    def nodes(self) -> FrozenSet[ChemicalSystem]:
        return self._nodes

    def to_dict(self) -> dict:
        ...

    @classmethod
    def from_dict(cls, d: dict):
        ...
