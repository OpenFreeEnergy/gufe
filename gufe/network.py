# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import FrozenSet, Iterable, Optional, Tuple

import networkx as nx
from .tokenization import GufeTokenizable

from .chemicalsystem import ChemicalSystem
from .transformations import Transformation


class AlchemicalNetwork(GufeTokenizable):
    """A network of `ChemicalSystem`s as nodes, `Transformation`s as edges.

    Attributes
    ----------
    edges : FrozenSet[Transformation]
        The edges of the network, given as a `frozenset` of `Transformation`s.
    nodes : FrozenSet[ChemicalSystem]
        The nodes of the network, given as a `frozenset` of `ChemicalSystem`s.
    name : Optional identifier for the network.

    """
    def __init__(
        self,
        edges: Iterable[Transformation] = None,
        nodes: Iterable[ChemicalSystem] = None,
        name: str = None,
    ):
        self._edges: FrozenSet[Transformation] = frozenset(edges) if edges else frozenset()
        self._nodes: FrozenSet[ChemicalSystem]

        self._name = name

        # possible to get more nodes via edges above,
        # so we merge these together
        if nodes is None:
            self._nodes = frozenset()
        else:
            self._nodes = frozenset(nodes)

        self._nodes = (
            self._nodes
            | frozenset(e.stateA for e in self._edges)
            | frozenset(e.stateB for e in self._edges)
        )

        self._graph = None

    def __eq__(self, other):
        return self.nodes == other.nodes and self.edges == other.edges

    @staticmethod
    def _generate_graph(edges, nodes):
        g = nx.MultiDiGraph()

        for transformation in edges:
            g.add_edge(
                transformation.stateA, transformation.stateB, object=transformation
            )

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

    @property
    def name(self):
        return self._name

    def _to_dict(self) -> dict:
        return {"nodes": sorted(self.nodes),
                "edges": sorted(self.edges),
                "name": self.name}

    @classmethod
    def _from_dict(cls, d: dict):
        return cls(nodes=frozenset(d['nodes']),
                   edges=frozenset(d['edges']),
                   name=d.get('name'))

    def _defaults(self):
        return super()._defaults()

    def to_graphml(self) -> str:
        """ """

    @classmethod
    def from_graphml(self, str):
        """ """
