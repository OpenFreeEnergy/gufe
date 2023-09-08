# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Iterable, Optional

import networkx as nx
from .tokenization import GufeTokenizable

from .chemicalsystem import ChemicalSystem
from .transformations import Transformation


class AlchemicalNetwork(GufeTokenizable):
    """A network with all the information needed for a simulation campaign.

    Nodes are :class:`.ChemicalSystem` instances and edges are
    :class:`.Transformation` instances.
    """
    def __init__(
        self,
        edges: Optional[Iterable[Transformation]] = None,
        nodes: Optional[Iterable[ChemicalSystem]] = None,
        name: Optional[str] = None,
    ):
        self._edges: frozenset[Transformation] = frozenset(edges) if edges else frozenset()
        self._nodes: frozenset[ChemicalSystem]

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
    def edges(self) -> frozenset[Transformation]:
        """
        Network edges as a ``frozenset`` of :class:`.Transformation` instances.
        """
        return self._edges

    @property
    def nodes(self) -> frozenset[ChemicalSystem]:
        """
        Network nodes as a ``frozenset`` of :class:`.ChemicalSystem` instances.
        """
        return self._nodes

    @property
    def name(self) -> Optional[str]:
        """Optional identifier for the network."""
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

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    def to_graphml(self) -> str:
        """ """
        raise NotImplementedError

    @classmethod
    def from_graphml(self, str):
        """ """
        raise NotImplementedError
