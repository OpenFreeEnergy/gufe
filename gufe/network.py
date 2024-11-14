# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import Generator, Iterable, Optional, Self

import networkx as nx
from .tokenization import GufeTokenizable

from .chemicalsystem import ChemicalSystem
from .transformations import Transformation


class AlchemicalNetwork(GufeTokenizable):
    _edges: frozenset[Transformation]
    _nodes: frozenset[ChemicalSystem]
    _name: Optional[str]
    _graph: Optional[nx.MultiDiGraph]  # lazily created and cached

    """A network with all the information needed for a simulation campaign.

    Nodes are :class:`.ChemicalSystem` instances and edges are
    :class:`.Transformation` instances.

    Parameters
    ----------
    edges : Optional[Iterable[Transformation]]
      the links between chemical states
    nodes : Optional[Iterable[ChemicalSystem]]
      the individual chemical states.  :class:`.ChemicalSystem` objects from
      Transformation objects in edges will be automatically extracted
    """
    def __init__(
        self,
        edges: Optional[Iterable[Transformation]] = None,
        nodes: Optional[Iterable[ChemicalSystem]] = None,
        name: Optional[str] = None,
    ):
        self._edges = frozenset(edges) if edges else frozenset()
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
    def _generate_graph(edges, nodes) -> nx.MultiDiGraph:
        g = nx.MultiDiGraph()

        for transformation in edges:
            g.add_edge(
                transformation.stateA, transformation.stateB, object=transformation
            )

        g.add_nodes_from(nodes)

        return nx.freeze(g)

    @property
    def graph(self) -> nx.MultiDiGraph:
        """A networkx representation of the AlchemicalNetwork

        Nodes are represented as :class:`.ChemicalSystem` objects and directed
        edges are represented as :class:`.Transformation` objects
        """
        if self._graph is None:
            self._graph = self._generate_graph(self._edges, self._nodes)
        return self._graph.copy(as_view=True)

    @property
    def edges(self) -> frozenset[Transformation]:
        """
        Network edges as a `frozenset` of :class:`.Transformation` instances.
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
    def _from_dict(cls, d: dict) -> Self:
        return cls(nodes=frozenset(d['nodes']),
                   edges=frozenset(d['edges']),
                   name=d.get('name'))

    @classmethod
    def _defaults(cls):
        return super()._defaults()

    def to_graphml(self) -> str:
        """Currently not implemented"""
        raise NotImplementedError

    @classmethod
    def _from_nx_graph(cls, nx_graph) -> Self:
        chemical_systems = [n for n in nx_graph.nodes()]
        tranformations = [e[2]['object'] for e in nx_graph.edges(data=True)]  # list of transformations
        return cls(nodes=chemical_systems, edges=tranformations)

    def connected_subgraphs(self) -> Generator[Self, None, None]:
        """Return all connected subgraphs of the alchemical network in order of size
        """
        node_groups = nx.weakly_connected_components(self.graph)
        for node_group in node_groups:
            nx_subgraph = self.graph.subgraph(node_group)
            alc_subgraph = self._from_nx_graph(nx_subgraph)
            yield(alc_subgraph)

    @classmethod
    def from_graphml(cls, str):
        """Currently not implemented"""
        raise NotImplementedError
