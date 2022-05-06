# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

from typing import FrozenSet, Iterable, Optional, Tuple

import networkx as nx
from openff.toolkit.utils.serialization import Serializable

from .chemicalsystem import ChemicalSystem
from .transformations import Transformation

from .executors.client import Client


class AlchemicalNetwork(Serializable):
    """A network of `ChemicalSystem`s as nodes, `Transformation`s as edges.

    Attributes
    ----------
    edges : FrozenSet[Transformation]
        The edges of the network, given as a `frozenset` of `Transformation`s.
    nodes : FrozenSet[ChemicalSystem]
        The nodes of the network, given as a `frozenset` of `ChemicalSystem`s.
    client : Optional[Client]
        A `Client`, which exposes results of evaluating `Transformation`s
        defined on the network via an `Executor`.

    """

    def __init__(
        self,
        edges: Iterable[Transformation] = None,
        nodes: Iterable[ChemicalSystem] = None,
        client: Optional[Client] = None
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
                       frozenset(e.initial for e in self._edges) |
                       frozenset(e.final for e in self._edges)
                      )

        self._graph = None

        self._client = client

    def __eq__(self, other):
        return self.nodes == other.nodes and self.edges == other.edges

    @staticmethod
    def _generate_graph(edges, nodes):
        g = nx.MultiDiGraph()

        for transformation in edges:
            g.add_edge(transformation.initial, transformation.final, object=transformation)

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

    @property
    def client(self) -> Client:
        return self._client

    @client.setter
    def client(self, client: Client):
        # TODO add validation checks
        self._client = client

    def estimate(self, transformations: Iterable[Transformation], estimator: str = None) -> Tuple[Tuple[float],Tuple[float],Tuple[float]]:
        """Get free energy estimates, uncertainties, and rates of convergence
        for the given transformations.

        Requires `client` to be defined on this `AlchemicalNetwork`.

        Paramters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.
        estimator : str
            Estimator to use in deriving transformation estimates using all
            other included transformations.

        Returns
        -------
        dG : Iterable[float]
            Free energy estimates for the given transformations, in order.
        ddG : Iterable[float]
            Uncertainties in dG for the given transformations, in order.
        rate_of_convergence : Iterable[float]
            Rate of convergence for dG for the given transformations, in order.
        """
        ...

    def estimate_dG(self, transformations: Iterable[Transformation], estimator: str = None) -> Tuple[float]:
        """Get free energy estimates for the given transformations.

        Requires `client` to be defined on this `AlchemicalNetwork`.

        Paramters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.
        estimator : str
            Estimator to use in deriving transformation estimates using all
            other included transformations.

        Returns
        -------
        dG : Iterable[float]
            Free energy estimates for the given transformations, in order.
        """
        ...

    def estimate_uncertainty(self, transformations: Iterable[Transformation], estimator: str = None) -> Tuple[float]:
        """Get free energy estimates for the given transformations.

        Requires `client` to be defined on this `AlchemicalNetwork`.

        Paramters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.
        estimator : str
            Estimator to use in deriving transformation estimates using all
            other included transformations.

        Returns
        -------
        ddG : Iterable[float]
            Uncertainties in dG for the given transformations, in order.
        """
        ...

    def estimate_rate_of_convergence(self, transformations: Iterable[Transformation], estimator: str = None) -> Tuple[float]:
        """Get free energy estimates for the given transformations.

        Requires `client` to be defined on this `AlchemicalNetwork`.

        Paramters
        ---------
        transformations : Iterable[Transformation]
            Transformations to retrieve estimates for.
        estimator : str
            Estimator to use in deriving transformation estimates using all
            other included transformations.

        Returns
        -------
        rate_of_convergence : Iterable[float]
            Rate of convergence for dG for the given transformations, in order.
        """
        ...
