# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from __future__ import annotations

import json
from collections.abc import Iterable
from itertools import chain
from typing import FrozenSet, Iterable, Optional, Union

import networkx as nx
import gufe
from gufe import SmallMoleculeComponent

from .mapping import LigandAtomMapping
from .tokenization import JSON_HANDLER, GufeTokenizable


class LigandNetwork(GufeTokenizable):
    """A directed graph connecting ligands according to their atom mapping.

    A network can be defined by specifying only edges, in which case the nodes are implicitly added.

    Parameters
    ----------
    edges : Iterable[LigandAtomMapping]
        Edges for this network, each specified as a LigandAtomMapping between two nodes.
    nodes : Iterable[SmallMoleculeComponent]
        Nodes for this network. Any nodes already included as a part of the 'edges' will be ignored.
        Nodes not already included in 'edges' will be added as isolated, unconnected nodes.
    """

    def __init__(
        self,
        edges: Iterable[LigandAtomMapping],
        nodes: Iterable[SmallMoleculeComponent] | None = None,
    ):
        if nodes is None:
            nodes = []

        self._edges = frozenset(edges)
        edge_nodes = set(chain.from_iterable((e.componentA, e.componentB) for e in edges))
        self._nodes = frozenset(edge_nodes) | frozenset(nodes)
        self._graph = None

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self) -> dict:
        return {"graphml": self.to_graphml()}

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls.from_graphml(dct["graphml"])

    @property
    def graph(self) -> nx.MultiDiGraph:
        """NetworkX graph for this network.

        This graph will have :class:`.SmallMoleculeComponent` objects as nodes and
        :class:`.LigandAtomMapping` objects as directed edges
        """
        if self._graph is None:
            graph = nx.MultiDiGraph()
            # set iterator order depends on PYTHONHASHSEED, sorting ensures
            # reproducibility
            for node in sorted(self._nodes):
                graph.add_node(node)
            for edge in sorted(self._edges):
                graph.add_edge(edge.componentA, edge.componentB, object=edge, **edge.annotations)

            self._graph = nx.freeze(graph)

        return self._graph

    @property
    def edges(self) -> frozenset[LigandAtomMapping]:
        """A read-only view of the edges of this network."""
        return self._edges

    @property
    def nodes(self) -> frozenset[SmallMoleculeComponent]:
        """A read-only view of the nodes of this network."""
        return self._nodes

    def _serializable_graph(self) -> nx.Graph:
        """Create a :mod:`networkx` graph with serializable attribute representations.

        This enables us to easily use different serialization approaches.
        """
        # sorting ensures that we always preserve order in files, so two
        # identical networks will show no changes if you diff their
        # serialized versions
        sorted_nodes = sorted(self.nodes, key=lambda m: (m.smiles, m.name))
        mol_to_label = {mol: f"mol{num}" for num, mol in enumerate(sorted_nodes)}

        edge_data = sorted(
            [
                (
                    mol_to_label[edge.componentA],
                    mol_to_label[edge.componentB],
                    json.dumps(list(edge.componentA_to_componentB.items())),
                    json.dumps(edge.annotations, cls=JSON_HANDLER.encoder),
                )
                for edge in self.edges
            ]
        )

        # from here, we just build the graph
        serializable_graph = nx.MultiDiGraph()
        for mol, label in mol_to_label.items():
            serializable_graph.add_node(label, moldict=json.dumps(mol.to_dict(), sort_keys=True))

        for molA, molB, mapping, annotation in edge_data:
            serializable_graph.add_edge(molA, molB, mapping=mapping, annotations=annotation)

        return serializable_graph

    @classmethod
    def _from_serializable_graph(cls, graph: nx.Graph):
        """Create network from :mod:`networkx` graph with serializable attributes.

        This is the inverse of ``_serializable_graph``.
        """
        label_to_mol = {
            node: SmallMoleculeComponent.from_dict(json.loads(d)) for node, d in graph.nodes(data="moldict")
        }

        edges = [
            LigandAtomMapping(
                componentA=label_to_mol[node1],
                componentB=label_to_mol[node2],
                componentA_to_componentB=dict(json.loads(edge_data["mapping"])),
                annotations=json.loads(
                    edge_data.get("annotations", "null"), cls=JSON_HANDLER.decoder
                ),  # work around old graphml files with missing edge annotations
            )
            for node1, node2, edge_data in graph.edges(data=True)
        ]

        return cls(edges=edges, nodes=label_to_mol.values())

    def to_graphml(self) -> str:
        """Return the GraphML string representing this network.

        This is the primary serialization mechanism for this class.

        Returns
        -------
        str
            String representing this network in GraphML format.
        """
        return "\n".join(nx.generate_graphml(self._serializable_graph()))

    @classmethod
    def from_graphml(cls, graphml_str: str) -> LigandNetwork:
        """Create from a GraphML string.

        Parameters
        ----------
        graphml_str : str
            GraphML string representation of a :class:`.Network`.

        Returns
        -------
        LigandNetwork
            New network from the GraphML.
        """
        return cls._from_serializable_graph(nx.parse_graphml(graphml_str))

    def enlarge_graph(self,
                      *,
                      edges: Optional[Iterable[LigandAtomMapping]] = None,
                      nodes: Optional[Iterable[SmallMoleculeComponent]] = None) -> LigandNetwork:
        """Create a new network with the given edges and nodes added.

        Parameters
        ----------
        edges : Iterable[:class:`.LigandAtomMapping`]
            Edges to append to this network.
        nodes : Iterable[:class:`.SmallMoleculeComponent`]
            Nodes to append to this network.

        Returns
        -------
        LigandNetwork
            A new network adding the given edges and nodes to this network.
        """
        if edges is None:
            edges = set()

        if nodes is None:
            nodes = set()

        return LigandNetwork(self.edges | set(edges), self.nodes | set(nodes))

    def trim_graph(self,
                     *,
                     edges: Optional[Iterable[LigandAtomMapping]] = None,
                     nodes: Optional[Iterable[SmallMoleculeComponent]] = None) -> LigandNetwork:
        """Create a new network with the given edges and nodes removed.

        Note that for removed ``nodes``, any edges that include them will also
        be removed.

        Parameters
        ----------
        edges : Iterable[:class:`.LigandAtomMapping`]
            Edges to drop from this network.
        nodes : Iterable[:class:`.SmallMoleculeComponent`]
            Nodes to drop from this network; all edges including these nodes
            will also be dropped.

        Returns
        -------
        LigandNetwork
            A new network with the given edges and nodes removed.
        """
        if edges is None:
            edges = list()

        if nodes is None:
            nodes = list()

        graph = self.graph.copy()
        graph.remove_nodes_from(nodes)
        graph.remove_edges_from([(edge.componentA, edge.componentB) for edge in edges])

        return LigandNetwork(edges=[obj for u, v, obj in graph.edges.data('object')],
                             nodes=graph.nodes)

    def _to_rfe_alchemical_network(
        self,
        components: dict[str, gufe.Component],
        leg_labels: dict[str, list[str]],
        protocol: gufe.Protocol,
        *,
        alchemical_label: str = "ligand",
        autoname=True,
        autoname_prefix="",
    ) -> gufe.AlchemicalNetwork:
        """Create an :class:`.AlchemicalNetwork` from this :class:`.LigandNetwork`.

        Parameters
        ----------
        components: dict[str, :class:`.Component`]
            Non-alchemical components (components that will be on both sides
            of a transformation).
        leg_labels: dict[str, list[str]]
            Mapping of the names for legs (the keys of this dict) to a list of
            the component names. The component names must be the same as used
            in the ``components`` dict.
        protocol: :class:`.Protocol`
            The protocol to apply.
        alchemical_label: str
            The label for the component undergoing an alchemical transformation
            (default ``'ligand'``).
        """
        transformations = []
        for edge in self.edges:
            for leg_name, labels in leg_labels.items():

                # define a helper func to avoid repeated code
                def sys_from_dict(component):
                    """
                    Input component alchemically changing. Other info taken
                    from the outer scope.
                    """
                    syscomps = {alchemical_label: component}
                    other_labels = set(labels) - {alchemical_label}
                    syscomps.update({label: components[label] for label in other_labels})

                    if autoname:
                        name = f"{component.name}_{leg_name}"
                    else:
                        name = ""

                    return gufe.ChemicalSystem(syscomps, name=name)

                sysA = sys_from_dict(edge.componentA)
                sysB = sys_from_dict(edge.componentB)
                if autoname:
                    prefix = f"{autoname_prefix}_" if autoname_prefix else ""
                    name = f"{prefix}{sysA.name}_{sysB.name}"
                else:
                    name = ""

                transformation = gufe.Transformation(sysA, sysB, protocol, mapping=edge, name=name)

                transformations.append(transformation)

        return gufe.AlchemicalNetwork(transformations)

    def to_rbfe_alchemical_network(
        self,
        solvent: gufe.SolventComponent,
        protein: gufe.ProteinComponent,
        protocol: gufe.Protocol,
        *,
        autoname: bool = True,
        autoname_prefix: str = "easy_rbfe",
        **other_components,
    ) -> gufe.AlchemicalNetwork:
        """Create an :class:`.AlchemicalNetwork` from this :class:`.LigandNetwork`.

        Parameters
        ----------
        protocol: Protocol
            The method to apply to edges.
        autoname: bool
            Whether to automatically name objects by the ligand name and state
            label.
        autoname_prefix: str
            Prefix for the autonaming; only used if autonaming is ``True``.
        other_components:
            Additional non-alchemical components; keyword will be the string
            label for the component.
        """
        components = {"protein": protein, "solvent": solvent, **other_components}
        leg_labels = {
            "solvent": ["ligand", "solvent"],
            "complex": (["ligand", "solvent", "protein"] + list(other_components)),
        }
        return self._to_rfe_alchemical_network(
            components=components,
            leg_labels=leg_labels,
            protocol=protocol,
            autoname=autoname,
            autoname_prefix=autoname_prefix,
        )

    # on hold until we figure out how to best hack in the PME/NoCutoff
    # switch
    # def to_rhfe_alchemical_network(self, *, solvent, protocol,
    #                                autoname=True,
    #                                autoname_prefix="easy_rhfe",
    #                                **other_components):
    #     leg_labels = {
    #         "solvent": ["ligand", "solvent"] + list(other_components),
    #         "vacuum": ["ligand"] + list(other_components),
    #     }
    #     return self._to_rfe_alchemical_network(
    #         components={"solvent": solvent, **other_components},
    #         leg_labels=leg_labels,
    #         protocol=protocol,
    #         autoname=autoname,
    #         autoname_prefix=autoname_prefix
    #     )

    def is_connected(self) -> bool:
        """Indicates whether all ligands in the network are (directly or indirectly)
        connected to each other.

        A ``False`` value indicates that either some nodes have no edges or
        that there are separate networks that do not link to each other.

        """
        return nx.is_weakly_connected(self.graph)
