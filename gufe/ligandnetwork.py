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
        """NetworkX graph for this network

        This graph will have :class:`.ChemicalSystem` objects as nodes and
        :class:`.Transformation` objects as directed edges
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
        """A read-only view of the edges of the Network"""
        return self._edges

    @property
    def nodes(self) -> frozenset[SmallMoleculeComponent]:
        """A read-only view of the nodes of the Network"""
        return self._nodes

    def _serializable_graph(self) -> nx.Graph:
        """
        Create NetworkX graph with serializable attribute representations.

        This enables us to use easily use different serialization
        approaches.
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
        """Create network from NetworkX graph with serializable attributes.

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
        """Return the GraphML string representing this Network

        This is the primary serialization mechanism for this class.

        Returns
        -------
        str :
            string representing this network in GraphML format
        """
        return "\n".join(nx.generate_graphml(self._serializable_graph()))

    @classmethod
    def from_graphml(cls, graphml_str: str) -> LigandNetwork:
        """Create from a GraphML string.

        Parameters
        ----------
        graphml_str : str
            GraphML string representation of a :class:`.Network`

        Returns
        -------
        LigandNetwork
            new network from the GraphML
        """
        return cls._from_serializable_graph(nx.parse_graphml(graphml_str))

    def enlarge_graph(self, *, edges=None, nodes=None) -> LigandNetwork:
        """
        Create a new network with the given edges and nodes added

        Parameters
        ----------
        edges : Iterable[:class:`.LigandAtomMapping`]
            edges to append to this network
        nodes : Iterable[:class:`.SmallMoleculeComponent`]
            nodes to append to this network

        Returns
        -------
        LigandNetwork
            a new network adding the given edges and nodes to this network
        """
        if edges is None:
            edges = set()

        if nodes is None:
            nodes = set()

        return LigandNetwork(self.edges | set(edges), self.nodes | set(nodes))

    def remove_edges(self, edges: Union[LigandAtomMapping, list[LigandAtomMapping]]) -> LigandNetwork:
        """Create a new copy of this network with some edges removed

        Note that this will not remove any nodes, potentially resulting in
        disconnected networks

        Parameters
        ----------
        edges : list[LigandAtomMapping] or LigandAtomMapping
          the edges to remove, these *must* be present in the network

        Returns
        -------
        network : LigandNetwork
        """
        if isinstance(edges, LigandAtomMapping):
            edges = [edges]

        to_remove = set(edges)
        current = set(self.edges)

        # check that all edges to remove are present
        if extras := to_remove - current:
            raise ValueError("Some edges weren't already present: "
                             f"{extras}")

        new_edges = current - to_remove
new_net = LigandNetwork(new_edges, self.nodes)
if not new_net.graph().is_connected():
    raise ValueError("The result is a disconnected Network!")
return new_net 
        return LigandNetwork(new_edges, self.nodes)

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
        """
        Parameters
        ----------
        components: dict[str, :class:`.Component`]
            non-alchemical components (components that will be on both sides
            of a transformation)
        leg_labels: dict[str, list[str]]
            mapping of the names for legs (the keys of this dict) to a list
            of the component names. The component names must be the same as
            used in the ``components`` dict.
        protocol: :class:`.Protocol`
            the protocol to apply
        alchemical_label: str
            the label for the component undergoing an alchemical
            transformation (default ``'ligand'``)
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
        """Convert the ligand network to an AlchemicalNetwork

        Parameters
        ----------
        protocol: Protocol
            the method to apply to edges
        autoname: bool
            whether to automatically name objects by the ligand name and
            state label
        autoname_prefix: str
            prefix for the autonaming; only used if autonaming is True
        other_components:
            additional non-alchemical components, keyword will be the string
            label for the component
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
        """Are all ligands in the network (indirectly) connected to each other

        A "False" value indicates that either some ligands have no edges or that
        there are separate networks that do not link to each other.
        """
        return nx.is_weakly_connected(self.graph)
