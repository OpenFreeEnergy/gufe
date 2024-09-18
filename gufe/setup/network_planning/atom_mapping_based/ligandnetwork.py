# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from __future__ import annotations

from itertools import chain
import json
import networkx as nx
from typing import FrozenSet, Iterable, Optional
import gufe

from gufe import SmallMoleculeComponent
from gufe.tokenization import GufeTokenizable

from gufe.setup.network_planning.network_plan import NetworkPlan
from gufe.setup.network_planning.atom_mapping_based.ligand_atom_mapping import LigandAtomMapping


class LigandNetwork(NetworkPlan):
    """A directed graph connecting many ligands according to their atom mapping

    Parameters
    ----------
    edges : Iterable[LigandAtomMapping]
        edges for this network
    nodes : Iterable[SmallMoleculeComponent]
        nodes for this network
    """
    def __init__(
        self,
        edges: Iterable[LigandAtomMapping],
        nodes: Optional[Iterable[SmallMoleculeComponent]] = None
    ):
        if nodes is None:
            nodes = []

        self._edges = frozenset(edges)
        edge_nodes = set(chain.from_iterable((e.componentA, e.componentB)
                                             for e in edges))
        self._nodes = frozenset(edge_nodes) | frozenset(nodes)
        self._graph = None

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self) -> dict:
        return {'graphml': self.to_graphml()}

    @classmethod
    def _from_dict(cls, dct: dict):
        return cls.from_graphml(dct['graphml'])

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
        mol_to_label = {mol: f"mol{num}"
                        for num, mol in enumerate(sorted_nodes)}

        edge_data = sorted([
            (
                mol_to_label[edge.componentA],
                mol_to_label[edge.componentB],
                json.dumps(list(edge.componentA_to_componentB.items()))
            )
            for edge in self.edges
        ])

        # from here, we just build the graph
        serializable_graph = nx.MultiDiGraph()
        for mol, label in mol_to_label.items():
            serializable_graph.add_node(label,
                                        moldict=json.dumps(mol.to_dict(),
                                                           sort_keys=True))

        for molA, molB, mapping in edge_data:
            serializable_graph.add_edge(molA, molB, mapping=mapping)

        return serializable_graph

    @classmethod
    def _from_serializable_graph(cls, graph: nx.Graph):
        """Create network from NetworkX graph with serializable attributes.

        This is the inverse of ``_serializable_graph``.
        """
        label_to_mol = {node: SmallMoleculeComponent.from_dict(json.loads(d))
                        for node, d in graph.nodes(data='moldict')}

        edges = [
            LigandAtomMapping(componentA=label_to_mol[node1],
                              componentB=label_to_mol[node2],
                              componentA_to_componentB=dict(json.loads(mapping)))
            for node1, node2, mapping in graph.edges(data='mapping')
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
            edges = set([])

        if nodes is None:
            nodes = set([])

        return LigandNetwork(self.edges | set(edges), self.nodes | set(nodes))

    def _to_rfe_alchemical_network(
        self,
        components: dict[str, gufe.Component],
        leg_labels: dict[str, list[str]],
        protocol: gufe.Protocol,
        *,
        alchemical_label: str = "ligand",
        autoname=True,
        autoname_prefix=""
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
                    syscomps.update({label: components[label]
                                     for label in other_labels})

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

                transformation = gufe.Transformation(sysA, sysB, protocol,
                                                     mapping=edge,
                                                     name=name)

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
        **other_components
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
        components = {
            'protein': protein,
            'solvent': solvent,
            **other_components
        }
        leg_labels = {
            "solvent": ["ligand", "solvent"],
            "complex": (["ligand", "solvent", "protein"]
                        + list(other_components)),
        }
        return self._to_rfe_alchemical_network(
            components=components,
            leg_labels=leg_labels,
            protocol=protocol,
            autoname=autoname,
            autoname_prefix=autoname_prefix
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