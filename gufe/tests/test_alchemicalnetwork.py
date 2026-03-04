# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import networkx as nx
import pytest

from gufe import AlchemicalNetwork

from .test_tokenization import GufeTokenizableTestsMixin


class TestAlchemicalNetwork(GufeTokenizableTestsMixin):
    cls = AlchemicalNetwork
    repr = None

    @pytest.fixture
    def instance(self, alchem_network_benzene_variants):
        return alchem_network_benzene_variants

    def test_init(self, alchem_network_benzene_variants):
        alnet = alchem_network_benzene_variants

        assert len(alnet.edges) == 12
        assert len(alnet.nodes) == 14

        # should be two unconnected subgraphs given that we defined no
        # connections between solvent and complex systems
        assert not nx.is_weakly_connected(alnet.graph)
        assert nx.number_weakly_connected_components(alnet.graph) == 2

    def test_hashable(self, alchem_network_benzene_variants):
        {alchem_network_benzene_variants}

    def test_connectivity(self, alchem_network_benzene_variants):
        alnet = alchem_network_benzene_variants

        # test connectivity from benzene nodes
        for node in alnet.nodes:
            if node.name == "benzene-solvent":
                edges = alnet.graph.edges(node)
                assert len(edges) == 6
            elif node.name == "benzene-complex":
                edges = alnet.graph.edges(node)
                assert len(edges) == 6
            else:
                edges = alnet.graph.edges(node)
                assert len(edges) == 0

    def test_connected_subgraphs_multiple_subgraphs(self, alchem_network_benzene_variants):
        """Identify two separate networks and one floating nodes as subgraphs."""
        # remove an edge to create a network w/ two subnetworks and one floating node
        edge_list = [e for e in alchem_network_benzene_variants.edges]
        alnet = alchem_network_benzene_variants.copy_with_replacements(edges=edge_list[:-1])

        subgraphs = [subgraph for subgraph in alnet.connected_subgraphs()]

        assert {len(subgraph.nodes) for subgraph in subgraphs} == {6, 7, 1}

        # which graph has the removed node is not deterministic, so we just
        # check that one graph is all-solvent and the other is all-protein
        for subgraph in subgraphs:
            components = [frozenset(n.components.keys()) for n in subgraph.nodes]
            if {"solvent", "protein", "ligand"} in components:
                assert set(components) == {frozenset({"solvent", "protein", "ligand"})}
            else:
                assert set(components) == {frozenset({"solvent", "ligand"})}

    def test_connected_subgraphs_one_subgraph(self, alchem_network_benzene_variants_solvent_only):
        """Return the same network if it only contains one connected component."""
        alnet = alchem_network_benzene_variants_solvent_only
        subgraphs = [subgraph for subgraph in alnet.connected_subgraphs()]
        assert subgraphs == [alnet]
