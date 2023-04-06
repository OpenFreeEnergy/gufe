# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
import networkx as nx

from gufe import AlchemicalNetwork, ChemicalSystem, Transformation

from .test_protocol import DummyProtocol, DummyProtocolResult
from .test_tokenization import GufeTokenizableTestsMixin


class TestAlchemicalNetwork(GufeTokenizableTestsMixin):

    cls = AlchemicalNetwork
    key = "AlchemicalNetwork-c065435095312e88fd0cdb834c6d9b49"
    repr = "<AlchemicalNetwork-c065435095312e88fd0cdb834c6d9b49>"

    @pytest.fixture
    def instance(self, benzene_variants_star_map):
        return benzene_variants_star_map

    def test_init(self, benzene_variants_star_map):
        alnet = benzene_variants_star_map

        assert len(alnet.edges) == 12
        assert len(alnet.nodes) == 14

        # should be two unconnected subgraphs given that we defined no
        # connections between solvent and complex systems
        assert not nx.is_weakly_connected(alnet.graph)
        assert nx.number_weakly_connected_components(alnet.graph) == 2

    def test_hashable(self, benzene_variants_star_map):
        {benzene_variants_star_map}

    def test_connectivity(self, benzene_variants_star_map):
        alnet = benzene_variants_star_map

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
