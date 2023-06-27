# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from typing import Iterable, NamedTuple
import pytest
import importlib.resources
import gufe
from gufe.tests.test_protocol import DummyProtocol
from gufe import SmallMoleculeComponent, LigandNetwork, LigandAtomMapping

from rdkit import Chem

from networkx import NetworkXError

from .test_tokenization import GufeTokenizableTestsMixin


def mol_from_smiles(smi):
    m = Chem.MolFromSmiles(smi)
    m.Compute2DCoords()

    return m

class _NetworkTestContainer(NamedTuple):
    """Container to facilitate network testing"""
    network: LigandNetwork
    nodes: Iterable[SmallMoleculeComponent]
    edges: Iterable[LigandAtomMapping]
    n_nodes: int
    n_edges: int


@pytest.fixture
def ligandnetwork_graphml():
    with importlib.resources.path('gufe.tests.data', 'ligand_network.graphml') as file:
        with open(file, 'r') as f:
            yield f.read()


@pytest.fixture
def mols():
    mol1 = SmallMoleculeComponent(mol_from_smiles("CCO"))
    mol2 = SmallMoleculeComponent(mol_from_smiles("CC"))
    mol3 = SmallMoleculeComponent(mol_from_smiles("CO"))
    return mol1, mol2, mol3


@pytest.fixture
def std_edges(mols):
    mol1, mol2, mol3 = mols
    edge12 = LigandAtomMapping(mol1, mol2, {0: 0, 1: 1})
    edge23 = LigandAtomMapping(mol2, mol3, {0: 0})
    edge13 = LigandAtomMapping(mol1, mol3, {0: 0, 2: 1})
    return edge12, edge23, edge13


@pytest.fixture
def simple_network(mols, std_edges):
    """Network with no edges duplicated and all nodes in edges"""
    network = LigandNetwork(std_edges)
    return _NetworkTestContainer(
        network=network,
        nodes=mols,
        edges=std_edges,
        n_nodes=3,
        n_edges=3,
    )


@pytest.fixture
def doubled_edge_network(mols, std_edges):
    """Network with more than one edge associated with a node pair"""
    mol1, mol2, mol3 = mols
    extra_edge = LigandAtomMapping(mol1, mol3, {0: 0, 1: 1})
    edges = list(std_edges) + [extra_edge]
    return _NetworkTestContainer(
        network=LigandNetwork(edges),
        nodes=mols,
        edges=edges,
        n_nodes=3,
        n_edges=4,
    )


@pytest.fixture
def singleton_node_network(mols, std_edges):
    """Network with nodes not included in any edges"""
    extra_mol = SmallMoleculeComponent(mol_from_smiles("CCC"))
    all_mols = list(mols) + [extra_mol]
    return _NetworkTestContainer(
        network=LigandNetwork(edges=std_edges, nodes=all_mols),
        nodes=all_mols,
        edges=std_edges,
        n_nodes=4,
        n_edges=3,
    )

@pytest.fixture
def real_molecules_network(benzene, phenol, toluene):
    """Small network with full mappings"""
    # benzene to phenol
    bp_mapping = {i: i for i in range(10)}
    bp_mapping.update({10: 12, 11:11})

    # benzene to toluene
    bt_mapping = {i: i + 4 for i in range(10)}
    bt_mapping.update({10: 2, 11: 14})

    network = gufe.LigandNetwork([
        gufe.LigandAtomMapping(benzene, toluene, bt_mapping),
        gufe.LigandAtomMapping(benzene, phenol, bp_mapping),
    ])
    return network

@pytest.fixture(params=['simple', 'doubled_edge', 'singleton_node'])
def network_container(
    request,
    simple_network,
    doubled_edge_network,
    singleton_node_network,
):
    """Fixture to allow parameterization of the network test"""
    network_dct = {
        'simple': simple_network,
        'doubled_edge': doubled_edge_network,
        'singleton_node': singleton_node_network,
    }
    return network_dct[request.param]


class TestLigandNetwork(GufeTokenizableTestsMixin):
    cls = LigandNetwork
    key = "LigandNetwork-8d9c3198d7fbfc29e73cb09911bccc7f"
    repr = "<LigandNetwork-8d9c3198d7fbfc29e73cb09911bccc7f>"

    @pytest.fixture
    def instance(self, simple_network):
        return simple_network.network

    def test_node_type(self, network_container):
        n = network_container.network

        assert all((isinstance(node, SmallMoleculeComponent) for node in n.nodes))

    def test_graph(self, network_container):
        # The NetworkX graph that comes from the ``.graph`` property should
        # have nodes and edges that match the Network container object.
        graph = network_container.network.graph
        assert len(graph.nodes) == network_container.n_nodes
        assert set(graph.nodes) == set(network_container.nodes)
        assert len(graph.edges) == network_container.n_edges
        # extract the AtomMappings from the nx edges
        mappings = [
            atommapping for _, _, atommapping in graph.edges.data('object')
        ]
        assert set(mappings) == set(network_container.edges)
        # ensure LigandAtomMapping stored in nx edge is consistent with nx edge
        for mol1, mol2, atommapping in graph.edges.data('object'):
            assert atommapping.componentA == mol1
            assert atommapping.componentB == mol2

    def test_graph_annotations(self, mols, std_edges):
        mol1, mol2, mol3 = mols
        edge12, edge23, edge13 = std_edges
        annotated = edge12.with_annotations({'foo': 'bar'})
        network = LigandNetwork([annotated, edge23, edge13])
        assert network.graph[mol1][mol2][0]['foo'] == 'bar'

    def test_graph_immutability(self, mols, network_container):
        # The NetworkX graph that comes from that ``.graph`` property should
        # be immutable.
        graph = network_container.network.graph
        mol = SmallMoleculeComponent(mol_from_smiles("CCCC"))
        mol_CC = mols[1]
        edge = LigandAtomMapping(mol, mol_CC, {1: 0, 2: 1})
        with pytest.raises(NetworkXError, match="Frozen graph"):
            graph.add_node(mol)
        with pytest.raises(NetworkXError, match="Frozen graph"):
            graph.add_edge(edge)

    def test_nodes(self, network_container):
        # The nodes reported by a ``Network`` should match expectations for
        # that network.
        network = network_container.network
        assert len(network.nodes) == network_container.n_nodes
        assert set(network.nodes) == set(network_container.nodes)

    def test_edges(self, network_container):
        # The edges reported by a ``Network`` should match expectations for
        # that network.
        network = network_container.network
        assert len(network.edges) == network_container.n_edges
        assert set(network.edges) == set(network_container.edges)

    def test_enlarge_graph_add_nodes(self, simple_network):
        # New nodes added via ``enlarge_graph`` should exist in the newly
        # created network
        new_mol = SmallMoleculeComponent(mol_from_smiles("CCCC"))
        network = simple_network.network
        new_network = network.enlarge_graph(nodes=[new_mol])
        assert new_network is not network
        assert new_mol not in network.nodes
        assert new_mol in new_network.nodes
        assert len(new_network.nodes) == len(network.nodes) + 1
        assert set(new_network.edges) == set(network.edges)

    def test_enlarge_graph_add_edges(self, mols, simple_network):
        # New edges added via ``enlarge_graph`` should exist in the newly
        # created network
        mol1, _, mol3 = mols
        extra_edge = LigandAtomMapping(mol1, mol3, {0: 0, 1: 1})
        network = simple_network.network
        new_network = network.enlarge_graph(edges=[extra_edge])
        assert new_network is not network
        assert extra_edge not in network.edges
        assert extra_edge in new_network.edges
        assert len(new_network.edges) == len(network.edges) + 1
        assert set(new_network.nodes) == set(network.nodes)

    def test_enlarge_graph_add_edges_new_nodes(self, mols, simple_network):
        # New nodes included implicitly by edges added in ``enlarge_graph``
        # should exist in the newly created network
        new_mol = SmallMoleculeComponent(mol_from_smiles("CCCC"))
        mol_CC = mols[1]
        extra_edge = LigandAtomMapping(new_mol, mol_CC, {1: 0, 2: 1})
        network = simple_network.network
        new_network = network.enlarge_graph(edges=[extra_edge])
        assert new_network is not network
        assert extra_edge not in network.edges
        assert extra_edge in new_network.edges
        assert new_mol not in network.nodes
        assert new_mol in new_network.nodes
        assert len(new_network.edges) == len(network.edges) + 1
        assert len(new_network.nodes) == len(network.nodes) + 1

    def test_enlarge_graph_add_duplicate_node(self, simple_network):
        # Adding a duplicate of an existing node should create a new network
        # with the same edges and nodes as the previous one.
        duplicate = SmallMoleculeComponent(mol_from_smiles("CC"))
        network = simple_network.network

        existing = network.nodes
        assert duplicate in existing  # matches by ==

        new_network = network.enlarge_graph(nodes=[duplicate])
        assert len(new_network.nodes) == len(network.nodes)
        assert set(new_network.nodes) == set(network.nodes)
        assert len(new_network.edges) == len(network.edges)
        assert set(new_network.edges) == set(network.edges)

    def test_enlarge_graph_add_duplicate_edge(self, mols, simple_network):
        # Adding a duplicate of an existing edge should create a new network
        # with the same edges and nodes as the previous one.
        mol1, _, mol3 = mols
        duplicate = LigandAtomMapping(mol1, mol3, {0: 0, 2: 1})
        network = simple_network.network

        existing = network.edges
        assert duplicate in existing  # matches by ==
        assert any(duplicate is edge for edge in existing)  # one edge *is* the duplicate

        new_network = network.enlarge_graph(edges=[duplicate])
        assert len(new_network.nodes) == len(network.nodes)
        assert set(new_network.nodes) == set(network.nodes)
        assert len(new_network.edges) == len(network.edges)
        assert set(new_network.edges) == set(network.edges)

    def test_serialization_cycle(self, simple_network):
        network = simple_network.network
        serialized = network.to_graphml()
        deserialized = LigandNetwork.from_graphml(serialized)
        reserialized = deserialized.to_graphml()
        assert serialized == reserialized
        assert network == deserialized

    def test_to_graphml(self, simple_network, ligandnetwork_graphml):
        assert simple_network.network.to_graphml() == ligandnetwork_graphml

    def test_from_graphml(self, simple_network, ligandnetwork_graphml):
        assert LigandNetwork.from_graphml(ligandnetwork_graphml) == simple_network.network

    def test_is_connected(self, simple_network):
        assert simple_network.network.is_connected()

    def test_is_not_connected(self, singleton_node_network):
        assert not singleton_node_network.network.is_connected()

    def test_to_rbfe_network(self):
    @pytest.mark.parametrize('with_cofactor', [True, False])
    def test_to_rbfe_alchemical_network(
        self,
        real_molecules_network,
        prot_comp,
        solv_comp,
        request,
        with_cofactor,
    ):
        # obviously, this particular set of ligands with this particular
        # protein makes no sense, but we should still be able to set it up
        if with_cofactor:
            others = {'cofactor': request.getfixturevalue('styrene')}
        else:
            others = {}

        protocol = DummyProtocol(DummyProtocol.default_settings())
        rbfe = real_molecules_network.to_rbfe_alchemical_network(
            solvent=solv_comp,
            protein=prot_comp,
            protocol=protocol,
            **others
        )

        expected_names = {
            'easy_rbfe_benzene_solvent_toluene_solvent',
            'easy_rbfe_benzene_complex_toluene_complex',
            'easy_rbfe_benzene_solvent_phenol_solvent',
            'easy_rbfe_benzene_complex_phenol_complex',
        }
        names = set(edge.name for edge in rbfe.edges)
        assert names == expected_names

        assert len(rbfe.edges) == 2 * len(real_molecules_network.edges)
        for edge in rbfe.edges:
            assert edge.protocol == protocol

            compsA = edge.stateA.components
            compsB = edge.stateB.components
            if 'solvent' in edge.name:
                labels = {'solvent', 'ligand'}
            elif 'complex' in edge.name:
                labels = {'solvent', 'ligand', 'protein'}
                if with_cofactor:
                    labels.add('cofactor')
            else:  # -no-cov-
                raise RuntimeError("Something went weird in testing. Unable "
                                   f"to get leg for edge {edge}")

            assert set(compsA) == labels
            assert set(compsB) == labels

            assert compsA['ligand'] != compsB['ligand']
            assert compsA['ligand'].name == 'benzene'
            assert compsA['solvent'] == compsB['solvent']
            # for things that might not always exist, use .get
            assert compsA.get('protein') == compsB.get('protein')
            assert compsA.get('cofactor') == compsB.get('cofactor')

            assert list(edge.mapping) == ['ligand']
            assert edge.mapping['ligand'] in real_molecules_network.edges


    def test_to_rbfe_alchemical_network_autoname(
        self,
        real_molecules_network,
        prot_comp,
        solv_comp
    ):
        pytest.skip()


    def test_to_rhfe_alchemical_network(self, real_molecules_network,
                                        solv_comp):

        others = {}
        protocol = DummyProtocol(DummyProtocol.default_settings())
        rhfe = real_molecules_network.to_rhfe_alchemical_network(
            solvent=solv_comp,
            protocol=protocol,
            **others
        )

        expected_names = {
            'easy_rhfe_benzene_vacuum_toluene_vacuum',
            'easy_rhfe_benzene_solvent_toluene_solvent',
            'easy_rhfe_benzene_vacuum_phenol_vacuum',
            'easy_rhfe_benzene_solvent_phenol_solvent',
        }
        names = set(edge.name for edge in rhfe.edges)
        assert names == expected_names

        assert len(rhfe.edges) == 2 * len(real_molecules_network.edges)
        for edge in rhfe.edges:
            assert edge.protocol == protocol

            compsA = edge.stateA.components
            compsB = edge.stateB.components

            if 'vacuum' in edge.name:
                labels = {'ligand'}
            elif 'solvent' in edge.name:
                labels = {'ligand', 'solvent'}
            else:  # -no-cov-
                raise RuntimeError("Something went weird in testing. Unable "
                                   f"to get leg for edge {edge}")

            labels |= set(others)

            assert set(compsA) == labels
            assert set(compsB) == labels

            assert compsA['ligand'] != compsB['ligand']
            assert compsA['ligand'].name == 'benzene'
            assert compsA.get('solvent') == compsB.get('solvent')

            assert list(edge.mapping) == ['ligand']
            assert edge.mapping['ligand'] in real_molecules_network.edges


def test_empty_ligand_network(mols):
    # issue #217
    n = LigandNetwork(edges=[], nodes=[mols[0]])

    assert len(n.edges) == 0
    assert len(n.nodes) == 1
