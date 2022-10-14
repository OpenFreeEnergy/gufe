# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest
import networkx as nx

from gufe import AlchemicalNetwork, ChemicalSystem, Transformation

from .test_protocol import DummyProtocol, DummyProtocolResult
from .test_tokenization import GufeTokenizableTestsMixin


@pytest.fixture
def benzene_variants_star_map(
    benzene,
    toluene,
    phenol,
    benzonitrile,
    anisole,
    benzaldehyde,
    styrene,
    prot_comp,
    solv_comp,
):

    variants = [toluene, phenol, benzonitrile, anisole, benzaldehyde, styrene]

    # define the solvent chemical systems and transformations between benzene and the others
    solvated_ligands = {}
    solvated_ligand_transformations = {}

    solvated_ligands["benzene"] = ChemicalSystem(
        {"solvent": solv_comp, "ligand": benzene}, name="benzene-solvent"
    )

    for ligand in variants:
        solvated_ligands[ligand.name] = ChemicalSystem(
            {"solvent": solv_comp, "ligand": ligand}, name=f"{ligand.name}-solvnet"
        )
        solvated_ligand_transformations[("benzene", ligand.name)] = Transformation(
            solvated_ligands["benzene"],
            solvated_ligands[ligand.name],
            protocol=DummyProtocol(settings=None),
            mapping=None,
        )

    # define the complex chemical systems and transformations between benzene and the others
    solvated_complexes = {}
    solvated_complex_transformations = {}

    solvated_complexes["benzene"] = ChemicalSystem(
        {"protein": prot_comp, "solvent": solv_comp, "ligand": benzene},
        name="benzene-complex",
    )

    for ligand in variants:
        solvated_complexes[ligand.name] = ChemicalSystem(
            {"protein": prot_comp, "solvent": solv_comp, "ligand": ligand},
            name=f"{ligand.name}-complex",
        )
        solvated_complex_transformations[("benzene", ligand.name)] = Transformation(
            solvated_complexes["benzene"],
            solvated_complexes[ligand.name],
            protocol=DummyProtocol(settings=None),
            mapping=None,
        )

    return AlchemicalNetwork(
        list(solvated_ligand_transformations.values())
        + list(solvated_complex_transformations.values())
    )


class TestAlchemicalNetwork(GufeTokenizableTestsMixin):

    cls = AlchemicalNetwork
    key = "AlchemicalNetwork-54aa4f41f09cb07ebe6febdc65f1c278"

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
