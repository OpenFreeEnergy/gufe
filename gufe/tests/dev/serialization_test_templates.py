#!/usr/bin/env python
# creates ligand_network.graphml


from rdkit import Chem
from rdkit.Chem import AllChem

from gufe import LigandAtomMapping, LigandNetwork, SmallMoleculeComponent


def mol_from_smiles(smiles: str) -> Chem.Mol:
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)

    return m


# network_template.graphml
mol1 = SmallMoleculeComponent(mol_from_smiles("CCO"))
mol2 = SmallMoleculeComponent(mol_from_smiles("CC"))
mol3 = SmallMoleculeComponent(mol_from_smiles("CO"))

edge12 = LigandAtomMapping(mol1, mol2, {0: 0, 1: 1})
edge23 = LigandAtomMapping(mol2, mol3, {0: 0})
edge13 = LigandAtomMapping(mol1, mol3, {0: 0, 2: 1})

network = LigandNetwork([edge12, edge23, edge13])

with open("ligand_network.graphml", "w") as fn:
    fn.write(network.to_graphml())
