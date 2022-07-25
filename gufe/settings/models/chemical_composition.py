"""
Pydantic models used for storing settings for chemical Composition.
"""

from .abstract import SettingsBaseModel
from pydantic import BaseModel
from typing import List, Optional

# abstract entity
class chemical_entity(SettingsBaseModel):
    exists: bool = False


# Solutes:
class solute(chemical_entity):
    num_compounds: int = 0


# Protein:
class protein_cofactors(chemical_entity):
    metals: bool = False
    organics: bool = False


class protein(chemical_entity):
    num_residues: int = 0
    biological_assembly: bool = False
    cofactors: protein_cofactors = protein_cofactors()


# Solvent:
class solvent(chemical_entity):
    num_solvent_molecules: int = 0
    num_ions: int = 0
    ion_types: List[str] = []
    solvent_type: str = None


# Membrane:
class membrane(chemical_entity):
    num_lipids: int = 0
    lipid_composition: List[str] = []


# Collecting Chem Comp
class chemical_composition(SettingsBaseModel):
    # What does this look like?
    # We take in smiles or file?
    # Or do we take in a SmallMoleculeComponent object?
    solute: solute = solute()
    protein: protein = protein()
    solvent: solvent = solvent()
    membrane: membrane = membrane()
