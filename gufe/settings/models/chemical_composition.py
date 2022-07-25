"""
Pydantic models used for storing settings for chemical Composition.
"""

from .abstract import SettingsBaseModel
from typing import List, Optional

# abstract entity
class Chemical_Entity(SettingsBaseModel):
    exists: bool = False


# Solutes:
class Solute(Chemical_Entity):
    num_compounds: int = 0


# Protein:
class Protein_Cofactors(Chemical_Entity):
    metals: bool = False
    organics: bool = False


class Protein(Chemical_Entity):
    num_residues: int = 0
    biological_assembly: bool = False
    cofactors: Protein_Cofactors = Protein_Cofactors()


# Solvent:
class Solvent(Chemical_Entity):
    num_solvent_molecules: int = 0
    num_ions: int = 0
    ion_types: List[str] = []
    solvent_type: str = None


# Membrane:
class Membrane(Chemical_Entity):
    num_lipids: int = 0
    lipid_composition: List[str] = []


# Collecting Chem Comp
class Chemical_Composition(SettingsBaseModel):
    # What does this look like?
    # We take in smiles or file?
    # Or do we take in a SmallMoleculeComponent object?
    solute: Solute = solute()
    protein: Protein = protein()
    solvent: Solvent = solvent()
    membrane: Membrane = membrane()
