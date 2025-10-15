from typing import Dict, Tuple, Union

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D

try:
    import py3Dmol
    from matplotlib import pyplot as plt
    from matplotlib.colors import rgb2hex
except ImportError:
    pass  # Don't throw  error, will happen later

from . import AtomMapping


def _get_max_dist_in_x(atom_mapping: AtomMapping) -> float:
    """helper function
        find the correct mol shift, so no overlap happens in vis

    Returns
    -------
    float
        maximal size of mol in x dimension
    """
    posA = atom_mapping.componentA.to_rdkit().GetConformer().GetPositions()
    posB = atom_mapping.componentB.to_rdkit().GetConformer().GetPositions()
    max_d = []

    for pos in [posA, posB]:
        d = np.zeros(shape=(len(pos), len(pos)))
        for i, pA in enumerate(pos):
            for j, pB in enumerate(pos[i:], start=i):
                d[i, j] = (pB - pA)[0]

        max_d.append(np.max(d))

    estm = float(np.round(max(max_d), 1))
    return estm if (estm > 5) else 5


def _translate(mol, shift: Union[Tuple[float, float, float], NDArray[np.float64]]):
    """
        shifts the molecule by the shift vector

    Parameters
    ----------
    mol : Chem.Mol
        rdkit mol that get shifted
    shift : Tuple[float, float, float]
        shift vector

    Returns
    -------
    Chem.Mol
        shifted Molecule (copy of original one)
    """
    mol = Chem.Mol(mol)
    conf = mol.GetConformer()
    for i, atom in enumerate(mol.GetAtoms()):
        x, y, z = conf.GetAtomPosition(i)
        point = Point3D(x + shift[0], y + shift[1], z + shift[2])
        conf.SetAtomPosition(i, point)
    return mol


def _add_spheres(view: py3Dmol.view, mol1: Chem.Mol, mol2: Chem.Mol, mapping: Dict[int, int]):
    """
        will add spheres according to mapping to the view. (inplace!)

    Parameters
    ----------
    view : py3Dmol.view
        view to be edited
    mol1 : Chem.Mol
        molecule 1 of the mapping
    mol2 : Chem.Mol
        molecule 2 of the mapping
    mapping : Dict[int, int]
        mapping of atoms from mol1 to mol2
    """
    # Get colourmap of size mapping
    cmap = plt.get_cmap("hsv", len(mapping))
    for i, pair in enumerate(mapping.items()):
        p1 = mol1.GetConformer().GetAtomPosition(pair[0])
        p2 = mol2.GetConformer().GetAtomPosition(pair[1])
        color = rgb2hex(cmap(i))
        view.addSphere(
            {
                "center": {"x": p1.x, "y": p1.y, "z": p1.z},
                "radius": 0.6,
                "color": color,
                "alpha": 0.8,
            }
        )
        view.addSphere(
            {
                "center": {"x": p2.x, "y": p2.y, "z": p2.z},
                "radius": 0.6,
                "color": color,
                "alpha": 0.8,
            }
        )
