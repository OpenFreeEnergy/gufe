# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
from __future__ import annotations

import json
from typing import Any

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem

from gufe.components import SmallMoleculeComponent
from gufe.visualization.mapping_visualization import draw_mapping

from ..tokenization import JSON_HANDLER
from . import AtomMapping


class LigandAtomMapping(AtomMapping):
    """
    Container for an atom mapping between two small molecule components.

    This is a specialized version of :class:`.AtomMapping` for
    :class:`.SmallMoleculeComponent` which stores the mapping as a dict of
    integers.
    """

    componentA: SmallMoleculeComponent
    componentB: SmallMoleculeComponent
    _annotations: dict[str, Any]
    _compA_to_compB: dict[int, int]

    def __init__(
        self,
        componentA: SmallMoleculeComponent,
        componentB: SmallMoleculeComponent,
        componentA_to_componentB: dict[int, int],
        annotations: dict[str, Any] | None = None,
    ):
        """
        Parameters
        ----------
        componentA, componentB : SmallMoleculeComponent
          the ligand molecules on either end of the mapping
        componentA_to_componentB : dict[int, int]
          correspondence of indices of atoms between the two ligands; the
          keys are indices in componentA and the values are indices in
          componentB.
          These are checked that they are within the possible indices of the
          respective components.
        annotations : dict[str, Any]
          Mapping of annotation identifier to annotation data. Annotations may
          contain arbitrary JSON-serializable data. Annotation identifiers
          starting with ``ofe-`` may have special meaning in other parts of
          OpenFE. ``score`` is a reserved annotation identifier.
        """
        super().__init__(componentA, componentB)

        # validate compA_to_compB
        nA = self.componentA.to_rdkit().GetNumAtoms()
        nB = self.componentB.to_rdkit().GetNumAtoms()
        for i, j in componentA_to_componentB.items():
            if not (0 <= i < nA):
                raise ValueError(f"Got invalid index for ComponentA ({i}); " f"must be 0 <= n < {nA}")
            if not (0 <= j < nB):
                raise ValueError(f"Got invalid index for ComponentB ({i}); " f"must be 0 <= n < {nB}")

        self._compA_to_compB = componentA_to_componentB

        if annotations is None:
            annotations = {}

        self._annotations = annotations

    @classmethod
    def _defaults(cls):
        return {}

    @property
    def componentA_to_componentB(self) -> dict[int, int]:
        return dict(self._compA_to_compB)

    @property
    def componentB_to_componentA(self) -> dict[int, int]:
        return {v: k for k, v in self._compA_to_compB.items()}

    @property
    def componentA_unique(self):
        return (i for i in range(self.componentA.to_rdkit().GetNumAtoms()) if i not in self._compA_to_compB)

    @property
    def componentB_unique(self):
        return (i for i in range(self.componentB.to_rdkit().GetNumAtoms()) if i not in self._compA_to_compB.values())

    @property
    def annotations(self):
        """Any extra metadata, for example the score of a mapping"""
        # return a copy (including copy of nested)
        return json.loads(
            json.dumps(self._annotations, cls=JSON_HANDLER.encoder),
            cls=JSON_HANDLER.decoder,
        )

    def _to_dict(self):
        """Serialize to dict"""
        return {
            "componentA": self.componentA,
            "componentB": self.componentB,
            "componentA_to_componentB": self._compA_to_compB,
            "annotations": json.dumps(self._annotations, sort_keys=True, cls=JSON_HANDLER.encoder),
        }

    @classmethod
    def _from_dict(cls, d: dict):
        """Deserialize from dict"""
        # the mapping dict gets mangled sometimes
        mapping = d["componentA_to_componentB"]
        fixed = {int(k): int(v) for k, v in mapping.items()}

        return cls(
            componentA=d["componentA"],
            componentB=d["componentB"],
            componentA_to_componentB=fixed,
            annotations=json.loads(d["annotations"], cls=JSON_HANDLER.decoder),
        )

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(componentA={self.componentA!r}, "
            f"componentB={self.componentB!r}, "
            f"componentA_to_componentB={self._compA_to_compB!r}, "
            f"annotations={self.annotations!r})"
        )

    def _ipython_display_(self, d2d=None):  # pragma: no-cover
        """
        Visualize atom mapping in a Jupyter Notebook.

        Parameters
        ---------
        d2d : :class:`rdkit.Chem.Draw.rdMolDraw2D.MolDraw2D`
            If desired specify an instance of a MolDraw2D object.
            Default ``None`` will use the MolDraw2DCairo backend.

        Returns
        -------
        Image: IPython.core.display.Image
            Image of the atom map
        """
        from IPython.display import Image, display

        return display(
            Image(
                draw_mapping(
                    self._compA_to_compB,
                    self.componentA.to_rdkit(),
                    self.componentB.to_rdkit(),
                    d2d,
                )
            )
        )

    def draw_to_file(self, fname: str, d2d=None):
        """
        Save atom map visualization to disk

        Parameters
        ---------
        d2d : :class:`rdkit.Chem.Draw.rdMolDraw2D.MolDraw2D`
            If desired specify an instance of a MolDraw2D object.
            Default ``None`` will write a .png file using the MolDraw2DCairo
            backend.

        fname : str
            Name of file to save atom map
        """
        data = draw_mapping(
            self._compA_to_compB,
            self.componentA.to_rdkit(),
            self.componentB.to_rdkit(),
            d2d,
        )
        if type(data) == bytes:
            mode = "wb"
        else:
            mode = "w"

        with open(fname, mode) as f:
            f.write(
                draw_mapping(
                    self._compA_to_compB,
                    self.componentA.to_rdkit(),
                    self.componentB.to_rdkit(),
                    d2d,
                )
            )

    def with_annotations(self, annotations: dict[str, Any]) -> LigandAtomMapping:
        """Create a new mapping based on this one with extra annotations.

        Parameters
        ----------
        annotations : dict[str, Any]
            Annotation update for this mapping. New annotation keys will be
            added to the annotations dict; existing keys will be replaced by
            the data provided here.
        """
        return self.__class__(
            componentA=self.componentA,
            componentB=self.componentB,
            componentA_to_componentB=self._compA_to_compB,
            annotations=dict(**self.annotations, **annotations),
        )

    def get_distances(self) -> NDArray[np.float64]:
        """Return the distances between pairs of atoms in the mapping"""
        dists = []
        molA = self.componentA.to_rdkit().GetConformer()
        molB = self.componentB.to_rdkit().GetConformer()
        for i, j in self._compA_to_compB.items():
            dA = molA.GetAtomPosition(i)
            dB = molB.GetAtomPosition(j)
            dists.append(dA.Distance(dB))

        return np.array(dists)

    def get_alchemical_charge_difference(self) -> int:
        """
        Return the difference in formal charge between stateA and stateB defined as (formal charge A - formal charge B)

        Parameters
        ----------
        mapping: LigandAtomMapping
            The mapping between the end states A and B.

        Returns
        -------
        int:
            The difference in formal charge between the end states.
        """
        molA = self.componentA.to_rdkit()
        molB = self.componentB.to_rdkit()

        charge_a = Chem.rdmolops.GetFormalCharge(molA)
        charge_b = Chem.rdmolops.GetFormalCharge(molB)
        return charge_a - charge_b
