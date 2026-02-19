# This code is part of gufe and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe
import pytest

pytest.importorskip("ipywidgets")
pytest.importorskip("py3Dmol")

from ipywidgets.widgets.widget_box import VBox
from py3Dmol import view

from gufe.visualization.mapping_visualization import display_mapping_3d, display_mappings_3d


def test_score_mappings_rmsd(stereo_chem_mapping) -> None:
    """
    Currently a smoke test
    """
    v = display_mapping_3d(stereo_chem_mapping)
    assert isinstance(v, view)


def test_view_mapping(stereo_chem_mapping) -> None:
    """
    Currently a smoke test
    """
    view = display_mappings_3d([stereo_chem_mapping, stereo_chem_mapping])
    assert isinstance(view, VBox)
