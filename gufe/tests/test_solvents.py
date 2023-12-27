import pytest
import numpy as np
from gufe import SolventComponent
from openff.units import unit

from .test_tokenization import GufeTokenizableTestsMixin


def test_defaults():
    s = SolventComponent()

    assert s.smiles == 'O'
    assert s.positive_ion == "Na+"
    assert s.negative_ion == "Cl-"
    assert s.ion_concentration == 0.15 * unit.molar


@pytest.mark.parametrize('pos, neg', [
    # test: charge dropping, case sensitivity
    ('Na', 'Cl'), ('Na+', 'Cl-'), ('na', 'cl'),
])
def test_hash(pos, neg):
    s1 = SolventComponent(positive_ion='Na', negative_ion='Cl')
    s2 = SolventComponent(positive_ion=pos, negative_ion=neg)

    assert s1 == s2
    assert hash(s1) == hash(s2)
    assert s2.positive_ion == 'Na+'
    assert s2.negative_ion == 'Cl-'


def test_neq():
    s1 = SolventComponent(positive_ion='Na', negative_ion='Cl')
    s2 = SolventComponent(positive_ion='K', negative_ion='Cl')

    assert s1 != s2


@pytest.mark.parametrize('conc', [0.0 * unit.molar, 1.75 * unit.molar])
def test_from_dict(conc):
    s1 = SolventComponent(positive_ion='Na', negative_ion='Cl',
                          ion_concentration=conc,
                          neutralize=False)

    assert SolventComponent.from_dict(s1.to_dict()) == s1


def test_conc():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl',
                         ion_concentration=1.75 * unit.molar)

    assert s.ion_concentration == unit.Quantity('1.75 M')


@pytest.mark.parametrize('conc,',
                         [1.22,  # no units, 1.22 what?
                          1.5 * unit.kg,  # probably a tad much salt
                          -0.1 * unit.molar])  # negative conc
def test_bad_conc(conc):
    errmsg = "ion_concentration must be"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(positive_ion='Na', negative_ion='Cl',
                             ion_concentration=conc)


def test_solvent_charge():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl',
                         ion_concentration=1.75 * unit.molar)

    assert s.total_charge is None


@pytest.mark.parametrize('pos, neg,', [
    ('Na', 'C'),
    ('F', 'I'),
])
def test_bad_inputs(pos, neg):
    errmsg = "Invalid"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(positive_ion=pos, negative_ion=neg)


@pytest.mark.parametrize('nsolv', [0, -42])
def test_negative_num_solvent(nsolv):
    errmsg = "num_solvent must be greater than zero"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(num_solvent=nsolv)


@pytest.mark.parametrize('box_vectors,solvent_padding', [
    [None, 1.2*unit.nanometer],
    [20 * np.identity(3) * unit.angstrom, None],
])
def test_defined_solvent_with_defined_box_error(box_vectors, solvent_padding):
    errmsg = "Cannot define the number of solvent molecules"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            num_solvent=42,
            solvent_density=850*unit.kilogram/unit.meter**3,
            solvent_padding=solvent_padding,
            box_vectors=box_vectors,
        )


def test_solvent_padding_box_vectors_error():
    errmsg = "cannot be defined alongside box_vectors"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            solvent_padding=1.2*unit.nanometer,
            box_vectors=20*np.identity(3)*unit.angstrom,
        )


def test_solvent_padding_not_distance_error():
    errmsg = "solvent_padding must be given in units of"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            solvent_padding=1.2*unit.molar,
        )


def test_negative_solvent_padding_error():
    errmsg = "solvent_padding must be positive"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            solvent_padding=-1*unit.nanometer
        )


def test_incompatible_density_error():
    errmsg = "solvent_density must be given in units compatible with g/mL"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            solvent_density=1*unit.nanometer
        )


def test_negative_density_error():
    errmsg = "solvent_density cannot be negative"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            solvent_density=-850*unit.kilogram/unit.meter**3,
        )


def test_unknown_box_type_error():
    errmsg = "Unknown box_shape passed"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            box_shape='rhombic'
        )


def test_no_units_box_vectors():
    errmsg = "box_vector must be defined as a unit.Quantity"
    with pytest.raises(ValueError, match=errmsg):
        _ = SolventComponent(
            solvent_padding=None,
            box_vectors=20 * np.identity(3),
        )


class TestSolventComponent(GufeTokenizableTestsMixin):

    cls = SolventComponent
    key = "SolventComponent-26b4034ad9dbd9f908dfc298ea8d449f"
    repr = "SolventComponent(name=O, Na+, Cl-, None, 0.15 molar, 1.2 nanometer, 1.2 nanometer, cube, None)"

    @pytest.fixture
    def instance(self):
        return SolventComponent(positive_ion='Na', negative_ion='Cl')
