import pytest

from gufe import SolventComponent
from openff.units import unit


def test_defaults():
    s = SolventComponent()

    assert s.smiles == 'O'
    assert s.positive_ion == "Na+"
    assert s.negative_ion == "Cl-"
    assert s.ion_concentration is None


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


def test_to_dict():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl')

    assert s.to_dict() == {'smiles': 'O',
                           'positive_ion': 'Na+',
                           'negative_ion': 'Cl-',
                           'neutralize': True,
                           'ion_concentration': None}


def test_from_dict():
    s1 = SolventComponent(positive_ion='Na', negative_ion='Cl',
                          ion_concentration=1.75 * unit.molar,
                          neutralize=False)

    assert SolventComponent.from_dict(s1.to_dict()) == s1


def test_conc():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl',
                         ion_concentration=1.75 * unit.molar)

    assert s.ion_concentration == unit.Quantity('1.75 M')


@pytest.mark.parametrize('conc,',
                         [1.22,  # no units, 1.22 what?
                          1.5 * unit.kg])  # probably a tad much salt
def test_bad_conc(conc):
    with pytest.raises(ValueError):
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
    with pytest.raises(ValueError):
        _ = SolventComponent(positive_ion=pos, negative_ion=neg)


@pytest.mark.parametrize('pos, neg', [
    ('Na', None), (None, 'Cl'), (None, None),
])
def test_conc_no_ions(pos, neg):
    # if you specify concentration you must also give ions
    with pytest.raises(ValueError):
        _ = SolventComponent(positive_ion=pos, negative_ion=neg,
                             ion_concentration=1.5 * unit.molar)
