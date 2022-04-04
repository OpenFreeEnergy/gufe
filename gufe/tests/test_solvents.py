import pytest

from gufe import SolventComponent
from openff.units import unit


def test_defaults():
    s = SolventComponent()

    assert s.smiles == 'O'
    assert s.ions is None
    assert s.ion_concentration is None


@pytest.mark.parametrize('other,', [
    # test: ordering, charge dropping, case sensitivity
    ('Cl', 'Na'), ('Na+', 'Cl-'), ('cl', 'na'),
])
def test_hash(other):
    s1 = SolventComponent(ions=('Na', 'Cl'))
    s2 = SolventComponent(ions=other)

    assert s1 == s2
    assert hash(s1) == hash(s2)


def test_neq():
    s1 = SolventComponent(ions=('Na', 'Cl'))
    s2 = SolventComponent(ions=('K', 'Cl'))

    assert s1 != s2


def test_to_dict():
    s = SolventComponent(ions=('Na', 'Cl'))

    assert s.to_dict() == {'smiles': 'O', 'ions': ('Cl', 'Na'),
                           'neutralize': True,
                           'ion_concentration': None}


def test_from_dict():
    s1 = SolventComponent(ions=('Na', 'Cl'),
                          ion_concentration=1.75 * unit.molar,
                          neutralize=False)

    assert SolventComponent.from_dict(s1.to_dict()) == s1


def test_conc():
    s = SolventComponent(ions=('Na', 'Cl'),
                         ion_concentration=1.75 * unit.molar)

    assert s.ion_concentration == unit.Quantity('1.75 M')


@pytest.mark.parametrize('conc,',
                         [1.22,  # no units, 1.22 what?
                          1.5 * unit.kg])  # probably a tad much salt
def test_bad_conc(conc):
    with pytest.raises(ValueError):
        _ = SolventComponent(ions=('Na', 'Cl'), ion_concentration=conc)


def test_solvent_charge():
    s = SolventComponent(ions=('Na', 'Cl'),
                         ion_concentration=1.75 * unit.molar)

    assert s.total_charge is None


@pytest.mark.parametrize('ions,', [
    ('Na',),
    ('Cl', 'Cl'),
    ('F', 'I')
])
def test_bad_inputs(ions):
    print(ions)
    with pytest.raises(ValueError):
        _ = SolventComponent(ions=ions)
