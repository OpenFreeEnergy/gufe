import pytest

from gufe import SolventComponent
from openff.units import unit

from .test_tokenization import GufeTokenizableTestsMixin


def test_defaults():
    s = SolventComponent()

    assert s.smiles == 'O'
    assert s.positive_ion == "Na+"
    assert s.negative_ion == "Cl-"
    assert s.ion_concentration == 0.0 * unit.molar


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

    assert s.to_dict() == {'__module__': 'gufe.components.solventcomponent',
                           '__qualname__': 'SolventComponent',
                           'smiles': 'O',
                           'positive_ion': 'Na+',
                           'negative_ion': 'Cl-',
                           'neutralize': True,
                           'ion_concentration': '0.0 molar'}


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


class TestSolventComponent(GufeTokenizableTestsMixin):

    cls = SolventComponent
    key = "SolventComponent-187d235ef3c2035d8505083c8ad7d0a0"

    @pytest.fixture
    def instance(self):
        return SolventComponent(positive_ion='Na', negative_ion='Cl')
