import pytest

from gufe import SolventComponent


def test_defaults():
    s = SolventComponent()

    assert s.smiles == 'O'
    assert s.ions == tuple()
    assert s.concentration is None


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
    s = SolventComponent(ions=('Na', 'Cl'))

    assert s != 42


def test_to_dict():
    s = SolventComponent(ions=['Na', 'Cl'])

    assert s.to_dict() == {'smiles': 'O', 'ions': ('Cl', 'Na'),
                           'concentration': None}
