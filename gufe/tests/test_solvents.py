from gufe import SolventComponent


def test_defaults():
    s = SolventComponent()

    assert s.smiles == 'O'
    assert s.ions == tuple()
    assert s.concentration is None


def test_hash():
    s1 = SolventComponent(ions=('Cl-', 'Na+'))
    s2 = SolventComponent(ions=('Na', 'Cl'))

    assert s1 == s2
    assert hash(s1) == hash(s2)


def test_to_dict():
    s = SolventComponent(ions=['Na', 'Cl'])

    assert s.to_dict() == {'smiles': 'O', 'ions': ('Cl', 'Na'),
                           'concentration': None}
