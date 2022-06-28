import pytest
import json
import hashlib

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


@pytest.mark.parametrize('conc', [None, 1.75 * unit.molar])
def test_from_dict(conc):
    s1 = SolventComponent(positive_ion='Na', negative_ion='Cl',
                          ion_concentration=conc,
                          neutralize=False)

    assert SolventComponent.from_dict(s1.to_dict()) == s1


def test_conc():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl',
                         ion_concentration=1.75 * unit.molar)

    assert s.ion_concentration == unit.Quantity('1.75 M')


@pytest.mark.parametrize('conc',
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


def test_to_storage_ready():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl',
                         ion_concentration=1.75 * unit.molar)
    results = s.to_storage_ready()
    assert len(results) == 1
    assert set(results) == {s}
    info = results[s]
    assert info.metadata[":module:"] == "gufe.solventcomponent"
    assert info.metadata[":class:"] == "SolventComponent"
    as_dict = json.loads(info.bytes_data.decode("utf-8"))

    # check that the bytes_data contains all necessary information
    assert SolventComponent.from_dict(as_dict) == s

    # check that we got the right hash
    expected_default_keys = ['smiles', 'neutralize']
    expected_md5 = hashlib.md5(
        json.dumps(
            {k: v for k, v in as_dict.items()
             if k not in expected_default_keys}
        ).encode('utf-8')
    ).hexdigest()
    assert info.md5 == expected_md5

    assert info.path == f"setup/components/{expected_md5[:10]}.json"


def test_from_storage_bytes():
    s = SolventComponent(positive_ion='Na', negative_ion='Cl',
                         ion_concentration=1.75 * unit.molar)
    serialized_bytes = s.to_storage_ready()[s].bytes_data
    loaded = SolventComponent.from_storage_bytes(serialized_bytes,
                                                 lambda x: None)
    assert loaded == s