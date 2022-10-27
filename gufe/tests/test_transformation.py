# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest

import io
import pathlib

from gufe.tests.test_tokenization import GufeTokenizableTestsMixin
from gufe.transformations import Transformation, NonTransformation
from gufe.transformations.transformation import ensure_filelike
from gufe.protocols.protocoldag import execute

from .test_protocol import DummyProtocol, DummyProtocolResult
from .test_tokenization import GufeTokenizableTestsMixin


@pytest.fixture
def absolute_transformation(solvated_ligand, solvated_complex):
    return Transformation(
        solvated_ligand,
        solvated_complex,
        protocol=DummyProtocol(settings=None),
        mapping=None,
    )


@pytest.fixture
def complex_equilibrium(solvated_complex):
    return NonTransformation(solvated_complex, protocol=DummyProtocol(settings=None))


@pytest.mark.parametrize('input_type', [
    "str", "path", "TextIO", "BytesIO", "StringIO"
])
def test_ensure_filelike(input_type, tmp_path):
    path = tmp_path / "foo.txt"
    # we choose to use bytes for pathlib.Path just to mix things up;
    # string filename or path can be either bytes or string, so we give one
    # to each
    use_bytes = input_type in {'path', 'BytesIO'}
    filelike = input_type not in {'str', 'path'}
    dumper = {
        'str': str(path),
        'path': path,
        'TextIO': open(path, mode='w'),
        'BytesIO': open(path, mode='wb'),
        'StringIO': io.StringIO(),
    }[input_type]

    if filelike:
        write_mode, read_mode = None, None
    else:
        if use_bytes:
            write_mode, read_mode = "wb", "rb"
        else:
            write_mode, read_mode = "w", "r"

    written = b"bar" if use_bytes else "bar"
    with ensure_filelike(dumper, mode=write_mode) as write_f:
        write_f.write(written)
        write_f.flush()

    if input_type == 'StringIO':
        dumper.seek(0)

    loader = {
        'str': str(path),
        'path': path,
        'TextIO': open(path, mode='r'),
        'BytesIO': open(path, mode='rb'),
        'StringIO': dumper,
    }[input_type]

    with ensure_filelike(loader, mode=read_mode) as read_f:
        loaded = read_f.read()

    assert loaded == written

    # we close pathlikes; do not close filelikes (by default)
    assert write_f.closed is (not filelike)
    assert read_f.closed is (not filelike)

    # ensure everything is closed before we finish
    write_f.close()
    read_f.close()

@pytest.mark.parametrize("input_type", ["TextIO", "BytesIO", "StringIO"])
def test_ensure_filelike_force_close(input_type, tmp_path):
    path = tmp_path / "foo.txt"
    dumper = {
        'TextIO': open(path, mode='w'),
        'BytesIO': open(path, mode='wb'),
        'StringIO': io.StringIO(),
    }[input_type]
    written = b"foo" if input_type == "BytesIO" else "foo"

    with ensure_filelike(dumper, force_close=True) as f:
        f.write(written)

    assert f.closed

@pytest.mark.parametrize("input_type", ["TextIO", "BytesIO", "StringIO"])
def test_ensure_filelike_mode_warning(input_type, tmp_path):
    path = tmp_path / "foo.txt"
    dumper = {
        'TextIO': open(path, mode='w'),
        'BytesIO': open(path, mode='wb'),
        'StringIO': io.StringIO(),
    }[input_type]

    with pytest.warns(UserWarning,
                      match="User-specified mode will be ignored"):
        _ = ensure_filelike(dumper, mode="w")

    dumper.close()

def test_ensure_filelike_default_mode():
    path = "foo.txt"
    loader = ensure_filelike(path)
    assert loader.mode == 'r'


class TestTransformation(GufeTokenizableTestsMixin):

    cls = Transformation
    key = "Transformation-eb839ae80d905adee778c32f516944ad"

    @pytest.fixture
    def instance(self, absolute_transformation):
        return absolute_transformation

    def test_init(self, absolute_transformation, solvated_ligand, solvated_complex):
        tnf = absolute_transformation

        assert tnf.stateA is solvated_ligand
        assert tnf.stateB is solvated_complex

    def test_protocol(self, absolute_transformation):
        tnf = absolute_transformation

        assert isinstance(tnf.protocol, DummyProtocol)

        protocoldag = tnf.create()
        protocoldagresult = execute(protocoldag)

        protocolresult = tnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        len(protocolresult.data) == 1

    def test_equality(self, absolute_transformation, solvated_ligand, solvated_complex):

        opposite = Transformation(
            solvated_complex, solvated_ligand, protocol=DummyProtocol(settings=None)
        )
        assert absolute_transformation != opposite

        different_protocol_settings = Transformation(
            solvated_ligand,
            solvated_complex,
            protocol=DummyProtocol(settings={"lol": True}),
        )
        assert absolute_transformation != different_protocol_settings

        identical = Transformation(
            solvated_ligand,
            solvated_complex,
            protocol=DummyProtocol(settings=None),
            mapping=None,
        )
        assert absolute_transformation == identical

    def test_dump_load_roundtrip(self, absolute_transformation):
        string = io.StringIO()
        absolute_transformation.dump(string)
        string.seek(0)
        recreated = Transformation.load(string)
        assert absolute_transformation == recreated


class TestNonTransformation(GufeTokenizableTestsMixin):

    cls = NonTransformation
    key = "NonTransformation-bc37c512fe411b4dbd38533c7233a5f3"

    @pytest.fixture
    def instance(self, complex_equilibrium):
        return complex_equilibrium

    def test_init(self, complex_equilibrium, solvated_complex):

        ntnf = complex_equilibrium

        assert ntnf.system is solvated_complex

    def test_protocol(self, complex_equilibrium):
        ntnf = complex_equilibrium

        assert isinstance(ntnf.protocol, DummyProtocol)

        protocoldag = ntnf.create()
        protocoldagresult = execute(protocoldag)

        protocolresult = ntnf.gather([protocoldagresult])

        assert isinstance(protocolresult, DummyProtocolResult)

        len(protocolresult.data) == 1

    def test_equality(self, complex_equilibrium, solvated_ligand, solvated_complex):

        different_protocol_settings = NonTransformation(
            solvated_complex, protocol=DummyProtocol(settings={"lol": True})
        )
        assert complex_equilibrium != different_protocol_settings

        identical = NonTransformation(
            solvated_complex, protocol=DummyProtocol(settings=None)
        )
        assert complex_equilibrium == identical

        different_system = NonTransformation(
            solvated_ligand, protocol=DummyProtocol(settings=None)
        )
        assert complex_equilibrium != different_system

    def test_dict_roundtrip(self):
        # TODO: need registration of `Protocol`s for this to work
        ...
