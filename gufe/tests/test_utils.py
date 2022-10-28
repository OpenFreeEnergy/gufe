# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import pytest

import io
import pathlib

from gufe.utils import ensure_filelike


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
