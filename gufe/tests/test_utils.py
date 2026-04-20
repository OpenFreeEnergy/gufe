# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import bz2
import gzip
import io
import lzma

import pytest

from gufe.utils import ensure_filelike, magic_open


@pytest.mark.parametrize("input_type", ["str", "path", "TextIO", "BytesIO", "StringIO"])
def test_ensure_filelike(input_type, tmp_path):
    path = tmp_path / "foo.txt"
    # we choose to use bytes for pathlib.Path just to mix things up;
    # string filename or path can be either bytes or string, so we give one
    # to each
    use_bytes = input_type in {"path", "BytesIO"}
    filelike = input_type not in {"str", "path"}
    dumper = {
        "str": str(path),
        "path": path,
        "TextIO": open(path, mode="w"),
        "BytesIO": open(path, mode="wb"),
        "StringIO": io.StringIO(),
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

    if input_type == "StringIO":
        dumper.seek(0)

    loader = {
        "str": str(path),
        "path": path,
        "TextIO": open(path),
        "BytesIO": open(path, mode="rb"),
        "StringIO": dumper,
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
        "TextIO": open(path, mode="w"),
        "BytesIO": open(path, mode="wb"),
        "StringIO": io.StringIO(),
    }[input_type]
    written = b"foo" if input_type == "BytesIO" else "foo"

    with ensure_filelike(dumper, force_close=True) as f:
        f.write(written)

    assert f.closed


@pytest.mark.parametrize("input_type", ["TextIO", "BytesIO", "StringIO"])
def test_ensure_filelike_mode_warning(input_type, tmp_path):
    path = tmp_path / "foo.txt"
    dumper = {
        "TextIO": open(path, mode="w"),
        "BytesIO": open(path, mode="wb"),
        "StringIO": io.StringIO(),
    }[input_type]

    with pytest.warns(UserWarning, match="User-specified mode will be ignored"):
        _ = ensure_filelike(dumper, mode="w")

    dumper.close()


def test_ensure_filelike_default_mode():
    path = "foo.txt"
    loader = ensure_filelike(path)
    assert loader.mode == "r"


PLAIN_TEXT = b"OpenFE can be quite useful on real projects\n"


@pytest.fixture
def plain_text_path(tmp_path):
    p = tmp_path / "test.txt"
    p.write_bytes(PLAIN_TEXT)
    return p


@pytest.fixture(params=["gz", "bz2", "xz"])
def compressed_path(tmp_path, request):
    suffix = request.param
    p = tmp_path / f"test.txt.{suffix}"
    openers = {"gz": gzip.open, "bz2": bz2.open, "xz": lzma.open}
    with openers[suffix](p, "wb") as f:
        f.write(PLAIN_TEXT)
    return p


class _NonSeekableBytesIO(io.BytesIO):
    def seekable(self) -> bool:
        return False

    def seek(self, *args, **kwargs):
        raise io.UnsupportedOperation("not seekable")


class _BadBinaryLike:
    def read(self, size=-1):
        return "not bytes"

    def seekable(self):
        return True

    def seek(self, offset, whence=0):
        return 0


class TestOpenTextStream:
    def test_plain_text_path(self, plain_text_path):
        with magic_open(plain_text_path) as f:
            assert f.read() == PLAIN_TEXT.decode()

    def test_compressed_path(self, compressed_path):
        """All three compression formats are detected via magic bytes."""
        with magic_open(compressed_path) as f:
            assert f.read() == PLAIN_TEXT.decode()

    def test_magic_bytes_not_extension(self, tmp_path):
        """Compression is detected from magic bytes, not file extension."""
        p = tmp_path / "test.txt"  # no .gz extension
        with gzip.open(p, "wb") as f:
            f.write(PLAIN_TEXT)
        with magic_open(p) as f:
            assert f.read() == PLAIN_TEXT.decode()

    def test_binary_stream(self, plain_text_path):
        with open(plain_text_path, "rb") as binary_stream:
            with magic_open(binary_stream) as f:
                assert f.read() == PLAIN_TEXT.decode()

    def test_compressed_binary_stream(self, compressed_path):
        with open(compressed_path, "rb") as binary_stream:
            with magic_open(binary_stream) as f:
                assert f.read() == PLAIN_TEXT.decode()

    def test_non_seekable_stream(self, compressed_path):
        """Non-seekable streams are buffered into BytesIO."""
        with open(compressed_path, "rb") as f:
            data = f.read()
        with magic_open(_NonSeekableBytesIO(data)) as stream:
            assert stream.read() == PLAIN_TEXT.decode()

    def test_caller_stream_not_closed(self, plain_text_path):
        """When a stream is passed in, it is not closed on context manager exit."""
        with open(plain_text_path, "rb") as binary_stream:
            with magic_open(binary_stream):
                pass
            assert not binary_stream.closed

    def test_path_stream_is_closed(self, plain_text_path):
        """When a path is passed in, the stream is closed on context manager exit."""
        with magic_open(plain_text_path) as f:
            pass
        assert f.closed

    def test_text_stream_is_passthrough(self, plain_text_path):
        with open(plain_text_path) as text_stream:
            text_stream.read(5)

            with magic_open(text_stream) as stream:
                assert stream is text_stream
                assert stream.read() == PLAIN_TEXT.decode()[5:]

            assert not text_stream.closed

    def test_non_seekable_plain_binary_stream(self):
        with magic_open(_NonSeekableBytesIO(PLAIN_TEXT)) as stream:
            assert stream.read() == PLAIN_TEXT.decode()

    def test_compressed_caller_stream_not_closed(self, compressed_path):
        with open(compressed_path, "rb") as binary_stream:
            with magic_open(binary_stream) as stream:
                assert stream.read() == PLAIN_TEXT.decode()
            assert not binary_stream.closed

    def test_non_text_stream_must_return_bytes(self):
        with pytest.raises(TypeError, match="must return bytes"):
            with magic_open(_BadBinaryLike()):
                pass

    def test_already_detached_wrapper_on_exit(self, plain_text_path):
        """Cleanup tolerates a caller-owned plain stream already detached by user code."""
        with open(plain_text_path, "rb") as binary_stream:
            with magic_open(binary_stream) as stream:
                assert isinstance(stream, io.TextIOWrapper)
                assert stream.read() == PLAIN_TEXT.decode()

                detached = stream.detach()
                assert detached is binary_stream

            # The context manager should not fail when cleanup sees an already
            # detached wrapper, and the caller-owned stream must remain open.
            assert not binary_stream.closed
            binary_stream.seek(0)
            assert binary_stream.read() == PLAIN_TEXT
