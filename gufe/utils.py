# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import bz2
import functools
import gzip
import io
import lzma
import warnings
from collections.abc import Callable, Generator
from contextlib import ExitStack, contextmanager
from os import PathLike
from typing import BinaryIO, TextIO

# Ordered longest-to-shortest so more-specific signatures win.
_MAGIC_SIGNATURES = [
    (b"\xfd\x37\x7a\x58\x5a\x00", lzma.open),  # xz
    (b"\x42\x5a\x68", bz2.open),  # bzip2
    (b"\x1f\x8b", gzip.open),  # gzip
]

# Keep track of how many bytes we need to read
_HEADER_READ_SIZE = max(len(magic) for magic, _ in _MAGIC_SIGNATURES)


def _detect_opener(header: bytes):
    """Return the decompression opener for *header*, or None for plain text."""
    for magic, opener in _MAGIC_SIGNATURES:
        if header.startswith(magic):
            return opener
    return None


def _detach_safely(stream: io.TextIOWrapper) -> None:
    """Detach a wrapper without closing the caller-owned binary stream."""
    try:
        stream.detach()
    except ValueError:
        # Already closed/detached.
        pass


@contextmanager
def open_text_stream(
    path_or_stream: str | PathLike | BinaryIO | TextIO,
) -> Generator[TextIO, None, None]:
    """Open a file path or stream as text.

    Paths and binary streams are inspected for gzip/bzip2/lzma magic bytes.
    Text streams are yielded unchanged and are assumed to already be decoded
    (and decompressed, if applicable).

    Caller-owned streams are never closed.
    """
    if isinstance(path_or_stream, io.TextIOBase):
        yield path_or_stream
        return

    with ExitStack() as stack:
        if isinstance(path_or_stream, (str, PathLike)):
            raw = stack.enter_context(open(path_or_stream, "rb"))
            raw_is_caller_owned = False
        else:
            raw = path_or_stream
            raw_is_caller_owned = True

        header = raw.read(_HEADER_READ_SIZE)
        if not isinstance(header, (bytes, bytearray, memoryview)):
            raise TypeError(
                "Binary streams passed to open_text_stream must return bytes. "
                "Text streams should be passed as text and will be yielded unchanged."
            )
        header = bytes(header)

        if raw.seekable():
            raw.seek(0)
            source = raw
            source_is_caller_owned = raw_is_caller_owned
        else:
            source = io.BytesIO(header + raw.read())
            source_is_caller_owned = False
            stack.callback(source.close)

        opener = _detect_opener(header)

        if opener is not None:
            stream = opener(source, "rt")
            stack.callback(stream.close)
        else:
            stream = io.TextIOWrapper(source)
            if source_is_caller_owned:
                stack.callback(_detach_safely, stream)
            else:
                stack.callback(stream.close)

        yield stream


class ensure_filelike:
    """Context manager to convert pathlike or filelike to filelike.

    This makes it so that methods can allow a range of user inputs.

    Parameters
    ----------
    fn : PathLike or FileLike
        The user input to normalize.
    mode : str or None
        The mode, if ``fn`` is pathlike. If ``fn`` is filelike, a warning
        will be emitted and the mode will be ignored.
    force_close : bool, default False
        Whether to forcibly close the stream on exit. For pathlike inputs,
        the stream will always be closed. Filelike inputs will close
        if this parameter is True.
    """

    def __init__(self, fn, mode=None, force_close=False):
        filelikes = (io.TextIOBase, io.RawIOBase, io.BufferedIOBase)
        if isinstance(fn, filelikes):
            if mode is not None:
                warnings.warn(
                    f"mode='{mode}' specified with {fn.__class__.__name__}. User-specified mode will be ignored."
                )
            self.to_open = None
            self.do_close = force_close
            self.context = fn
        else:
            if mode is None:
                mode = "r"
            self.to_open = fn
            self.do_close = True
            self.context = None

        self.mode = mode

    def __enter__(self):
        if self.to_open is not None:
            self.context = open(self.to_open, mode=self.mode)

        return self.context

    def __exit__(self, type, value, traceback):
        if self.do_close:
            self.context.close()


# taken from openfe who shamelessly borrowed from openff.toolkit
def requires_package(package_name: str) -> Callable:
    """
    Helper function to denote that a function requires some optional
    dependency. A function decorated with this decorator will raise
    ``MissingDependencyError`` if the package is not found by
    ``importlib.import_module()``.

    Parameters
    ----------
    package_name : str
        The directory path to enter within the context
    Raises
    ------
    MissingDependencyError
    """

    def test_import_for_require_package(function: Callable) -> Callable:
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            import importlib

            try:
                importlib.import_module(package_name)
            except (ImportError, ModuleNotFoundError):
                raise ImportError(function.__name__ + " requires package: " + package_name)
            except Exception as e:
                raise e

            return function(*args, **kwargs)

        return wrapper

    return test_import_for_require_package
