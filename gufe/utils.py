# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import bz2
import functools
import gzip
import io
import lzma
import warnings
from collections.abc import Callable, Generator
from contextlib import contextmanager
from os import PathLike
from typing import IO, TextIO

# Ordered longest-to-shortest so more-specific signatures win.
# xz needs 6 bytes, bzip2 needs 3 (BZh), gzip needs 2.
# Reference: https://en.wikipedia.org/wiki/List_of_file_signatures
_MAGIC_SIGNATURES: list[tuple[bytes, object]] = [
    (b"\xfd\x37\x7a\x58\x5a\x00", lzma.open),  # xz  (6 bytes)
    (b"\x42\x5a\x68", bz2.open),  # bzip2 (3 bytes -- BZh)
    (b"\x1f\x8b", gzip.open),  # gzip (2 bytes)
]

# Minimum bytes to read for unambiguous magic-byte detection.
_HEADER_READ_SIZE = 6


def _detect_opener(header: bytes):
    """Return the decompression opener for *header*, or ``None`` for plain text."""
    for magic, opener in _MAGIC_SIGNATURES:
        if header[: len(magic)] == magic:
            return opener
    return None


@contextmanager
def open_text_stream(
    path_or_stream: str | PathLike | IO,
) -> Generator[TextIO, None, None]:
    """Open a file path or stream as text, transparently decompressing if needed.

    Supports gzip, bzip2, and lzma/xz compression, detected via magic bytes rather than file extension.
    Works with both file paths and streams, including non-seekable streams.

    Parameters
    ----------
    path_or_stream
        A file path or an already-open **binary** stream.
        Text streams are not supported as there is no way to inspect the magic bytes after decoding.

    Returns
    -------
    Generator[TextIO, None, None]
        Must be used as a context manager.  If a file path was provided, the
        underlying file will be closed on exit.  If an already-open stream was
        provided, it will not be closed on exit.

    Raises
    ------
    ValueError
        If a text-mode stream is passed.
        Streams must be opened in binary mode to allow magic byte inspection.

    Notes
    -----
    For seekable streams the stream is rewound after the header read so the
    full content is available to the decompressor without loading it into memory.
    For non-seekable streams the entire content is buffered into ``io.BytesIO``.
    Caller-owned streams are never closed when the context manager exits.
    """
    if isinstance(path_or_stream, (str, PathLike)):
        raw = open(path_or_stream, "rb")
        close_raw = True
    else:
        raw = path_or_stream
        close_raw = False
        # if (hasattr(raw, "mode") and "b" not in raw.mode) or isinstance(raw, io.TextIOBase):
        #    raise ValueError(
        #        "Streams must be opened in binary mode ('rb'), not text mode. "
        #        "open_text_stream will handle decoding."
        #    )

    try:
        header = raw.read(_HEADER_READ_SIZE)

        if raw.seekable():
            raw.seek(0)
            source = raw
            # Transfer ownership: clear close_raw so the safety-net finally
            # block does not double-close the handle.
            own_source = close_raw
            close_raw = False
        else:
            # Non-seekable: buffer everything, then close the original handle
            # if we own it so we never leak a file descriptor.
            remainder = raw.read()
            source = io.BytesIO(header + remainder)
            own_source = False  # BytesIO requires no explicit close
            if close_raw:
                raw.close()
                close_raw = False

        opener = _detect_opener(header)

        try:
            if opener is not None:
                # The compression wrappers do *not* close the underlying stream
                # when given a file object, so we remain responsible for
                # closing `source` if we own it.
                stream = opener(source, "rt")
            else:
                # Plain binary stream -- wrap in a text decoder.
                stream = io.TextIOWrapper(source)
        except BaseException:
            if own_source and not source.closed:
                source.close()
            raise

        try:
            yield stream
        finally:
            if opener is not None:
                stream.close()
                if own_source and not source.closed:
                    source.close()
            elif own_source:
                # TextIOWrapper.close() also closes the underlying source.
                stream.close()
            else:
                # Don't own the source: detach so TextIOWrapper doesn't
                # close the caller's handle on garbage collection.
                stream.detach()
    finally:
        # Safety net: close the file we opened if an exception fired before
        # ownership was transferred to `source`.
        if close_raw and not raw.closed:
            raw.close()


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
