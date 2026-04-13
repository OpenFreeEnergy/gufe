# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import bz2
import functools
import gzip
import io
import lzma
import warnings
from collections.abc import Callable
from contextlib import nullcontext
from os import PathLike
from typing import IO, TextIO

# Magic bytes for compression format detection
# Reference: https://en.wikipedia.org/wiki/List_of_file_signatures
MAGIC_BYTES = {
    b"\x1f\x8b": lambda f: gzip.open(f, "rt"),  # gzip
    b"\x42\x5a": lambda f: bz2.open(f, "rt"),  # bzip2
    b"\xfd\x37": lambda f: lzma.open(f, "rt"),  # xz/lzma
}


def open_text_stream(path_or_stream: str | PathLike | IO) -> TextIO:
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
    TextIO
        Should be used as a context manager to ensure the stream is properly closed.
        If a file path was provided, the underlying file will be closed on exit.
        If an already-open stream was provided, it will not be closed on exit.

    Raises
    ------
    ValueError
        If a text-mode stream is passed.
        Streams must be opened in binary mode to allow magic byte inspection.

    Notes
    -----
    For seekable streams, the stream is rewound to the start after reading the magic bytes, avoiding loading the file into memory.
    For non-seekable streams, the entire content is buffered into a ``io.BytesIO`` object.
    """
    if isinstance(path_or_stream, (str, PathLike)):
        f = open(path_or_stream, "rb")
        own_file = True
    else:
        f = path_or_stream
        own_file = False
        if hasattr(f, "mode") and "b" not in f.mode:
            raise ValueError(
                "Streams must be opened in binary mode ('rb'), not text mode. open_text_stream will handle decoding."
            )
        if isinstance(f, io.TextIOBase):
            raise ValueError(
                "Streams must be opened in binary mode ('rb'), not text mode. open_text_stream will handle decoding."
            )

    # read the header bytes, then reconstruct the stream so the
    # parser sees the full content regardless of seekability
    header = f.read(2)

    if f.seekable():
        f.seek(0)
        if own_file:
            buffered = f
        else:
            # copy into BytesIO so closing our stream never touches the caller's handle
            buffered = io.BytesIO(f.read())
    else:
        remainder = f.read()
        buffered = io.BytesIO(header + remainder)

    # Check to see if we need to decompress
    # If we do, then opener will be the function we need
    # If we don't, then opener will be None
    opener = MAGIC_BYTES.get(header)

    if opener:
        stream = opener(buffered)
    elif isinstance(buffered, (io.RawIOBase, io.BufferedIOBase)):
        # only wrap in TextIOWrapper if it's still a binary stream
        stream = io.TextIOWrapper(buffered)
    else:
        # already a text stream, pass through as-is
        stream = buffered

    return stream if own_file else nullcontext(stream)


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
