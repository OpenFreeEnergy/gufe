# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import bz2
import functools
import gzip
import io
import lzma
import warnings
from collections.abc import Callable
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
        An open text stream, decompressed if necessary.
    """
    if isinstance(path_or_stream, (str, PathLike)):
        f = open(path_or_stream, "rb")
    else:
        f = path_or_stream
        if hasattr(f, "mode") and "b" not in f.mode:
            raise ValueError(
                "Streams must be opened in binary mode ('rb'), not text mode. open_text_stream will handle decoding."
            )

    # read the header bytes, then reconstruct the stream so the
    # parser sees the full content regardless of seekability
    header = f.read(2)

    if f.seekable():
        f.seek(0)
        buffered = f
    else:
        remainder = f.read()
        buffered = io.BytesIO(header + remainder)

    opener = MAGIC_BYTES.get(header)
    if opener:
        return opener(buffered)

    # wrap in TextIOWrapper if still binary
    if isinstance(buffered, (io.RawIOBase, io.BufferedIOBase)):
        return io.TextIOWrapper(buffered)
    return buffered


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
