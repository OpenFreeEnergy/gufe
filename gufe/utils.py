# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import bz2
import functools
import gzip
import io
import lzma
import warnings
from collections.abc import Callable, Iterator
from contextlib import ExitStack, contextmanager
from os import PathLike
from typing import BinaryIO, TextIO, cast

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
def magic_open(
    path_or_stream: str | PathLike | BinaryIO | TextIO,
) -> Iterator[TextIO]:
    """
    Open a file path or stream as text, transparently decompressing when possible.

    File paths and binary streams are inspected for compression using magic bytes rather than file extensions.
    gzip, bzip2, and lzma/xz are supported.
    Text streams are yielded unchanged and are assumed to already be decoded and, if applicable, already decompressed.

    Parameters
    ----------
    path_or_stream : str or os.PathLike or typing.BinaryIO or typing.TextIO input source to open as text.

        - If a file path is provided, the file is opened in binary mode and inspected for compression.
        - If a binary stream is provided, its leading bytes are inspected for compression.
        - If a text stream is provided, it is yielded unchanged without any compression detection.

    Yields
    ------
    TextIO
        A text-mode stream suitable for iteration and line-based reading.

        - For compressed inputs, this is a text wrapper around the appropriate decompressor.
        - For uncompressed binary inputs, this is an ``io.TextIOWrapper``.
        - For text-mode inputs, this is the original stream.

    Raises
    ------
    TypeError
        Raised if a non-text stream is provided whose ``read()`` method does not return a bytes-like object.

    Notes
    -----
    Compression detection is only possible for file paths and binary streams.
    For text streams, decoding has already occurred, so the original byte signature is no longer available.

    Seekable binary streams are rewound after header inspection so the full content remains available to the decompressor or text wrapper.

    Non-seekable binary streams are fully buffered into an ``io.BytesIO`` object so that the header bytes consumed during detection are preserved.

    Streams opened by this function are closed when the context manager exits.
    Caller-owned streams are never closed by this function.

    Examples
    --------
    Open a plain text file by path:

    >>> with magic_open("data.txt") as f:
    ...     first_line = f.readline()

    Open a compressed file by path:

    >>> with magic_open("data.txt.gz") as f:
    ...     text = f.read()

    Pass an already-open binary stream:

    >>> with open("data.txt.gz", "rb") as raw:
    ...     with magic_open(raw) as f:
    ...         text = f.read()

    Pass an already-open text stream:

    >>> with open("data.txt", "rt") as f:
    ...     with magic_open(f) as text_stream:
    ...         first_line = text_stream.readline()
    """

    if isinstance(path_or_stream, io.TextIOBase):
        # mypy doesn't treat io.TextIOBase and typing.TextIO as the same thing
        # cast JUST tell's mypy "this is TextIO type" IT DOES NOT (repeat) DOES NOT
        # change the type of path_or_stream, just returns it
        yield cast(TextIO, path_or_stream)
        return

    with ExitStack() as stack:
        raw: BinaryIO
        if isinstance(path_or_stream, (str, PathLike)):
            raw = stack.enter_context(open(path_or_stream, "rb"))
            raw_is_caller_owned = False
        else:
            # cast tells mypy that path_or_stream is a BinaryIO object
            # it DOES NOT (repeat) DOES NOT change the type of path_or_stream
            # it only is for type checking and doesn't modify path_or_stream at
            # run time, just returns it
            raw = cast(BinaryIO, path_or_stream)
            raw_is_caller_owned = True

        header = raw.read(_HEADER_READ_SIZE)
        if not isinstance(header, (bytes, bytearray, memoryview)):
            raise TypeError(
                "Binary streams passed to magic_open must return bytes. "
                "Text streams should be passed as text and will be yielded unchanged."
            )
        header = bytes(header)

        # If we can re-wind the stream to 0 after reading in the header
        # we do. This saves us from having to load the whole file into
        # memory. This branch doesn't affect the ownership of the source
        if raw.seekable():
            raw.seek(0)
            source = raw
            source_is_caller_owned = raw_is_caller_owned
        # If we can't re-wind the stream to 0, then we read in the whole
        # stream and concatenate it with the header we already read. This
        # changes the ownership of the stream to "us" so we want to make
        # sure that we close the stream when we are done
        else:
            source = io.BytesIO(header + raw.read())
            source_is_caller_owned = False
            stack.callback(source.close)

        opener = _detect_opener(header)

        stream: TextIO
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
