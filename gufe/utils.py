# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import io
import warnings

from os import PathLike, rmdir
import pathlib

import logging
_logger = logging.getLogger(__name__)


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
                    f"mode='{mode}' specified with {fn.__class__.__name__}."
                    " User-specified mode will be ignored."
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


def delete_empty_dirs(root: PathLike, delete_root: bool = True):
    """Delete all empty directories.

    Repeats so that directories that only contained empty directories also
    get deleted.
    """
    root = pathlib.Path(root)

    def find_empty_dirs(directory):
        if not (paths := list(directory.iterdir())):
            return [directory]
        directories = [p for p in paths if p.is_dir()]
        return sum([find_empty_dirs(d) for d in directories], [])

    while root.exists() and (empties := find_empty_dirs(root)):
        if empties == [root] and not delete_root:
            return
        for directory in empties:
            _logger.debug(f"Removing '{directory}'")
            rmdir(directory)
