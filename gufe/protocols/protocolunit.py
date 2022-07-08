# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""

import abc
import uuid
from os import PathLike
from copy import copy
from typing import Iterable, List, Dict, Any, Optional

from dask.base import tokenize, normalize_token

from .base import ProtocolUnitKey, ProtocolUnitMixin
from .results import ProtocolUnitResult


class ProtocolUnit(abc.ABC, ProtocolUnitMixin):
    """A unit of work within a ProtocolDAG."""

    def __init__(
        self,
        *,
        name: Optional[str] = None,
        pure: bool = False,
        **inputs
    ):
        """Create an instance of a ProtocolUnit.

        Parameters
        ----------
        **inputs 
            Keyword arguments, which an include other `ProtocolUnit`s on which this
            `ProtocolUnit` is dependent.

        """
        self._name = name
        self._pure = pure

        # token is generated via hashing when possible, falling back to `uuid4`
        # if not;
        # key is the combination of `<classname>-<token>`, and is used for
        # referencing this object
        self._token = None

        self._key = None

        self._inputs = self._keyencode_dependencies(inputs, ProtocolUnit)

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    #def __hash__(self):
    #    return hash(
    #        (self.__class__.__name__, self._settings, frozenset(self._kwargs.items()))
    #    )

    def _tokenize(self):
        # inspired by implementation of tokenization in dask
        if self._name is not None:
            return self._name
        elif self._pure:
            # tokenize on inputs
            return tokenize(self)
        else:
            # tokenize with uuid
            return uuid.uuid4()

    def __dask_tokenize__(self):
        return list(map(normalize_token, 
            [self._settings,
             self._inputs,
             self._name,
             self._pure]))

    @property
    def settings(self):
        return self._settings

    @property
    def name(self):
        # set name to key if not manually set; used for display
        if self._name is None:
            self._name = self.key()
        return self._name

    @property
    def inputs(self):
        return copy(self._inputs)

    @property
    def pure(self):
        return self._pure

    @property
    def token(self):
        if self._token is None:
            self._token = self._tokenize()
        return self._token

    @property
    def key(self):
        if self._key is None:
            prefix = type(self).__name__
            self._key = ProtocolUnitKey(f"{prefix}-{self.token}")
        return self._key

    def execute(self, block=True, **inputs) -> ProtocolUnitResult:
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        block : bool
            If `True`, block until execution completes; otherwise run in its own thread.

        """

        # process inputs that point to ProtocolUnits

        if block:
            outputs = self._execute(**inputs)
            result = ProtocolUnitResult(
                name=self._name, key=self.key, pure=self.pure, inputs=inputs, outputs=outputs
            )

        else:
            # TODO: wrap in a thread; update status
            ...

        return result

    @abc.abstractmethod
    def _execute(
        self, dependency_results: Iterable[ProtocolUnitResult]
    ) -> Dict[str, Any]:
        ...

    def get_artifacts(self) -> Dict[str, PathLike]:
        """Return a dict of file-like artifacts produced by this
        `ProtocolUnit`.

        """
        ...
