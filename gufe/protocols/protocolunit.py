# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""

import abc
from os import PathLike
from typing import Iterable, List, Dict, Any, Optional

from .results import (
    ProtocolUnitResult, ProtocolUnitCompletion, ProtocolUnitFailure,
)


class ProtocolUnit(abc.ABC):
    """A unit of work computable by"""

    def __init__(
        self,
        settings: Optional["ProtocolSettings"] = None,  # type: ignore
        name: Optional[str] = None,
        **kwargs
    ):

        self._settings = settings
        self._kwargs = kwargs
        self._name = name

    def __hash__(self):
        return hash(
            (self.__class__.__name__, self._settings, frozenset(self._kwargs.items()))
        )

    @property
    def settings(self):
        return self._settings

    @property
    def name(self):
        return self._name

    def execute(
        self, dependency_results: Iterable[ProtocolUnitResult], block=True
    ) -> ProtocolUnitResult:
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        dependency_results : Iterable[ProtocolUnitResult]
            The `ProtocolUnitResult`s from the `ProtocolUnit`s this unit is dependent on.
        block : bool
            If `True`, block until execution completes; otherwise run in its own thread.
        """
        if block:
            try:
                out = self._execute(dependency_results)
                result = ProtocolUnitCompletion(
                    name=self._name, dependencies=dependency_results, **out
                )
            except Exception as e:
                result = ProtocolUnitFailure(
                    name=self._name,
                    depedencies=dependency_results,
                    exception=e,
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
