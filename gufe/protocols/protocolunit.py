# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""

import abc
import uuid
from os import PathLike
from copy import copy
from typing import Iterable, List, Dict, Any, Optional, Union

from dask.base import normalize_token

from ..base import GufeTokenizable
from .base import ProtocolUnitMixin
from .results import (
    ProtocolUnitResult, ProtocolUnitFailure,
)


class ProtocolUnit(GufeTokenizable, ProtocolUnitMixin):
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
        name : str
            Custom name to give this 
        **inputs 
            Keyword arguments, which an include other `ProtocolUnit`s on which this
            `ProtocolUnit` is dependent.

        """
        self._name = name
        self._pure = pure

        # we likely want to key-encode dependencies upfront
        # this saves having a full chain of dependencies for computing key for
        # this object; effectively makes key calculation happen 
        # approximately once per `ProtocolUnit`
        # although if the key is cached, no need to recalculate

        # so perhaps best to leave this be for simplicity?
        self._inputs = inputs

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        if self._pure:
            # we use `dask` tokenization components for this
            # allows for deterministic handling of e.g. pandas DataFrames
            return list(map(normalize_token,
                            sorted(self.to_dict(include_defaults=False).items(),
                                   key=str)))
        else:
            # tokenize with uuid
            return uuid.uuid4()

    def _defaults(self):
        # not used by `ProtocolUnit`s
        return {}

    def _to_dict(self):
        return {'inputs': self.inputs,
                'name': self.name,
                'pure': self.pure}

    def _from_dict(cls, dct: Dict):
        cls(name=dct['name'],
            pure=dct['pure'],
            **dct['inputs'])

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

    def execute(self, **inputs) -> Union[ProtocolUnitResult, ProtocolUnitFailure]:
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        **inputs
            Keyword arguments giving the named inputs to `_execute`.
            These can include `ProtocolUnitResult`s from `ProtocolUnit`s this
            unit is dependent on.

        """
        try:
            outputs = self._execute(**inputs)
            result = ProtocolUnitResult(
                name=self._name, key=self.key, pure=self.pure, inputs=inputs, outputs=outputs
            )

        except Exception as e:
            result = ProtocolUnitFailure(
                name=self._name,
                key=self.key,
                pure=self.pure,
                inputs=inputs,
                outputs=dict(),
                exception=e,
            )

        return result

    @abc.abstractmethod
    def _execute(self, **inputs) -> Dict[str, Any]:
        ...
