# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""
from __future__ import annotations

import abc
import datetime
import sys
import tempfile
import traceback
import uuid
from collections.abc import Iterable
from copy import copy
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from ..tokenization import TOKENIZABLE_REGISTRY, GufeKey, GufeTokenizable
from .errors import ExecutionInterrupt


@dataclass
class Context:
    """Data class for passing around execution context components to
    `ProtocolUnit._execute`.

    """

    scratch: PathLike
    shared: PathLike


def _list_dependencies(inputs, cls):
    deps = []
    for key, value in inputs.items():
        if isinstance(value, dict):
            for k, v in value.items():
                if isinstance(v, cls):
                    deps.append(v)
        elif isinstance(value, list):
            for i in value:
                if isinstance(i, cls):
                    deps.append(i)
        elif isinstance(value, cls):
            deps.append(value)
    return deps


class ProtocolUnitResult(GufeTokenizable):
    """
    Successful result of a single :class:`ProtocolUnit` execution.
    """

    def __init__(
        self,
        *,
        name: str | None = None,
        source_key: GufeKey,
        inputs: dict[str, Any],
        outputs: dict[str, Any],
        start_time: datetime.datetime | None = None,
        end_time: datetime.datetime | None = None,
    ):
        """
        Parameters
        ----------
        name : Optional[str]
            Name of the `ProtocolUnit` that produced this `ProtocolUnitResult`.
        source_key : GufeKey
            Key of the `ProtocolUnit` that produced this `ProtocolUnitResult`
        inputs : Dict[str, Any]
            Inputs to the `ProtocolUnit` that produced this
            `ProtocolUnitResult`. Includes any `ProtocolUnitResult`s this
            `ProtocolUnitResult` is dependent on.
        outputs : Dict[str, Any]
            Outputs from the `ProtocolUnit._execute` that generated this
            `ProtocolUnitResult`.
        start_time, end_time: datetime.datetime
            The start and end time for executing this Unit
        """
        self._name = name
        self._source_key = source_key
        self._inputs = inputs
        self._outputs = outputs
        self._start_time = start_time
        self._end_time = end_time

        # for caching
        self._dependencies = None

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        # tokenize with uuid
        return uuid.uuid4().hex

    @classmethod
    def _defaults(cls):
        return {}

    def _to_dict(self):
        return {
            "name": self.name,
            "_key": self.key,
            "source_key": self.source_key,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "start_time": self.start_time,
            "end_time": self.end_time,
        }

    @classmethod
    def _from_dict(cls, dct: dict):
        key = dct.pop("_key")
        obj = cls(**dct)
        obj._set_key(key)
        return obj

    @property
    def name(self):
        return self._name

    @property
    def source_key(self):
        return self._source_key

    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs

    @property
    def dependencies(self) -> list[ProtocolUnitResult]:
        """All results that this result was dependent on"""
        if self._dependencies is None:
            self._dependencies = _list_dependencies(self._inputs, ProtocolUnitResult)
        return self._dependencies  # type: ignore

    @staticmethod
    def ok() -> bool:
        return True

    @property
    def start_time(self) -> datetime.datetime | None:
        """The time execution of this Unit began"""
        return self._start_time

    @property
    def end_time(self) -> datetime.datetime | None:
        """The time at which execution of this Unit finished"""
        return self._end_time


class ProtocolUnitFailure(ProtocolUnitResult):
    """Failed result of a single :class:`ProtocolUnit` execution."""

    def __init__(
        self,
        *,
        name=None,
        source_key,
        inputs,
        outputs,
        _key=None,
        exception,
        traceback,
        start_time: datetime.datetime | None = None,
        end_time: datetime.datetime | None = None,
    ):
        """
        Parameters
        ----------
        name : Optional[str]
            Name of the `ProtocolUnit` that produced this `ProtocolUnitResult`.
        source_key : GufeKey
            Key of the `ProtocolUnit` that produced this `ProtocolUnitResult`
        inputs : Dict[str, Any]
            Inputs to the `ProtocolUnit` that produced this
            `ProtocolUnitResult`. Includes any `ProtocolUnitResult` this
            `ProtocolUnitResult` was dependent on.
        outputs : Dict[str, Any]
            Outputs from the `ProtocolUnit._execute` that generated this
            `ProtocolUnitResult`.
        exception : Tuple[str, Tuple[Any, ...]]
            A tuple giving details on the exception raised during `ProtocolUnit`
            execution. The first element gives the type of exception raised; the
            second element is a tuple giving the exception's `args` values.
        traceback : str
            The traceback given by the exception.
        """
        self._exception = exception
        self._traceback = traceback
        super().__init__(
            name=name,
            source_key=source_key,
            inputs=inputs,
            outputs=outputs,
            start_time=start_time,
            end_time=end_time,
        )

    def _to_dict(self):
        dct = super()._to_dict()
        dct.update({"exception": self.exception, "traceback": self.traceback})
        return dct

    @property
    def exception(self) -> tuple[str, tuple[Any, ...]]:
        return self._exception

    @property
    def traceback(self) -> str:
        return self._traceback

    @staticmethod
    def ok() -> bool:
        return False


class ProtocolUnit(GufeTokenizable):
    """A unit of work within a ProtocolDAG."""

    _dependencies: list[ProtocolUnit] | None

    def __init__(self, *, name: str | None = None, **inputs):
        """Create an instance of a ProtocolUnit.

        Parameters
        ----------
        name : str
            Custom name to give this
        **inputs
            Keyword arguments, which can include other `ProtocolUnit`s on which
            this `ProtocolUnit` is dependent. Should be either `gufe` objects
            or JSON-serializables.

        """
        self._name = name
        self._inputs = inputs

        # for caching
        self._dependencies = None

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        # tokenize with uuid
        return uuid.uuid4().hex

    @classmethod
    def _defaults(cls):
        # not used by `ProtocolUnit`s
        return {}

    def _to_dict(self):
        return {"inputs": self.inputs, "name": self.name, "_key": self.key}

    @classmethod
    def _from_dict(cls, dct: dict):
        _key = dct.pop("_key")

        obj = cls(name=dct["name"], **dct["inputs"])
        obj._set_key(_key)
        return obj

    @property
    def name(self) -> str | None:
        """
        Optional name for the `ProtocolUnit`.
        """
        return self._name

    @property
    def inputs(self) -> dict[str, Any]:
        """
        Inputs to the `ProtocolUnit`.

        Includes any `ProtocolUnit` instances this `ProtocolUnit` depends on.
        """
        return copy(self._inputs)

    @property
    def dependencies(self) -> list[ProtocolUnit]:
        """All units that this unit is dependent on (parents)"""
        if self._dependencies is None:
            self._dependencies = _list_dependencies(self._inputs, ProtocolUnit)
        return self._dependencies  # type: ignore

    def execute(
        self, *, context: Context, raise_error: bool = False, **inputs
    ) -> ProtocolUnitResult | ProtocolUnitFailure:
        """Given `ProtocolUnitResult` s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        context : Context
            Execution context for this `ProtocolUnit`; includes e.g. ``shared``
            and ``scratch`` `Path` s.
        raise_error : bool
            If True, raise any errors instead of catching and returning a
            `ProtocolUnitFailure` default False
        **inputs
            Keyword arguments giving the named inputs to `_execute`.
            These can include `ProtocolUnitResult` objects from `ProtocolUnit`
            objects this unit is dependent on.

        """
        result: ProtocolUnitResult | ProtocolUnitFailure
        start = datetime.datetime.now()

        try:
            outputs = self._execute(context, **inputs)
            result = ProtocolUnitResult(
                name=self.name,
                source_key=self.key,
                inputs=inputs,
                outputs=outputs,
                start_time=start,
                end_time=datetime.datetime.now(),
            )

        except KeyboardInterrupt or ExecutionInterrupt:
            # we always want to raise in these situations
            raise
        except Exception as e:
            if raise_error:
                raise

            result = ProtocolUnitFailure(
                name=self._name,
                source_key=self.key,
                inputs=inputs,
                outputs=dict(),
                exception=(e.__class__.__qualname__, e.args),
                traceback=traceback.format_exc(),
                start_time=start,
                end_time=datetime.datetime.now(),
            )

        return result

    @staticmethod
    @abc.abstractmethod
    def _execute(ctx: Context, **inputs) -> dict[str, Any]:
        """Method to override in custom `ProtocolUnit` subclasses.

        A `Context` is always given as its first argument, which provides execution
        context components like filesystem scratch paths. Every other argument
        is provided as **inputs corresponding to the keyword arguments given to the
        `ProtocolUnit` subclass on instantiation.

        An example of a subclass implementation signature might be:

        >>> class MyProtocolUnit(ProtocolUnit):
        >>>     @staticmethod
        >>>     def _execute(ctx, *, initialization: ProtocolUnitResult, settings, some_arg=2, **inputs):
        >>>        ...

        where instantiation with the subclass `MyProtocolUnit` would look like:

        >>> unit = MyProtocolUnit(settings=settings_dict,
                                  initialization=another_protocolunit,
                                  some_arg=7,
                                  another_arg="five")

        Inside of `_execute` above:
        - `settings`, and `some_arg`, would have their values set as given
        - `initialization` would get the `ProtocolUnitResult` that comes from
          `another_protocolunit`'s own execution.
        - `another_arg` would be accessible via `inputs['another_arg']`

        This allows protocol developers to define how `ProtocolUnit`s are
        chained together, with their outputs exposed to `ProtocolUnit`s
        dependent on them.

        The return dict should contain objects that are either `gufe` objects
        or JSON-serializables.

        """
        ...
