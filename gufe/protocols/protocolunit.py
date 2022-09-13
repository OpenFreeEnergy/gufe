# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""
from __future__ import annotations


import abc
from dataclasses import dataclass
import sys
import traceback
import uuid
from os import PathLike
from pathlib import Path
from copy import copy
from typing import Iterable, Tuple, List, Dict, Any, Optional, Union
import tempfile

from ..tokenization import GufeTokenizable, GufeKey, normalize


@dataclass
class Context:
    """Data class for passing around execution context components to
    `ProtocolUnit._execute`.

    """
    scratch: PathLike
    shared: PathLike


class ProtocolUnitResultBase(GufeTokenizable):
    def __init__(self, *, 
            name: Optional[str] = None, 
            source_key: GufeKey, 
            dependencies: Dict[str, GufeKey],
            outputs: Dict[str, Any],
            _key: GufeKey = None
        ):
        """Generate a `ProtocolUnitResult`.

        Parameters
        ----------
        name : Optional[str]
            Name of the `ProtocolUnit` that produced this `ProtocolUnitResult`.
        source_key : GufeKey
            Key of the `ProtocolUnit` that produced this `ProtocolUnitResult`
        dependencies: Dict[str, GufeKey]
            The key of any `ProtocolUnitResult`s this `ProtocolUnitResult` was dependent on.
        outputs : Dict[str, Any]
            Outputs from the `ProtocolUnit._execute` that generated this
            `ProtocolUnitResult`.
        _key : GufeKey
            Used by deserialization to set UUID-based key for this
            `ProtocolUnitResult` before creation.
        """
        if _key is not None:
            self._key = _key
            
        self._name = name
        self._source_key = source_key
        self._outputs = outputs
        self._dependencies = dependencies

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        # tokenize with uuid
        return uuid.uuid4()

    def _defaults(self):
        return {}

    def _to_dict(self):
        return {'name': self.name,
                '_key': self.key,
                'source_key': self.source_key,
                'dependencies': self.dependencies,
                'outputs': self.outputs}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

    @property
    def name(self):
        return self._name

    @property
    def source_key(self):
        return self._source_key

    @property
    def outputs(self):
        return self._outputs

    @property
    def dependencies(self) -> Dict[str, GufeKey]:
        return self._dependencies


class ProtocolUnitResult(ProtocolUnitResultBase):
    """Result for a single `ProtocolUnit` execution.

    Attributes
    ----------
    name : Optional[str]
        Name of the `ProtocolUnit` that produced this `ProtocolUnitResult`.
    source_key : GufeKey
        Key of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    inputs : Dict[str, Any]
        Inputs to the `ProtocolUnit` that produced this
        `ProtocolUnitResult`. Includes any `ProtocolUnitResult`s this
        `ProtocolUnitResult` is dependent on.
    outputs : Dict[str, GufeKey]
        Outputs from the `ProtocolUnit._execute` that generated this
        `ProtocolUnitResult`.
    dependencies : List[ProtocolUnitResult]
        A list of the `ProtocolUnitResult`s depended upon.

    """

    def ok(self) -> bool:
        return True


class ProtocolUnitFailure(ProtocolUnitResultBase):
    """Failed result for a single `ProtocolUnit` execution.

    Attributes
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
    dependencies : List[ProtocolUnitResult]
        A list of the `ProtocolUnitResult`s depended upon.
    exception : Tuple[str, Tuple[Any, ...]]
        A tuple giving details on the exception raised during `ProtocolUnit`
        execution. The first element gives the type of exception raised; the
        second element is a tuple giving the exception's `args` values.
    traceback : str
        The traceback given by the exception.

    """

    def __init__(self, *, name=None, source_key, dependencies, outputs, _key=None, exception, traceback):
        self._exception = exception
        self._traceback = traceback
        super().__init__(name=name, source_key=source_key, dependencies=dependencies, outputs=outputs, _key=_key)

    def _to_dict(self):
        dct = super()._to_dict()
        dct.update({'exception': self.exception,
                    'traceback': self.traceback})
        return dct

    @property
    def exception(self) -> Tuple[str, Tuple[Any, ...]]:
        return self._exception

    @property
    def traceback(self) -> str:
        return self._traceback

    def ok(self) -> bool:
        return False


class ProtocolUnit(GufeTokenizable):
    """A unit of work within a ProtocolDAG.

    Attributes
    ----------
    name : Optional[str]
        Optional name for the `ProtocolUnit`.
    inputs : Dict[str, Any]
        Inputs to the `ProtocolUnit`. Includes any `ProtocolUnit`s this
        `ProtocolUnit` is dependent on.

    """

    def __init__(
        self,
        *,
        name: Optional[str] = None,
        dependencies: Optional[Dict[str, ProtocolUnit]] = None,
        **inputs,
    ):
        """Create an instance of a ProtocolUnit.

        Parameters
        ----------
        name : str
            Custom name to give this
        dependencies : dict of ProtocolUnits
            the preceeding ProtocolUnit(s) that this Unit requires results from.  The `.key`
            attribute is taken from each and is available via the `.dependencies` attribute
        **inputs 
            Other named arguments, which must all be json-serializable.
        """
        self._name = name
        self._inputs = inputs

        if dependencies is None:
            dependencies = {}
        # Convert dependencies to the key that identifies them
        self._dependencies = {k: v.key for k, v in dependencies.items()}

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        # tokenize with uuid
        return uuid.uuid4()

    def _defaults(self):
        # not used by `ProtocolUnit`s
        return {}

    def _to_dict(self):
        return {'inputs': self.inputs,
                'dependencies': self.dependencies,
                'name': self.name,
                '_key': self.key}

    @classmethod
    def _from_dict(cls, dct: Dict):
        _key = dct.pop('_key')

        obj = cls(name=dct['name'],
                  **dct['inputs'])
        # deps are already tokenized, so hack on
        obj._dependencies = dct['dependencies']
        obj._key = _key

        return obj

    @property
    def name(self):
        return self._name

    @property
    def inputs(self):
        return copy(self._inputs)

    @property
    def dependencies(self) -> Dict[str, GufeKey]:
        """A mapping of {identifier: key} for each dependency

        The key is the `.key` attribute for the original ProtocolUnit
        """
        return self._dependencies     # type: ignore

    def execute(self, *, 
            shared: PathLike, 
            scratch: Optional[PathLike] = None,
            dependencies: Dict[str, ProtocolUnitResult]) -> Union[ProtocolUnitResult, ProtocolUnitFailure]:
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        shared : PathLike
           Path to scratch space that persists across whole DAG execution, but
           is removed after. Used by some `ProtocolUnit`s to pass file contents
           to dependent `ProtocolUnit`s.
        scratch : Optional[PathLike]
            Path to scratch space that persists during execution of this
            `ProtocolUnit`, but removed after.
        dependencies : Dict[str, ProtocolUnitResult]
            Keyword arguments giving the named inputs to `_execute`.
            These can include `ProtocolUnitResult`s from `ProtocolUnit`s this
            unit is dependent on.
        """
        result: Union[ProtocolUnitResult, ProtocolUnitFailure]

        if scratch is None:
            scratch_tmp = tempfile.TemporaryDirectory()
            scratch_ = Path(scratch_tmp.name)
        else:
            scratch_ = Path(scratch)

        context = Context(shared=shared,
                          scratch=scratch_)

        try:
            outputs = self._execute(context, **dependencies)
            result = ProtocolUnitResult(
                name=self.name,
                source_key=self.key,
                dependencies=dependencies,
                outputs=outputs
            )

        except Exception as e:
            result = ProtocolUnitFailure(
                name=self._name,
                source_key=self.key,
                inputs=dependencies,
                outputs=dict(),
                exception=(str(type(e)), e.args),
                traceback=traceback.format_exc()
            )

        # TODO: change this part once we have clearer ideas on how to inject
        # persistent storage use
        if scratch is None:
            scratch_tmp.cleanup()

        return result

    @abc.abstractmethod
    def _execute(self, ctx: Context, **dependencies) -> Dict[str, Any]:
        """Method to override in custom `ProtocolUnit` subclasses.

        A `Context` is always given as its first argument, which provides execution
        context components like filesystem scratch paths.

        Next dependencies are provided, as ProtocolUnitResults, corresponding to the
        `dependencies` `ProtocolUnit` arguments given on instantiation.

        An example of a subclass implementation signature might be:

        >>> class MyProtocolUnit(ProtocolUnit):
        >>>     @staticmethod
        >>>     def _execute(ctx, *, initialization: ProtocolUnitResult):
        >>>        ...

        where instantiation with the subclass `MyProtocolUnit` would look like:

        >>> unit = MyProtocolUnit(settings=settings_dict, 
                                  some_arg=7,
                                  another_arg="five",
                                  dependencies={'initialization': init})

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
