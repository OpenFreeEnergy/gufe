# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""

import abc
from dataclasses import dataclass
import traceback
import uuid
from os import PathLike
from pathlib import Path
from copy import copy
from typing import Iterable, List, Dict, Any, Optional, Union
import tempfile

from ..tokenization import GufeTokenizable, GufeKey, normalize


@dataclass
class Context:
    """Data class for passing around execution context components to
    `ProtocolUnit._execute`.

    """
    unit_scratch: PathLike
    dag_scratch: PathLike


class ProtocolUnitMixin:

    def _list_dependencies(self, cls):
        deps = [] 
        for key, value in self.inputs.items():
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


class ProtocolUnitResultBase(GufeTokenizable, ProtocolUnitMixin):
    def __init__(self, *, 
            name: Optional[str] = None, 
            source_key: GufeKey, 
            inputs: Dict[str, Any], 
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
        inputs : Dict[str, Any]
            Inputs to the `ProtocolUnit` that produced this
            `ProtocolUnitResult`. Includes any `ProtocolUnitResult`s this
            `ProtocolUnitResult` is dependent on.
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
        self._inputs = inputs
        self._outputs = outputs

        # for caching
        self._dependencies = None

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
                'inputs': self.inputs,
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
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs

    @property
    def dependencies(self) -> List["ProtocolUnitResult"]:
        """Generate a list of all `ProtocolUnitResult`s dependent on.

        """
        if self._dependencies is None:
            self._dependencies = self._list_dependencies(ProtocolUnitResultBase)
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
    outputs : Dict[str, Any]
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
    exception : Exception
        The exception raised during `ProtocolUnit` exectuion.
    traceback : str
        The traceback given by the exception.

    """

    def __init__(self, *, name=None, source_key, inputs, outputs, _key=None, exception, traceback):
        self._exception = exception
        self._traceback = traceback
        super().__init__(name=name, source_key=source_key, inputs=inputs, outputs=outputs, _key=_key)

    def _to_dict(self):
        dct = super()._to_dict()
        dct.update({'exception': self.exception,
                    'traceback': self.traceback})
        return dct

    @property
    def exception(self) -> Exception:
        return self._exception

    @property
    def traceback(self) -> str:
        return self._traceback

    def ok(self) -> bool:
        return False


class ProtocolUnit(GufeTokenizable, ProtocolUnitMixin):
    """A unit of work within a ProtocolDAG.

    Attributes
    ----------
    name : Optional[str]
        Optional name for the `ProtocolUnit`.
    inputs : Dict[str, Any]
        Inputs to the `ProtocolUnit`. Includes any `ProtocolUnit`s this
        `ProtocolUnit` is dependent on.
    dependencies : List[ProtocolUnit]
        A list of the `ProtocolUnit`s depended upon.

    """

    def __init__(
        self,
        *,
        name: Optional[str] = None,
        _key: GufeKey = None,
        **inputs
    ):
        """Create an instance of a ProtocolUnit.

        Parameters
        ----------
        name : str
            Custom name to give this 
        _key : GufeKey
            Used by deserialization to set UUID-based key for this
            `ProtocolUnit` before creation.
        **inputs 
            Keyword arguments, which can include other `ProtocolUnit`s on which
            this `ProtocolUnit` is dependent. For easy serializability, should
            be composed of JSON-serializable types where possible.

        """
        if _key is not None:
            self._key = _key

        self._name = name
        self._inputs = inputs

        # for caching
        self._dependencies = None

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
                'name': self.name,
                '_key': self.key}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(name=dct['name'],
                   _key=dct['_key'],
                   **dct['inputs'])

    @property
    def name(self):
        return self._name

    @property
    def inputs(self):
        return copy(self._inputs)

    @property
    def dependencies(self) -> List["ProtocolUnit"]:
        """Generate a list of all `ProtocolUnit`s dependent on.

        """
        if self._dependencies is None:
            self._dependencies = self._list_dependencies(ProtocolUnit)
        return self._dependencies

    def execute(self, *, 
            dag_scratch: PathLike, 
            unit_scratch: Optional[PathLike] = None, 
            **inputs) -> Union[ProtocolUnitResult, ProtocolUnitFailure]:
        """Given `ProtocolUnitResult`s from dependencies, execute this `ProtocolUnit`.

        Parameters
        ----------
        dag_scratch : PathLike
           Path to scratch space that persists across whole DAG execution, but
           is removed after. Used by some `ProtocolUnit`s to pass file contents
           to dependent `ProtocolUnit`s.
        unit_scratch : Optional[PathLike]
            Path to scratch space that persists during execution of this
            `ProtocolUnit`, but removed after.
            
        **inputs
            Keyword arguments giving the named inputs to `_execute`.
            These can include `ProtocolUnitResult`s from `ProtocolUnit`s this
            unit is dependent on.

        """
        result: Union[ProtocolUnitResult, ProtocolUnitFailure]

        if unit_scratch is None:
            unit_scratch_tmp = tempfile.TemporaryDirectory()
            unit_scratch_ = Path(unit_scratch_tmp.name)
        else:
            unit_scratch_ = Path(unit_scratch)

        context = Context(dag_scratch=dag_scratch,
                          unit_scratch=unit_scratch_)

        try:
            outputs = self._execute(context, **inputs)
            result = ProtocolUnitResult(
                name=self.name, source_key=self.key, inputs=inputs, outputs=outputs
            )

        except Exception as e:
            result = ProtocolUnitFailure(
                name=self._name,
                source_key=self.key,
                inputs=inputs,
                outputs=dict(),
                exception=e,
                traceback=traceback.format_exc(),
            )

        # TODO: change this part once we have clearer ideas on how to inject
        # persistent storage use
        if unit_scratch is None:
            unit_scratch_tmp.cleanup()

        return result

    @staticmethod
    @abc.abstractmethod
    def _execute(ctx: Context, **inputs) -> Dict[str, Any]:
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

        """
        ...
