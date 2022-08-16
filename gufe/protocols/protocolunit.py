# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

"""The `ProtocolUnit` class should be subclassed for all units to be used as
part of a `ProtocolDAG`.

"""

import abc
import traceback
import uuid
from os import PathLike
from copy import copy
from typing import Iterable, List, Dict, Any, Optional, Union
import tempfile

from ..tokenize import GufeTokenizable, GufeKey, normalize


class ProtocolUnitResultBase(GufeTokenizable):
    inputs: Dict[str, Any]   # `ProtocolUnitResult`s this `ProtocolUnitResult` is dependent on

    name: Optional[str]      # name of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    pure: bool               # whether `ProtocolUnit` that produced this `ProtocolUnitResult` was a function purely of its inputs
    source_key: GufeKey             # key of the `ProtocolUnit` that produced this `ProtocolUnitResult`

    outputs: Dict[str, Any]  # outputs is a dict returned by a `ProtocolUnit`'s `_execute` method

    def __init__(self, *, 
            name: Optional[str] = None, 
            source_key: GufeKey, 
            pure: bool, 
            inputs: Dict[str, Any], 
            outputs: Dict[str, Any]
        ):
        """Generate a `ProtocolUnitResult`.

        Parameters
        ----------
        name : Optional[str]
            Name of the `ProtocolUnit` that produced this `ProtocolUnitResult`.
        source_key : GufeKey
            Key of the `ProtocolUnit` that produced this `ProtocolUnitResult`
        pure : bool
            If `True`, this `ProtocolUnitResult` is purely a function of the
            inputs that produced it.
        inputs : Dict[str, Any]
            Inputs to the `ProtocolUnit` that produced this
            `ProtocolUnitResult`. Includes any `ProtocolUnitResult`s this
            `ProtocolUnitResult` is dependent on.
        outputs
            Outputs from the `ProtocolUnit._execute` that generated this
            `ProtocolUnitResult`.
        """
            
        self._name = name
        self._source_key = source_key
        self._pure = pure
        self._inputs = inputs
        self._outputs = outputs

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        if self._pure:
            return normalize(self.to_dict(include_defaults=False))
        else:
            # tokenize with uuid
            return uuid.uuid4()

    def _defaults(self):
        # not used by `ProtocolDAG`
        return {}

    def _to_dict(self):
        return {'name': self.name,
                'source_key': self.source_key,
                'pure': self.pure,
                'inputs': self.inputs,
                'outputs': self.outputs}

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(**dct)

    @property
    def name(self):
        """

        """
        # set name to source_key if not manually set; used for display
        if self._name is None:
            self._name = self.source_key
        return self._name

    @property
    def source_key(self):
        return self._source_key

    @property
    def pure(self):
        return self._pure

    @property
    def inputs(self):
        return self._inputs

    @property
    def outputs(self):
        return self._outputs


class ProtocolUnitResult(ProtocolUnitResultBase):
    """Result for a single `ProtocolUnit` execution.

    Immutable upon creation.

    """

    def ok(self) -> bool:
        return True


class ProtocolUnitFailure(ProtocolUnitResultBase):
    exception: Exception

    def __init__(self, *, name=None, source_key, pure, inputs, outputs, exception, traceback):
        self._exception = exception
        self._traceback = traceback
        super().__init__(name=name, source_key=source_key, pure=pure, inputs=inputs, outputs=outputs)

    @property
    def exception(self):
        return self._exception

    @property
    def traceback(self):
        return self._traceback

    def ok(self) -> bool:
        return False


class ProtocolUnit(GufeTokenizable):
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
            `ProtocolUnit` is dependent. For serializability, should be composed of

        """
        self._name = name
        self._pure = pure
        self._inputs = inputs

        self.dag_scratch = None
        self.unit_scratch = None

    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

    def _gufe_tokenize(self):
        if self._pure:
            return normalize(self.to_dict(include_defaults=False))
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

    @classmethod
    def _from_dict(cls, dct: Dict):
        return cls(name=dct['name'],
                   pure=dct['pure'],
                   **dct['inputs'])

    @property
    def name(self):
        # set name to key if not manually set; used for display
        if self._name is None:
            self._name = self.key
        return self._name

    @property
    def inputs(self):
        return copy(self._inputs)

    @property
    def pure(self):
        return self._pure

    def execute(self, *, 
            dag_scratch: PathLike, 
            unit_scratch: PathLike = None, 
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

        self.dag_scratch = dag_scratch
        if unit_scratch is None:
            unit_scratch_tmp = tempfile.TemporaryDirectory()
            self.unit_scratch = unit_scratch_tmp.name
        else:
            self.unit_scratch = unit_scratch

        try:
            outputs = self._execute(**inputs)
            result = ProtocolUnitResult(
                name=self.name, source_key=self.key, pure=self.pure, inputs=inputs, outputs=outputs
            )

        except Exception as e:
            result = ProtocolUnitFailure(
                name=self._name,
                source_key=self.key,
                pure=self.pure,
                inputs=inputs,
                outputs=dict(),
                exception=e,
                traceback=traceback.format_exc(),
            )

        # TODO: change this part once we have clearer ideas on how to inject
        # persistent storage use
        if unit_scratch is None:
            unit_scratch_tmp.cleanup()
        self.unit_scratch = None
        self.dag_scratch = None

        return result

    @abc.abstractmethod
    def _execute(self, **inputs) -> Dict[str, Any]:
        ...
