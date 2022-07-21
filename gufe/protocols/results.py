# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import uuid
from typing import List, Optional, Any, Iterable, Dict

import networkx as nx

from pydantic import BaseModel, PrivateAttr

from .base import ProtocolUnitKey, ProtocolUnitMixin


class ProtocolUnitResultBase(BaseModel, ProtocolUnitMixin):
    class Config:
        arbitrary_types_allowed = True
        allow_mutation = False

    inputs: Dict[str, Any]

    name: Optional[str]      # name of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    pure: bool               # whether `ProtocolUnit` that produced this `ProtocolUnitResult` was a function purely of its inputs
    key: ProtocolUnitKey     # key of the `ProtocolUnit` that produced this `ProtocolUnitResult`

    outputs: Dict[str, Any]  # outputs is a dict returned by a `ProtocolUnit`'s `_execute` method


class ProtocolUnitResult(ProtocolUnitResultBase):
    """Result for a single `ProtocolUnit` execution.

    Immutable upon creation.

    """
    def __init__(self, *, name=None, key, pure, inputs, outputs):
            
        inputs = self._keyencode_dependencies(inputs, ProtocolUnitResult)
        super().__init__(name=name, key=key, pure=pure, inputs=inputs, outputs=outputs)

    def ok(self) -> bool:
        return True


class ProtocolUnitFailure(ProtocolUnitResultBase):
    exception: Exception

    def ok(self) -> bool:
        return False


class ProtocolDAGResult(BaseModel):
    """Result for a single `ProtocolDAG` execution.

    There may be many of these in a given `ResultStore` for a given `Transformation`.
    Data elements from these objects are combined by `Protocol.gather` into a
    `ProtocolResult`.

    Attributes
    ----------
    name : str
        Unique identifier for this `ProtocolDAGResult`.
    graph : nx.DiGraph
        The `ProtocolUnit`s, with dependencies set, as a networkx `DiGraph`.
        Each `ProtocolUnit` features its `ProtocolUnitCompletion` as a `result` attribute.

    """

    class Config:
        arbitrary_types_allowed = True
        allow_mutation = False

    graph: nx.DiGraph
    name: Optional[str]

    @property
    def protocol_units(self):
        return [pu for pu in self.graph.nodes]

    @property
    def protocol_unit_results(self):
        return list(nx.get_node_attributes(self.graph, "result").values())

    def ok(self) -> bool:
        return True


class ProtocolDAGFailure(ProtocolDAGResult):

    def ok(self) -> bool:
        return False

    @property
    def protocol_unit_failures(self):
        return [r for r in nx.get_node_attributes(self.graph, "result").values() if not r.ok()]

