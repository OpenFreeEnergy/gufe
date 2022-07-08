# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import uuid
from typing import List, Optional, Any, Iterable, Dict

import networkx as nx

from pydantic import BaseModel, PrivateAttr

from .base import ProtocolUnitToken, ProtocolUnitMixin

class ProtocolUnitResult(BaseModel, ProtocolUnitMixin):
    """Result for a single `ProtocolUnit` execution.

    Immutable upon creation.

    """

    class Config:
        arbitrary_types_allowed = True
        allow_mutation = False

    inputs: Dict[str,Any]

    name: Optional[str]      # name of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    pure: bool               # whether `ProtocolUnit` that produced this `ProtocolUnitResult` was a function purely of its inputs
    token: ProtocolUnitToken # token of the `ProtocolUnit` that produced this `ProtocolUnitResult`

    outputs: Any

    def __init__(self, *, name=None, token, pure, inputs, outputs):
            
        inputs = self._tokenize_dependencies(inputs, ProtocolUnitResult)
        super().__init__(name=name, token=token, inputs=inputs, outputs=outputs, )


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
        Each `ProtocolUnit` features its `ProtocolUnitResult` as a `result` attribute.

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
