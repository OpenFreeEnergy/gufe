# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import uuid
from typing import List, Optional, Any, Iterable

import networkx as nx

from pydantic import BaseModel, PrivateAttr


class ProtocolUnitResult(BaseModel):
    """Result for a single `ProtocolUnit` execution.

    Immutable upon creation.

    """

    class Config:
        extra = "allow"
        allow_mutation = False

    _dependencies: List[str] = PrivateAttr()
    _uuid: str = PrivateAttr()

    name: Optional[
        str
    ]  # name of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    data: Any  # should likely be fleshed out, currently a free-for-all

    def __init__(
        self, dependencies: Optional[Iterable["ProtocolUnitResult"]] = None, **data
    ):

        if dependencies is None:
            dependencies = []

        super().__init__(**data)

        self._uuid = str(uuid.uuid4())

        dep_uuids = [dep._uuid for dep in dependencies]
        self._dependencies = dep_uuids


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

    graph: nx.DiGraph
    name: Optional[str]

    @property
    def protocol_units(self):
        return [pu for pu in self.graph.nodes]

    @property
    def protocol_unit_results(self):
        return list(nx.get_node_attributes(self.graph, "result").values())
