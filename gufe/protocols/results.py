# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/gufe

import abc
import uuid
from typing import List, Optional, Any, Iterable

import networkx as nx

from pydantic import BaseModel, PrivateAttr


class ProtocolUnitResult(BaseModel, abc.ABC):
    """Result for a single `ProtocolUnit` execution.

    Immutable upon creation.

    """
    class Config:
        extra = "allow"
        allow_mutation = False
        arbitrary_types_allowed = True

    _dependencies: List[str] = PrivateAttr()
    _uuid: str = PrivateAttr()

    name: Optional[
        str
    ]  # name of the `ProtocolUnit` that produced this `ProtocolUnitResult`
    data: Any  # should likely be fleshed out, currently a free-for-all

    def __init__(
        self, dependencies: Optional[Iterable["ProtocolUnitResult"]] = None,
            **data
    ):
        """
        Parameters
        ----------
        dependencies : Optional[Iterable["ProtocolUnitResult"]]
          The `ProtocolUnitResult`s from the predecessors of the `ProtocolUnit`
          that creates this `ProtocolUnitResult`.
        **data
          All other keyword arguments are retained as attributes of this
          `ProtocolUnitResult`.
        """

        if dependencies is None:
            dependencies = []

        super().__init__(**data)

        self._uuid = str(uuid.uuid4())

        dep_uuids = [dep._uuid for dep in dependencies]
        self._dependencies = dep_uuids

    def ok(self) -> bool:
        return True


class ProtocolUnitFailure(ProtocolUnitResult):
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

