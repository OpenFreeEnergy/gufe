# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

from typing import Optional, Iterable, List, Dict, Any

import pytest
import networkx as nx

from gufe.chemicalsystem import ChemicalSystem
from gufe.mapping import Mapping
from gufe.protocols import (
    Protocol,
    ProtocolDAG,
    ProtocolUnit,
    ProtocolResult,
    ProtocolDAGResult,
    ProtocolUnitResult,
)


class InitializeUnit(ProtocolUnit):
    def _execute(self, dependency_results):

        return dict(
            data="initialized",
        )


class SimulationUnit(ProtocolUnit):
    def _execute(self, dependency_results):

        output = [r.data for r in dependency_results]
        output.append("running_md_{}".format(self._kwargs["window"]))

        return dict(
            data=output,
            window=self._kwargs["window"],  # extra attributes allowed
        )


class FinishUnit(ProtocolUnit):
    def _execute(self, dependency_results):

        output = [r.data for r in dependency_results]
        output.append("assembling_results")

        return dict(
            data=output,
        )


class DummyProtocolResult(ProtocolResult):
    def get_estimate(self):
        ...

    def get_uncertainty(self):
        ...

    def get_rate_of_convergence(self):
        ...


class DummyProtocol(Protocol):

    _results_cls = DummyProtocolResult

    @classmethod
    def get_default_settings(cls):
        return {}

    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[Mapping] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
    ) -> nx.DiGraph:

        alpha = InitializeUnit(
            self.settings,
            stateA=stateA,
            stateB=stateB,
            mapping=mapping,
            start=extend_from,
        )

        simulations: List[ProtocolUnit] = [
            SimulationUnit(self.settings, window=i) for i in range(20)
        ]

        omega = FinishUnit(self.settings, name="the end")

        dag = nx.DiGraph()
        for sim in simulations:
            dag.add_edge(sim, alpha)
            dag.add_edge(omega, sim)

        return dag

    def _gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> Dict[str, Any]:

        outputs = []
        for pdr in protocol_dag_results:
            for pur in pdr.protocol_unit_results:
                if pur.name == "the end":
                    outputs.append(pur.data)

        return dict(data=outputs)


class TestProtocol:
    def test_init(self):
        protocol = DummyProtocol(settings=None)

    def test_create(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings=None)
        dag = protocol.create(stateA=solvated_ligand, stateB=solvated_complex)

    def test_create_execute_gather(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings=None)
        dag = protocol.create(
            stateA=solvated_ligand, stateB=solvated_complex, name="a dummy run"
        )
        dagresult = dag.execute()

        finishresult = [
            dagresult.graph.nodes[u]["result"]
            for u in dagresult.graph.nodes
            if u.__class__.__name__ == "FinishUnit"
        ][0]
        simulationresults = [
            dagresult.graph.nodes[u]["result"]
            for u in dagresult.graph.nodes
            if u.__class__.__name__ == "SimulationUnit"
        ]

        # check that we have dependency information in results
        assert set(finishresult._dependencies) == {u._uuid for u in simulationresults}

        # check that we have as many units as we expect in resulting graph
        assert len(dagresult.graph) == 22

        # gather aggregated results of interest
        protocolresult = protocol.gather([dagresult])

        assert len(protocolresult.data) == 1
        assert len(protocolresult.data[0]) == 20 + 1
