# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

from typing import Optional, Iterable

import pytest

from gufe.protocols import (Protocol, ProtocolDAG, ProtocolUnit, ProtocolResult,
 ProtocolDAGResult, ProtocolUnitResult)


class InitializeUnit(ProtocolUnit):

    def _execute(self, dependency_results):

        return dict(
                data="initialized",
                )


class SimulationUnit(ProtocolUnit):

    def _execute(self, dependency_results):

        output = [r.data for r in dependency_results]
        output.append("running_md_{}".format(self._kwargs['window']))

        return dict(
                data=output,
                window=self._kwargs['window'], # extra attributes allowed
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

    def _create(self,
            initial: "ChemicalSystem", 
            final: "ChemicalSystem",
            mapping: Optional["Mapping"] = None,
            extend_from: Optional[ProtocolDAGResult] = None,
            settings: Optional["ProtocolSettings"] = None
        ) -> "ProtocolDAG":

        if settings is None:
            settings = self.settings

        alpha = InitializeUnit(
                self.settings, initial=initial, final=final, mapping=mapping, start=extend_from)

        simulations = [SimulationUnit(self.settings, alpha, window=i)  for i in range(20)]

        omega = FinishUnit(self.settings, *simulations)

        return ProtocolDAG([alpha, omega] + simulations)

    def _gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]) -> ProtocolResult:

        outputs = []
        for pdr in protocol_dag_results:
            for pu in pdr.units:
                if pu.name == "FinishUnit":
                    outputs.append(pu.data)

        return dict(data=outputs)


class TestProtocol:

    def test_init(self):
        protocol = DummyProtocol(settings={})

    def test_create(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings={})
        dag = protocol.create(initial=solvated_ligand, final=solvated_complex)

    def test_create_execute_gather(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings={})
        dag = protocol.create(initial=solvated_ligand, final=solvated_complex)
        dagresult = dag.execute()

        assert len(dagresult.units) == 22

        protocolresult = protocol.gather([dagresult])

        assert len(protocolresult.data) == 1
        assert len(protocolresult.data[0]) == 20 + 1
