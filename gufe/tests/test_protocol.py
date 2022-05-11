# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

from typing import Optional, Iterable

import pytest

from gufe.protocols import Protocol, ProtocolDAG, ProtocolUnit
from gufe.protocols.results import ProtocolResult, ProtocolDAGResult, ProtocolUnitResult


class InitializeUnit(ProtocolUnit):

    def _execute(self, dependency_results):

        self._result = ProtocolUnitResult(output="initialized")
        return self._result

    def _results(self):
        return self._result


class SimulationUnit(ProtocolUnit):

    def _execute(self, dependency_results):

        import pdb
        pdb.set_trace()

        self._output = "\n".join([r.output for r in dependency_results])
        self._output += "\nrunning_md"

        return self._output

    def _results(self):
        return ProtocolUnitResult(output=self._output)


class FinishUnit(ProtocolUnit):

    def _execute(self, dependency_results):

        self._output = "\n".join([r.output for r in dependency_results])
        self._output += "\nassembling_results"

        return self._output

    def _results(self):
        return ProtocolUnitResult(output=self._output)


class DummyProtocol(Protocol):

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

    def _gather(self, protocol_dag_results: Iterable[ProtocolDAGResult]):
        ...
        


class TestProtocol:

    def test_init(self):
        protocol = DummyProtocol(settings={})

    def test_create(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings={})
        dag = protocol.create(initial=solvated_ligand, final=solvated_complex)

    def test_create_execute(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings={})
        dag = protocol.create(initial=solvated_ligand, final=solvated_complex)
        dagresult = dag.execute()

        protocol.gather(dagresult)
