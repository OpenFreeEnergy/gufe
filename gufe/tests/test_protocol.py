# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

from typing import Optional, Iterable, List, Dict, Any
from collections import defaultdict

import pytest
import networkx as nx
import numpy as np

from gufe.chemicalsystem import ChemicalSystem
from gufe.mapping import Mapping
from gufe.protocols import (
    Protocol,
    ProtocolDAG,
    ProtocolUnit,
    ProtocolResult,
    ProtocolDAGResult,
    ProtocolDAGFailure,
    ProtocolUnitResult,
    ProtocolUnitFailure,
)

from .test_tokenization import GufeTokenizableTestsMixin


class InitializeUnit(ProtocolUnit):

    @staticmethod
    def _execute(ctx, *, settings, stateA, stateB, mapping, start, **inputs):
        return dict(
            log="initialized",
        )


class SimulationUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx, *, initialization: ProtocolUnitResult, **inputs):

        output = [initialization.outputs['log']]
        output.append("running_md_{}".format(inputs["window"]))

        return dict(
            log=output,
            window=inputs["window"],
            key_result=(100 - (inputs["window"] - 10)**2),
            unit_scratch=ctx.unit_scratch,
            dag_scratch=ctx.dag_scratch
        )


class FinishUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx, *, simulations, **inputs):

        output = [s.outputs['log'] for s in simulations]
        output.append("assembling_results")

        key_results = {s.inputs['window']: s.outputs['key_result'] for s in simulations}

        return dict(
            log=output,
            key_results=key_results
        )


class DummyProtocolResult(ProtocolResult):
    def get_estimate(self):
        # we'll pretend that the free energy estimate here is the sum of the
        # product of neighboring simulation window `key_result`s

        dgs = []
        for sample in self.data['key_results']:
            windows = sorted(sample.keys())
            dg = 0
            for i, j in zip(windows[:-1], windows[1:]):
                dg += sample[i] * sample[j]

            dgs.append(dg)

        return np.mean(dg)

    def get_uncertainty(self):
        ...

    def get_rate_of_convergence(self):
        ...


class DummyProtocol(Protocol):

    result_cls = DummyProtocolResult

    @classmethod
    def get_default_settings(cls):
        return {}

    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[Mapping] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
    ) -> List[ProtocolUnit]:

        # convert protocol inputs into starting points for independent simulations
        alpha = InitializeUnit(
            name="the beginning",
            settings=self.settings,
            stateA=stateA,
            stateB=stateB,
            mapping=mapping,
            start=extend_from,
            some_dict={'a': 2, 'b': 12},
            pure=True,          # states that this unit is purely a function of its explicit arguments;
                                # allows scheduler to deduplicate based on inputs
        )

        # create several units that would each run an independent simulation
        simulations: List[ProtocolUnit] = [
            SimulationUnit(settings=self.settings, name=f"sim {i}", window=i, initialization=alpha) for i in range(21)
        ]

        # gather results from simulations, finalize outputs
        omega = FinishUnit(settings=self.settings, name="the end", simulations=simulations)

        # return all `ProtocolUnit`s we created
        return [alpha, *simulations, omega]

    def _gather(
        self, protocol_dag_results: Iterable[ProtocolDAGResult]
    ) -> Dict[str, Any]:

        outputs = defaultdict(list)
        for pdr in protocol_dag_results:
            for pur in pdr.protocol_unit_results:
                if pur.name == "the end":
                    outputs['logs'].append(pur.outputs['log'])
                    outputs['key_results'].append(pur.outputs['key_results'])

        return dict(data=outputs)


class BrokenSimulationUnit(SimulationUnit):
    @staticmethod
    def _execute(ctx, **inputs):
        raise ValueError("I have failed my mission", {'data': 'lol'})


class BrokenProtocol(DummyProtocol):
    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[Mapping] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
    ) -> nx.DiGraph:

        # convert protocol inputs into starting points for independent simulations
        alpha = InitializeUnit(
            settings=self.settings,
            stateA=stateA,
            stateB=stateB,
            mapping=mapping,
            start=extend_from,
        )

        # create several units that would each run an independent simulation
        simulations: List[ProtocolUnit] = [
            SimulationUnit(settings=self.settings, name=f"sim {i}", window=i, initialization=alpha) for i in range(21)
        ]

        # introduce a broken ProtocolUnit
        simulations.append(BrokenSimulationUnit(settings=self.settings, window=21, name="problem child", initialization=alpha))

        # gather results from simulations, finalize outputs
        omega = FinishUnit(settings=self.settings, name="the end", simulations=simulations)

        # return all `ProtocolUnit`s we created
        return [alpha, *simulations, omega]


class TestProtocol:
    def test_init(self):
        protocol = DummyProtocol(settings=None)

    def test_create(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings=None)
        dag = protocol.create(stateA=solvated_ligand, stateB=solvated_complex)

    def test_create_execute(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings=None)
        dag = protocol.create(
            stateA=solvated_ligand, stateB=solvated_complex, name="a dummy run"
        )
        dagresult: ProtocolDAGResult = dag.execute()

        assert dagresult.ok()

        # the FinishUnit will always be the last to execute
        finishresult = dagresult.protocol_unit_results[-1]
        assert finishresult.name == "the end"

        # gather SimulationUnits
        simulationresults = [dagresult.unit_to_result(pu)
                             for pu in dagresult.protocol_units
                             if isinstance(pu, SimulationUnit)]

        # check that we have dependency information in results
        assert set(finishresult.inputs['simulations']) == {u for u in simulationresults}

        # check that we have as many units as we expect in resulting graph
        assert len(dagresult.graph) == 23
        
        # check that dag_scratch directory the same for all simulations
        assert len(set(i.outputs['dag_scratch'] for i in simulationresults)) == 1

        # check that unit_scratch directory is different for all simulations
        assert len(set(i.outputs['unit_scratch'] for i in simulationresults)) == len(simulationresults)

    def test_create_execute_failure(self, solvated_ligand, solvated_complex):
        protocol = BrokenProtocol(settings=None)

        dag = protocol.create(
            stateA=solvated_ligand, stateB=solvated_complex, name="a broken dummy run"
        )

        dagfailure: ProtocolDAGFailure = dag.execute()

        assert not dagfailure.ok()
        assert isinstance(dagfailure, ProtocolDAGFailure)

        failed_units = dagfailure.protocol_unit_failures

        assert len(failed_units) == 1
        assert failed_units[0].name == "problem child"
        assert failed_units[0].exception.args[1]['data'] == "lol"
        assert isinstance(failed_units[0], ProtocolUnitFailure)

        succeeded_units = dagfailure.protocol_unit_results

        assert len(succeeded_units) > 0

    def test_create_execute_gather(self, solvated_ligand, solvated_complex):
        protocol = DummyProtocol(settings=None)
        dag = protocol.create(
            stateA=solvated_ligand, stateB=solvated_complex, name="a dummy run"
        )
        dagresult: ProtocolDAGResult = dag.execute()

        assert dagresult.ok()

        # gather aggregated results of interest
        protocolresult = protocol.gather([dagresult])

        assert len(protocolresult.data['logs']) == 1
        assert len(protocolresult.data['logs'][0]) == 21 + 1

        assert protocolresult.get_estimate() == 105336


class TestProtocolUnit(GufeTokenizableTestsMixin):

    cls = SimulationUnit

    @pytest.fixture
    def instance(self, solvated_complex, solvated_ligand):

        # convert protocol inputs into starting points for independent simulations
        alpha = InitializeUnit(
            name="the beginning",
            settings={},
            stateA=solvated_complex,
            stateB=solvated_ligand,
            mapping=None,
            start=None,
            some_dict={'a': 2, 'b': 12},
            pure=True,
        )

        return SimulationUnit(name=f"simulation", initialization=alpha, pure=True)

    def test_purity_behavior(self, instance):
        ser = instance.to_dict()
        deser = self.cls.from_dict(ser)

        # instance and deser should be both equal and identical
        deser == instance
        deser is instance

        # if we make two non-pure versions, the keys for these won't be a
        # function of their intputs, but a uuid
        ser['pure'] = False
        deser_impure_1 = self.cls.from_dict(ser)
        deser_impure_2 = self.cls.from_dict(ser)

        # they are equal in contents
        assert deser_impure_1 == deser_impure_2

        # but they are not the same object
        assert deser_impure_1 is not deser_impure_2


