# This code is part of OpenFE and is licensed under the MIT license.
# For details, see https://github.com/OpenFreeEnergy/openfe

import itertools
from typing import Optional, Iterable, List, Dict, Any
from collections import defaultdict

import pytest
import networkx as nx
import numpy as np

from gufe.chemicalsystem import ChemicalSystem
from gufe.mapping import ComponentMapping
from gufe.protocols import (
    Protocol,
    ProtocolDAG,
    ProtocolUnit,
    ProtocolResult,
    ProtocolDAGResult,
    ProtocolUnitResult,
    ProtocolUnitFailure,
)

from gufe.protocols.protocoldag import execute

from .test_tokenization import GufeTokenizableTestsMixin


class InitializeUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx, *, settings, stateA, stateB, mapping, start, **inputs):
        return dict(
            log="initialized",
        )


class SimulationUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx, *, initialization, **inputs):
        output = [initialization.outputs['log']]
        output.append("running_md_{}".format(inputs["window"]))

        return dict(
            log=output,
            window=inputs["window"],
            key_result=(100 - (inputs["window"] - 10)**2),
            scratch=ctx.scratch,
            shared=ctx.shared
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
    def _default_settings(cls):
        return {}

    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[dict[str, ComponentMapping]] = None,
        extend_from: Optional[ProtocolDAGResult] = None,
    ) -> List[ProtocolUnit]:

        # rip apart `extend_from` if needed to feed into `InitializeUnit`
        if extend_from is not None:
            # this is an example; wouldn't want to pass in whole ProtocolDAGResult into
            # any ProtocolUnits below, since this could create dependency hell;
            # instead, extract what's needed from it for starting point here
            starting_point = extend_from.protocol_unit_results[-1].outputs['final_positions']
        else:
            starting_point = None

        # convert protocol inputs into starting points for independent simulations
        alpha = InitializeUnit(
            name="the beginning",
            settings=self.settings,
            stateA=stateA,
            stateB=stateB,
            mapping=mapping,
            start=starting_point,
            some_dict={'a': 2, 'b': 12})

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

        return dict(outputs)


class BrokenSimulationUnit(SimulationUnit):
    @staticmethod
    def _execute(ctx, **inputs):
        raise ValueError("I have failed my mission", {'data': 'lol'})


class BrokenProtocol(DummyProtocol):
    def _create(
        self,
        stateA: ChemicalSystem,
        stateB: ChemicalSystem,
        mapping: Optional[dict[str, ComponentMapping]] = None,
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


class TestProtocol(GufeTokenizableTestsMixin):

    cls = DummyProtocol
    key = "DummyProtocol-5660965464c9afdaac0ac4486a9566b3"
    
    @pytest.fixture
    def instance(self):
        return DummyProtocol(settings=None)

    @pytest.fixture
    def protocol_dag(self, solvated_ligand, vacuum_ligand):
        protocol = DummyProtocol(settings=None)
        dag = protocol.create(
            stateA=solvated_ligand, stateB=vacuum_ligand, name="a dummy run"
        )
        dagresult: ProtocolDAGResult = execute(dag)

        return protocol, dag, dagresult

    @pytest.fixture
    def protocol_dag_broken(self, solvated_ligand, vacuum_ligand):
        protocol = BrokenProtocol(settings=None)
        dag = protocol.create(
            stateA=solvated_ligand, stateB=vacuum_ligand, name="a broken dummy run"
        )

        dagfailure: ProtocolDAGResult = execute(dag)

        return protocol, dag, dagfailure

    def test_dag_execute(self, protocol_dag):
        protocol, dag, dagresult = protocol_dag

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
        
        # check that shared directory the same for all simulations
        assert len(set(i.outputs['shared'] for i in simulationresults)) == 1

        # check that scratch directory is different for all simulations
        assert len(set(i.outputs['scratch'] for i in simulationresults)) == len(simulationresults)

    def test_dag_execute_failure(self, protocol_dag_broken):
        protocol, dag, dagfailure = protocol_dag_broken

        assert not dagfailure.ok()
        assert isinstance(dagfailure, ProtocolDAGResult)

        failed_units = dagfailure.protocol_unit_failures

        assert len(failed_units) == 1
        assert failed_units[0].name == "problem child"

        # parse exception arguments
        assert failed_units[0].exception[1][1]['data'] == "lol"
        assert isinstance(failed_units[0], ProtocolUnitFailure)

        succeeded_units = dagfailure.protocol_unit_results

        assert len(succeeded_units) > 0

    def test_create_execute_gather(self, protocol_dag):
        protocol, dag, dagresult = protocol_dag

        assert dagresult.ok()

        # gather aggregated results of interest
        protocolresult = protocol.gather([dagresult])

        assert len(protocolresult.data['logs']) == 1
        assert len(protocolresult.data['logs'][0]) == 21 + 1

        assert protocolresult.get_estimate() == 105336

    class ProtocolDAGTestsMixin(GufeTokenizableTestsMixin):
        
        def test_protocol_units(self, instance):
            # ensure that protocol units are given in-order based on DAG
            # dependencies
            checked = []
            for pu in instance.protocol_units:
                assert set(pu.dependencies).issubset(checked) 
                checked.append(pu)

        def test_graph(self, instance):
            assert isinstance(instance.graph, nx.DiGraph)

            # walk the nodes, check dependencies as given by each node against
            # edges in the graph
            for node in instance.graph.nodes:

                # check that each dep is represented by an edge
                for dep in node.dependencies:
                    assert isinstance(instance.graph.edges[node, dep], dict)

                # check that each edge corresponds to a known dependency
                for neighbor in instance.graph.neighbors(node):
                    assert neighbor in node.dependencies

        def test_key_stable(self, instance):
            # for the DAG system, keys for `ProtocolUnit`s are based on UUIDs,
            # so keys aren't stable up through `ProtocolDAG`s
            pass

    class TestProtocolDAG(ProtocolDAGTestsMixin):
        cls = ProtocolDAG
        key = "..."
        
        @pytest.fixture
        def instance(self, protocol_dag):
            protocol, dag, dagresult = protocol_dag
            return dag

    class TestProtocolDAGResult(ProtocolDAGTestsMixin):
        cls = ProtocolDAGResult
        key = "..."

        @pytest.fixture
        def instance(self, protocol_dag):
            protocol, dag, dagresult = protocol_dag
            assert dagresult.ok()
            return dagresult

        def test_protocol_unit_results(self, instance: ProtocolDAGResult):
            # ensure that protocolunitresults are given in-order based on DAG
            # dependencies
            checked: List[ProtocolUnitResult] = []
            for pur in instance.protocol_unit_results:
                assert set(pur.dependencies).issubset(checked)
                checked.append(pur)

        def test_result_graph(self, instance: ProtocolDAGResult):
            assert isinstance(instance.result_graph, nx.DiGraph)

            # walk the nodes, check dependencies as given by each node against
            # edges in the graph
            for node in instance.result_graph.nodes:

                # check that each dep is represented by an edge
                for dep in node.dependencies:
                    assert isinstance(instance.result_graph.edges[node, dep], dict)

                # check that each edge corresponds to a known dependency
                for neighbor in instance.result_graph.neighbors(node):
                    assert neighbor in node.dependencies

        def test_unit_to_result(self, instance: ProtocolDAGResult):
            # check that every unit has a result that we can retrieve
            for pu in instance.protocol_units:
                pur: ProtocolUnitResult = instance.unit_to_result(pu)
                assert pur.source_key == pu.key

        def test_result_to_unit(self, instance: ProtocolDAGResult):
            for pur in instance.protocol_unit_results:
                pu: ProtocolUnit = instance.result_to_unit(pur)
                assert pu.key == pur.source_key

        def test_protocol_unit_failures(self, instance: ProtocolDAGResult):
            assert len(instance.protocol_unit_failures) == 0

    class TestProtocolDAGResultFailure(ProtocolDAGTestsMixin):
        cls = ProtocolDAGResult
        key = "..."

        @pytest.fixture
        def instance(self, protocol_dag_broken):
            protocol, dag, dagfailure = protocol_dag_broken
            assert not dagfailure.ok()
            return dagfailure

        def test_protocol_unit_failures(self, instance: ProtocolDAGResult):
            # protocolunitfailures should have no dependents
            for puf in instance.protocol_unit_failures:

                assert all([puf not in pu.dependencies
                            for pu in instance.protocol_unit_results])

            for node in instance.result_graph.nodes:
                with pytest.raises(KeyError):
                    instance.result_graph.edges[node, puf]

        def test_protocol_unit_failure_exception(self, instance: ProtocolDAGResult):
            for puf in instance.protocol_unit_failures:
                isinstance(puf.exception, ValueError)

        def test_protocol_unit_failure_traceback(self, instance: ProtocolDAGResult):
            for puf in instance.protocol_unit_failures:
                assert "I have failed my mission" in puf.traceback

    class TestProtocolUnit(GufeTokenizableTestsMixin):
        cls = SimulationUnit
        key = "..."
    
        @pytest.fixture
        def instance(self, vacuum_ligand, solvated_ligand):
    
            # convert protocol inputs into starting points for independent simulations
            alpha = InitializeUnit(
                name="the beginning",
                settings={},
                stateA=vacuum_ligand,
                stateB=solvated_ligand,
                mapping=None,
                start=None,
                some_dict={'a': 2, 'b': 12},
            )
    
            return SimulationUnit(name=f"simulation", initialization=alpha)

        def test_key_stable(self, instance):
            # for the DAG system, keys for `ProtocolUnit`s are based on UUIDs,
            # so keys aren't stable up through `ProtocolDAG`s
            pass

class NoDepUnit(ProtocolUnit):
    @staticmethod
    def _execute(ctx, **inputs) -> Dict[str, Any]:
        return {'local': inputs['val'] ** 2}


class NoDepResults(ProtocolResult):
    def get_estimate(self):
        return sum(self.data['vals'])

    def get_uncertainty(self):
        return len(self.data['vals'])

    def get_rate_of_convergence(self):
        return 0.0


class NoDepsProtocol(Protocol):
    """A protocol without dependencies"""
    result_cls = NoDepResults

    @classmethod
    def _defaults(cls):
        return {}

    @classmethod
    def _default_settings(cls):
        return {}

    def _create(
            self,
            stateA: ChemicalSystem,
            stateB: ChemicalSystem,
            mapping: Optional[dict[str, ComponentMapping]] = None,
            extend_from: Optional[ProtocolDAGResult] = None,
    ) -> List[ProtocolUnit]:
        return [NoDepUnit(settings=self.settings,
                          val=i)
                for i in range(3)]

    def _gather(self, dag_results):
        return {
            'vals': list(itertools.chain.from_iterable(
                (d.outputs['local'] for d in dag.protocol_unit_results) for dag in dag_results)),
        }


class TestNoDepProtocol:
    def test_create(self):
        p = NoDepsProtocol()

        dag = p.create(None, None)

        assert len(dag.protocol_units) == 3

    def test_gather(self):
        p = NoDepsProtocol()

        dag = p.create(None, None)

        dag_result = execute(dag)

        assert dag_result.ok()

        result = p.gather([dag_result])

        assert result.get_estimate() == 0 + 1 + 4
        assert result.get_uncertainty() == 3
